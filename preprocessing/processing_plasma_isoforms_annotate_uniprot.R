setwd("/Users/lab/Documents/githubspot/plasma_proteoforms/")
#Downgrade dbplyr
#devtools::install_version("dbplyr", version = "2.3.4")
library(biomaRt)
library(curl)
library(Rcpi)
library(Biostrings)
library(purrr)
library(dplyr)
library(WebGestaltR)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)

#peptides <- read.delim("data/peptide-olink_cors.csv", sep = ",", row.names = "X", header = T)
#peptides$cluster <- readRDS("leiden_1_partition_peptides_bypeptide.RDS")
#peptides <- readRDS(file =  "~/Documents/githubspot/plasma_proteoforms/data/peptides.RDS")

#or skip the next part if you already ran it and just want to make plots
#plasma_data <- readRDS("output/plasma_data_FASTA.RDS")

plasma_data <- read.delim("data/LCP114_PSM_table_filtered.txt", header = T)
plasma_data <- unique(plasma_data[,colnames(plasma_data) %in% c("Peptide", "Gene.Name", "Master.protein", "Peptide.label")])
#peptides$Peptide <- plasma_data$Peptide
partition<- readRDS("leiden_1_partition_Genesymb_peptides_bypeptide.RDS")
partition$'Peptide.label' <- partition$peptide
plasma_data <- left_join(plasma_data, partition, by = "Peptide.label")
####################################################
####################################################
#mart <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')
# uniprot_ids <- getBM(attributes = c("hgnc_symbol",#"pdb",
#                                     #"uniprot_gn_symbol",
#                                     "uniprotswissprot"
# ),
# values = gsub("_[0-9]", "", rownames(peptides)),
# mart = mart)

#filter which did not get a match
#uniprot_ids <- uniprot_ids[uniprot_ids$uniprotswissprot != "",]

#uniprot_ids <- data.frame("hgnc_symbol" = plasma_data$Gene.Name, "uniprotswissprot" = sub("^sp\\|([^|]+)\\|.*", "\\1", plasma_data$Master.protein))
#get fasta sequences for each hit

clean_substrings <- function(text) {
  parts <- strsplit(text, ";")[[1]]  # Split the string by ";"
  cleaned_parts <- sapply(parts, function(part) gsub("\\|.*", "", part))  # Clean each part
  return(paste(cleaned_parts, collapse = ";"))  # Rejoin parts with ";"
}

plasma_data$uniprotswissprot <- gsub("sp\\|", "", plasma_data$Master.protein)

plasma_data$uniprotswissprot <- sapply(plasma_data$uniprotswissprot, clean_substrings)

#filter to one uniprot map matches only
plasma_data <- plasma_data[!grepl(";", plasma_data$uniprotswissprot),]
# 
# get_FASTA <- function(uniprot_iois, just_first = TRUE){
#   uniprot_iois$fasta <- 0
#   for (i in 1:length(uniprot_iois$uniprotswissprot)) {
#     uniprot_ioi <- uniprot_iois$uniprotswissprot[i]
#     
#     #add the file type
#     FASTA_path <- paste0(tempfile(), ".fasta")
#     FASTA_file <- curl_fetch_disk(url = paste0('https://rest.uniprot.org/uniprotkb/', uniprot_ioi, '.fasta'),
#                                   handle = new_handle(),
#                                   path = FASTA_path)
#     
#     AASeq <- readFASTA(file = FASTA_file[["content"]])
#     AASeq <- flatten(AASeq)
#     
#     uniprot_iois$fasta[i] <- AASeq
#     
#     unlink(FASTA_file)
#     
#     # while(request_num <= uniprot_ids_options & alphafold_file[["status_code"]] == 404){
#     #   request_num <-  request_num +1
#     #   uniprot_ioi <- uniprot_ids[uniprot_ids$hgnc_symbol == ioi,]$uniprot_gn_id[request_num]
#     #   alphafold_file <- curl_fetch_disk(url = paste0("https://alphafold.ebi.ac.uk/files/AF-",uniprot_ioi,"-F1-model_v4.pdb"),
#     #                                     handle = new_handle(),
#     #                                     path = alphafold_path)
#     # }
#   }
#   
#   # if(just_first == TRUE){
#   #   
#   #   uniprot_iois <- unique(uniprot_iois)
#   #   
#   # }
#   
#   
#   return(uniprot_iois)
# }

get_FASTA <- function(uniprot_ids, just_first = TRUE) {
  results <- vector("list", length(uniprot_ids))
  names(results) <- uniprot_ids
  
  for (i in seq_along(uniprot_ids)) {
    uniprot_id <- uniprot_ids[i]
    
    # Prepare the URL and temporary file path
    url <- paste0('https://rest.uniprot.org/uniprotkb/', uniprot_id, '.fasta')
    fasta_path <- tempfile(fileext = ".fasta")
    
    # Download FASTA file
    request <- curl_fetch_disk(url, path = fasta_path)
    
    if (file.info(fasta_path)$size > 0 && grepl("^>", readLines(fasta_path, n = 1))) {
      # Read and process the FASTA file if it's valid
      fasta_data <- readDNAStringSet(fasta_path)
      results[[i]] <- toString(fasta_data)
    } else {
      results[[i]] <- NA  # Handle cases where no valid FASTA data is retrieved
    }
    
    # Clean up the temporary file
    unlink(fasta_path)
  }
  
  # Optionally, reduce to unique results if just_first is TRUE
  if (just_first) {
    results <- unique(results)
  }
  
  return(results)
}

fasta <- get_FASTA(unique(plasma_data$uniprotswissprot))
names(fasta) <- unique(plasma_data$uniprotswissprot)


fasta_df <- enframe(fasta, name = "uniprotswissprot", value = "fasta")
plasma_data <- left_join(plasma_data, unique(fasta_df), by = "uniprotswissprot")
#need the col name to line up
plasma_data$hgnc_symbol <- plasma_data$Gene.Name
#save out
saveRDS(fasta, "output/FASTA_list_obj.RDS")
saveRDS(plasma_data, "output/plasma_data_FASTA.RDS")
plasma_data <- readRDS("output/plasma_data_FASTA.RDS")

mart <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')  
View(listAttributes(mart))
domains <- getBM(attributes = c(#"hgnc_symbol",#"pdb",
                                "uniprotswissprot",
                                #"interpro_short_description"
                                "interpro_description"
                                #"biogrid"
),
values =unique(plasma_data$uniprotswissprot),
mart = mart)

domains <- getBM(attributes = c(#"hgnc_symbol",#"pdb",
                                #"uniprot_gn_symbol",
                                 "uniprotswissprot",
                                #"interpro_short_description",
                                "interpro_description",
                                "interpro_start",
                                "interpro_end"
                                #"biogrid"
),
values =distinct(domains),
mart = mart)

# biogrid <- getBM(attributes = c("hgnc_symbol",#"pdb",
#                                 #"uniprot_gn_symbol",
#                                 #"interpro"
#                                 "biogrid"
# ),
# values = unique(plasma_data$Gene.Name),
# mart = mart)

#remove missing
domains <- domains[complete.cases(domains) &
                     trimws(domains$uniprotswissprot) != "",]

saveRDS(domains, "output/interpro_domains.RDS")
#saveRDS(biogrid, "output/biogrid_interactions.RDS")

# #Make peptide seq list
# peptide_seq_list <- list()
# for(i in unique(plasma_data$Gene.Name)){
#   
#   peptide_seq_list[[i]] <- data.frame(
#     
#     as.data.frame(gsub("[0-9\\.\\+]", "",
#                        plasma_data[plasma_data$Gene.Name == i,]$Peptide)),
#     #as.data.frame(vertex_attr(peptide_mappings[[i]])[["name"]]),
#     as.data.frame(plasma_data[plasma_data$Gene.Name == i,]$cluster))
#   colnames(peptide_seq_list[[i]]) <- c("peptides", "id", #"peptides_full", 
#                                        "membership" )
#   
#   
# } #end for
# saveRDS(peptide_seq_list, "output/peptide_seq_list.RDS")



















