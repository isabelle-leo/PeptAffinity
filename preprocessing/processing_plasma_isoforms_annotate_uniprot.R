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
library(tidyverse)
library(ComplexHeatmap)
library(jsonlite)
library(stringr)



plasma_data <- readRDS("~/Documents/GitHub/PeptOlink/PeptOlink/data/plasma_data_FASTA_isoforms.RDS")
isoforms <- unique(unlist(strsplit(plasma_data$UniProt.MS, ";")))
mart <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')  
View(listAttributes(mart))
# domains <- getBM(attributes = c(#"hgnc_symbol",#"pdb",
#                                 "uniprotswissprot",
#                                 #"interpro_short_description"
#                                 "interpro_description"
#                                 #"biogrid"
# ),
# values =isoforms,
# mart = mart)

domains <- getBM(attributes = c(#"hgnc_symbol",#"pdb",
                                #"uniprot_gn_symbol",
                                 "uniprotswissprot",
                                #"interpro_short_description",
                                "interpro_description",
                                "interpro_start",
                                "interpro_end"
                                #"biogrid"
),
values =unique(isoforms),
mart = mart)

domains_isoform <- getBM(attributes = c(#"hgnc_symbol",#"pdb",
  #"uniprot_gn_symbol",
  "uniprot_isoform",
  #"interpro_short_description",
  "interpro_description",
  "interpro_start",
  "interpro_end"
  #"biogrid"
),
values =unique(isoforms),
mart = mart)

colnames(domains_isoform) <- c("uniprotswissprot",
                               "interpro_description",
                               "interpro_start",
                               "interpro_end")
domains <- rbind(domains, domains_isoform)
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
#standardize names
colnames(domains) <- c("uniprot",
  "description",
  "start",
  "end")


saveRDS(domains, "output/interpro_domains.RDS")


prosite <- getBM(attributes = c(#"hgnc_symbol",#"pdb",
  #"uniprot_gn_symbol",
  "uniprotswissprot",
  #"interpro_short_description",
  "scanprosite",
  "scanprosite_start",
  "scanprosite_end"
  #"biogrid"
),
values =unique(isoforms),
mart = mart)

prosite_isoform <- getBM(attributes = c(#"hgnc_symbol",#"pdb",
  #"uniprot_gn_symbol",
  "uniprot_isoform",
  #"interpro_short_description",
  "scanprosite",
  "scanprosite_start",
  "scanprosite_end"
  #"biogrid"
),
values =unique(isoforms),
mart = mart)

colnames(prosite_isoform) <- c(  "uniprotswissprot",
                                 #"interpro_short_description",
                                 "scanprosite",
                                 "scanprosite_start",
                                 "scanprosite_end")
prosite<- rbind(prosite, prosite_isoform)

#remove missing
prosite <- prosite[complete.cases(prosite) &
                     trimws(prosite$uniprotswissprot) != "",]

#standardize names
colnames(prosite) <- c("uniprot",
                       "scanprosite",
                       "start",
                       "end")

#replace description using prosite info
doc_lines <- readLines("~/Downloads/prosite.doc")

ps_entries <- list()
current_ps <- NULL
title_next <- FALSE
title_line <- NULL

for (i in seq_along(doc_lines)) {
  line <- doc_lines[i]
  
  if (str_detect(line, "^\\{PS\\d{5};")) {
    # Save previous entry
    if (!is.null(current_ps)) {
      ps_entries[[current_ps$id]] <- list(
        id = current_ps$id,
        name = current_ps$name,
        title = title_line
      )
      title_line <- NULL
    }
    # Start new entry
    parts <- str_match(line, "\\{(PS\\d{5});\\s*(.+)\\}")
    current_ps <- list(id = parts[2], name = parts[3])
    title_next <- FALSE
    
  } else if (str_detect(line, "^\\*{3,}")) {
    # When first line of asterisks found, look for title in next line
    title_next <- TRUE
    
  } else if (title_next && str_detect(line, "^\\*.*\\*$")) {
    # Actual title is this line between asterisks
    title_line <- str_trim(str_remove_all(line, "^\\*+|\\*+$"))
    title_next <- FALSE
  }
}

# Save final entry
if (!is.null(current_ps)) {
  ps_entries[[current_ps$id]] <- list(
    id = current_ps$id,
    name = current_ps$name,
    title = title_line
  )
}

prosite$description <- sapply(
  prosite$scanprosite,
  function(id) {
    entry <- ps_entries[[id]]
    entry$title %||% entry$name
  }
)


saveRDS(prosite, "output/prosite_motifs.RDS")




