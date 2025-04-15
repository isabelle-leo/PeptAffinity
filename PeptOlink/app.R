library(shiny)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(curl)
library(Biostrings)
library(Rcpi)
library(shinyjs)
library(classInt)
library(ggpubr)
library(circlize)
library(shinyWidgets)
library(shinycssloaders)
library(shinythemes)
library(bslib)
library(shinyBS)
library(sysfonts)
library(showtext)
library(data.table)
library(viridis)
library(heatmaply)
library(purrr)
library(plotly)
library(RColorBrewer)
library(NGLVieweR)

#   _____ ______ _______ _    _ _____  
#  / ____|  ____|__   __| |  | |  __ \ 
# | (___ | |__     | |  | |  | | |__) |
#  \___ \|  __|    | |  | |  | |  ___/ 
#  ____) | |____   | |  | |__| | |     
# |_____/|______|  |_|   \____/|_|     

font_add_google("Open Sans", "open-sans")  
showtext_auto()  # Enable showtext globally

theme_minimal(base_family = "open-sans")

plasma_data <- readRDS("data/plasma_data_FASTA_isoforms.RDS")
plasma_dt <- as.data.table(plasma_data)
setkey(plasma_dt, Gene.Name)


color_breaks <- seq(-1, 1, length.out = 100)
#correlation_palette <- viridis(100, option = "D")  
#correlation_palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))(100)
# correlation_palette_function <- function(x) {
#   # Ensure x is within the allowed range
#   x <- pmax(pmin(x, 1), -1)
#   # For x below the first median, use the first categorical color
#   if(x <= -0.35) return("grey")
#   # For x above the last median, use the last categorical color
#   if(x >= 0.85) return("#206e8c")
#   
#   # For values in between, determine which segment to interpolate
#   if(x < 0.4) {
#     # Transition from 'No correlation' to 'Weak correlation'
#     t <- (x - -0.35) / (0.4 - -0.35)
#     col <- colorRamp(c("grey", "#9ccfe7"))(t)
#   } else if(x < 0.6) {
#     # Transition from 'Weak correlation' to 'Moderate correlation'
#     t <- (x - 0.4) / (0.6 - 0.4)
#     col <- colorRamp(c("#9ccfe7", "#629db8"))(t)
#   } else {
#     # Transition from 'Moderate correlation' to 'Strong correlation'
#     t <- (x - 0.6) / (0.85 - 0.6)
#     col <- colorRamp(c("#629db8", "#206e8c"))(t)
#   }
#   rgb(col[1], col[2], col[3], maxColorValue = 255)
# }
# correlation_palette <- sapply(color_breaks, correlation_palette_function)

correlation_palette_function <- function(x) {
  x <- pmax(pmin(x, 1), -1)
  if(x < 0) {
    pink_pal <- colorRampPalette(c("#e60045","#ff1a5e", "#ff5c8d", "#ff80a6", "#ffb3c9", "#ffe6ed"))(100)
    idx <- round((x - (-1)) / (0 - (-1)) * 99) + 1
    return(pink_pal[idx])
  } else if (0 <= x && x < .3) {
    neut_pal <- colorRampPalette(c("#ffe6ed", "gray95", viridis(100, option = "mako")[100]))(100)
    idx <- round((x - 0) / (0.3 - 0) * 99) + 1
    return(neut_pal[idx])
  }else {
    idx <- round(100 - (x - .3) / (1 - .3) * 99) + 1
    return(viridis(100, option = "mako")[idx])
  }
}

correlation_palette <- sapply(color_breaks, correlation_palette_function)

color_scale_medians <- c("No correlation" = -.1,  
             "Weak correlation" = 0.4, 
             "Moderate correlation" = 0.6,
             "Strong correlation" = 0.85)

categorical_colors <- c(
  "No correlation" = correlation_palette_function(color_scale_medians["No correlation"]),
  "Weak correlation" = correlation_palette_function(color_scale_medians["Weak correlation"]),
  "Moderate correlation" = correlation_palette_function(color_scale_medians["Moderate correlation"]),
  "Strong correlation" = correlation_palette_function(color_scale_medians["Strong correlation"])
)

# For summary scatter plots
mako_colors <- mako(5)

my_theme <- bs_theme(
  version = 4,
  bootswatch = "flatly",
  base_font = font_google("Open Sans"),
  heading_font = font_google("Open Sans"),
  primary = "#FF876F",
  fg = '#4F0433', 
  bg = "#fff",
  "font-size-base" = "0.9rem"
)

# Load the data
correlation_long_filt <- read.csv('data/peptide_cors_overlappingProteins_filt.csv') |>
  mutate(gene_symbol = Gene.Name.MS) |>
  dplyr::rename(correlation = Correlation)

#  ______ _    _ _   _  _____ _______ _____ ____  _   _  _____ 
# |  ____| |  | | \ | |/ ____|__   __|_   _/ __ \| \ | |/ ____|
# | |__  | |  | |  \| | |       | |    | || |  | |  \| | (___  
# |  __| | |  | | . ` | |       | |    | || |  | | . ` |\___ \ 
# | |    | |__| | |\  | |____   | |   _| || |__| | |\  |____) |
# |_|     \____/|_| \_|\_____|  |_|  |_____\____/|_| \_|_____/ 

get_color_vector <- function (colors,
                              vec,
                              gray = "#494e57") {
  
  n <- length(unique(vec))
  
  suppl_needed <- max(n - length(colors), 0)
  col_suppl <- c(colors[1:min(n, length(colors))], rep(gray, suppl_needed))
  
  # count appearances of each group
  vec_accum <- table(vec) %>% sort(decreasing = TRUE)
  
  names(col_suppl) <- names(vec_accum) %>% as.character()
  
  # add a gray missing value
  col_suppl <- c(col_suppl, c("0" = gray))
  
  return(col_suppl)
}

get_alphafold_file <- function(ioi){
  # Build the AlphaFold URL directly using the isoform
  alphafold_path <- paste0(tempfile(), ".pdb")
  url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", ioi, "-F1-model_v4.pdb")
  
  alphafold_file <- curl_fetch_disk(url = url,
                                    handle = new_handle(),
                                    path = alphafold_path)
  
  if (alphafold_file[["status_code"]] == 404) {
    stop("No AlphaFold structure was found for the provided isoform.")
  }
  
  return(alphafold_file)
}


get_FASTA_fromalphafold <- function(alphafold_file){
  # Extract the uniprot isoform id from the AlphaFold file URL
  uniprot_ioi <- gsub("https://alphafold.ebi.ac.uk/files/AF-", "", alphafold_file[["url"]])
  uniprot_ioi <- gsub("-F1-model_v4.pdb", "", uniprot_ioi)
  
  # Define temporary path and fetch the FASTA file from UniProt
  FASTA_path <- paste0(tempfile(), ".fasta")
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_ioi, ".fasta")
  FASTA_file <- curl_fetch_disk(url = url,
                                handle = new_handle(),
                                path = FASTA_path)
  
  if (FASTA_file[["status_code"]] != 200) {
    stop("Failed to fetch FASTA file from UniProt.")
  }
  
  # Read the FASTA sequence
  AASeq <- readFASTA(file = FASTA_path)
  fasta_seq <- toString(AASeq[[1]])
  
  # Clean up temporary file
  unlink(FASTA_path)
  
  return(fasta_seq)
}



get_peptide_and_correlation_numeric <- function(fasta_list, peptide_seq_list, exclude = TRUE){
  peptide_indices <- list()
  for( i in 1:length(peptide_seq_list$peptides)){
    peptide_indices[[rownames(peptide_seq_list[i,])]] <- data.frame(regexpr(peptide_seq_list$peptides[i], fasta_list, fixed=TRUE),
                                                                    peptide_seq_list$peptide[i],
                                                                    peptide_seq_list$N[i],
                                                                    peptide_seq_list$target_correlation[i])
    names(peptide_indices[[rownames(peptide_seq_list[i,])]]) <- c("start_position", "peptide", "N", "target_correlation")
    
    if(exclude == TRUE){
      
      if(peptide_indices[[rownames(peptide_seq_list[i,])]][1] == nchar(fasta_list)){
        
        peptide_indices[[rownames(peptide_seq_list[i,])]] <- NULL
        
      }
      
    }
  }
  
  return(peptide_indices)
}




# Function to compute WCSS and BCSS
compute_wcss_bcss <- function(data) {
  total_mean <- mean(data$correlation, na.rm = TRUE)
  
  data_summarized <- data %>%
    group_by(gene_symbol, cluster) %>%
    summarise(
      cluster_mean = mean(correlation, na.rm = TRUE),
      WCSS = sum((correlation - cluster_mean)^2, na.rm = TRUE),
      replicate_count = mean(replicate_count, na.rm = TRUE),
      .groups = 'drop'
    ) 
  
  data_summarized <- data_summarized %>%
    group_by(gene_symbol) %>%
    summarise(
      WCSS = sum(WCSS, na.rm = TRUE),
      BCSS = sum(replicate_count * (cluster_mean - total_mean)^2, na.rm = TRUE),
      .groups = 'drop'
    )
  
  return(data_summarized)
}




compute_jenks_clusters_simple <- function(data, value_column = "correlation", n_classes = 3) {
  data <- data %>%
    group_by(gene_symbol) %>%
    mutate(
      jenks_class = classIntervals(
        var = get(value_column),
        n = n_classes,
        style = "jenks"
      )$classif %>%
        cut(get(value_column), breaks = ., labels = c("Low", "Middle", "High"))
    ) %>%
    ungroup()
  return(data)
}



compute_jenks_clusters <- function(data, value_column = "correlation", n_classes = 3) {
  values <- data[[value_column]]
  
  # Ensure unique values for Jenks computation
  values <- unique(values)
  
  # Compute Jenks breaks
  breaks <- tryCatch(
    {
      if (length(values) > n_classes) {
        classIntervals(values, n = n_classes, style = "jenks")$brks
      } else {
        stop("Insufficient unique values for Jenks")
      }
    },
    error = function(e) {
      # Fallback to quantile-based breaks
      message("Fallback to quantile-based breaks")
      quantile(values, probs = seq(0, 1, length.out = n_classes + 1))
    }
  )
  
  # Ensure breaks are unique
  if (length(unique(breaks)) != length(breaks)) {
    breaks <- seq(min(values), max(values), length.out = n_classes + 1)
    message("Adjusted to evenly spaced breaks due to non-unique values")
  }
  
  # Generate labels dynamically based on the number of breaks
  labels <- c("Low", "Middle", "High")[1:(length(breaks) - 1)]
  
  # Assign classes based on the computed breaks
  data <- data %>%
    mutate(
      jenks_class = cut(
        .data[[value_column]],
        breaks = breaks,
        labels = labels,
        include.lowest = TRUE
      )
    )
  
  return(data)
}

interpro_plot <- function(ioi, interpro_file, uniprot_ids, peptide_seq_list_file = NULL, fasta_file = NULL, abundance_file = NULL, abundance_column = "quant",
                          correlation_palette = c("red", "gray", "blue"), 
                          interpro_colors = c("#FF5376", "#72AFD9", "#E3D26F", "#A288E3", "#1B5299", "#68D8D6", "#B78DA3")) {
  
  interpro_ioi <- readRDS(interpro_file)
  interpro_ioi <- interpro_ioi[interpro_ioi$hgnc_symbol == ioi & 
                                 !is.null(interpro_ioi$interpro_description) & 
                                 !is.na(interpro_ioi$interpro_description),]
  
  ioi_alphafold <- tryCatch({get_alphafold_file(ioi, uniprot_ids)},
                            error = function(e){0})
  
  # if(!is.null(fasta_file)){
  #   fasta_list <- readRDS(fasta_file)
  #   fasta_list <- fasta_list[[which(names(fasta_list) %in% c(ioi, unique(uniprot_ids[uniprot_ids$hgnc_symbol == ioi,]$uniprotswissprot)))]][[1]] 
  # } else {
  #   # Get FASTA from alphafold
  #   fasta_list <- get_FASTA_fromalphafold(ioi_alphafold)
  # }
  
  # Get fasta seq from plasma_data_FASTA_by_gene.RDS (pass as fasta_file)
  fasta_list <- unique(fasta_file$fasta[fasta_file$Gene.Name == ioi])
  
  peptide_seq_list <- readRDS(peptide_seq_list_file)
  peptide_seq_list <- peptide_seq_list[[ioi]]
  peptide_indices <- tryCatch({get_peptide_and_correlation_numeric(fasta_list, peptide_seq_list)},
                              error = function(e){0})
  
  interpro_options <- unique(interpro_ioi$interpro_description)
  
  # Initialize meta_hmap with two columns: Amino acid residue, Olink target correlation
  meta_hmap <- matrix(nrow = nchar(fasta_list[[1]]), ncol = 2)
  colnames(meta_hmap) <- c("Amino acid residue", "Olink target correlation")
  meta_hmap[,'Amino acid residue'] <- as.numeric(1:nchar(fasta_list[[1]]))
  
  
  # Function to process peptides and handle overlaps (as before)
  process_peptides <- function(col) {
    overlap_peptides <- list()
    for (k in seq_along(peptide_indices)) {
      start <- as.numeric(peptide_indices[[k]][["start_position"]])
      end <- start + nchar(peptide_indices[[k]][["peptide"]]) - 1
      added_to_overlap <- FALSE
      cat("k =", k, "peptide =", peptide_indices[[k]][["peptide"]], "\n")
      
      before_len <- nchar(fasta_list)
      after_sub  <- gsub(paste0(peptide_indices[[k]][["peptide"]], ".*"), "", fasta_list)
      cat("Before length =", before_len, 
          "After length =", nchar(after_sub), 
          "start_position =", nchar(after_sub) + 1, 
          "\n"
      )
      
      if (is.na(start) || is.na(end) || start < 1 || end > nrow(meta_hmap) || start > end) {
        # skip this peptide if it lacks valid mapping
        # some lack valid mapping!!!!
        next
      }
      if(all(is.na(meta_hmap[start:end, col]))){
        meta_hmap[start:end, col] <- rep(peptide_indices[[k]][["target_correlation"]], times = length(start:end))
      } else {
        added_to_overlap = TRUE
      }
      
      if (added_to_overlap) {
        overlap_peptides[[length(overlap_peptides) + 1]] <- peptide_indices[[k]]
      }
    }
    return(list(overlap_peptides=overlap_peptides, meta_hmap=meta_hmap))
  }
  
  # Manage overlaps, starting with the correlation column (2)
  col <- 2
  repeat {
    result <- process_peptides(col)
    current_overlaps <- result$overlap_peptides
    meta_hmap <- result$meta_hmap
    if (length(current_overlaps) == 0) {
      break
    }
    peptide_indices <- current_overlaps
    # If more overlaps remain, add another correlation column
    col <- ncol(meta_hmap) + 1
    meta_hmap <- cbind(meta_hmap, NA)
    colnames(meta_hmap)[col] <- paste("Correlation", col - 2)
  }
  
  meta_hmap <- as.data.frame(meta_hmap, stringsAsFactors = FALSE)
  
  # Add domain columns
  for(domain in interpro_options) {
    meta_hmap[, domain] <- NA
  }
  
  # Assign domain annotations
  for(i in meta_hmap[,'Amino acid residue']) {
    aa_idx <- which(meta_hmap[,'Amino acid residue'] == i)
    for(j in seq_len(nrow(interpro_ioi))) {
      if(i >= interpro_ioi$interpro_start[j] & i <= interpro_ioi$interpro_end[j]) {
        dname <- interpro_ioi$interpro_description[j]
        meta_hmap[aa_idx, dname] <- dname
      }
    }
  }
  
  # We now have: 
  # "Amino acid residue", "Olink target correlation", maybe more "Correlation X" columns, and domain columns.
  
  # Convert this to numeric
  meta_hmap[,"Olink target correlation"] <- as.numeric(meta_hmap[,"Olink target correlation"])
  
  cor_data <- meta_hmap[,grepl("correlation", colnames(meta_hmap))|grepl("Correlation", colnames(meta_hmap)), drop=FALSE]
  # Convert cor_data into a numeric matrix for heatmaply
  cor_mat <- as.matrix(cor_data)
  rownames(cor_mat) <- meta_hmap[,"Amino acid residue"]
  
  # Determine domain feature category per residue
  domain_matrix <- meta_hmap[,interpro_options, drop=FALSE]
  domains_per_res <- apply(domain_matrix, 1, function(x) sum(!is.na(x)))
  
  feature_category <- character(nrow(meta_hmap))
  for (r in seq_len(nrow(meta_hmap))) {
    n_doms <- domains_per_res[r]
    if(n_doms == 0) {
      feature_category[r] <- "No Feature"
    } else if(n_doms == 1) {
      # Extract the single non-NA domain
      dom_name <- domain_matrix[r,][!is.na(domain_matrix[r,])]
      feature_category[r] <- dom_name
    } else {
      feature_category[r] <- "Multiple Features"
    }
  }
  
  feature_category <- factor(feature_category, levels = c("No Feature", sort(setdiff(unique(feature_category), c("No Feature","Multiple Features"))), "Multiple Features"))
  
  # # Assign colors to domains consistently
  # # Ensure we have enough colors
  # all_categories <- levels(feature_category)
  # n_cats <- length(all_categories)
  # if(length(interpro_colors) < n_cats) {
  #   interpro_colors <- colorRampPalette(interpro_colors)(n_cats)
  # }
  # cat_color_mapping <- setNames(interpro_colors[1:n_cats], all_categories)
  # 
  # Prepare hover text
  # Hover: residue number, correlation, feature category
  hover_text <- matrix(nrow = dim(cor_mat)[1], ncol = dim(cor_mat)[2])
  
  colnames(hover_text) <- colnames(cor_mat)
  
  for(j in colnames(cor_mat)){
    
    hover_text[,j] <- paste0("Residue: ", meta_hmap[,"Amino acid residue"],
                             "<br>Correlation: ", round(cor_mat[,j], 3),
                             "<br>Feature: ", feature_category)
    
  }
  
  # row_side_colors: we can provide a data frame with one column: feature_category
  # heatmaply will display these as side colors.
  row_side_colors_df <- data.frame(Feature = feature_category, stringsAsFactors = FALSE)
  rownames(row_side_colors_df) <- rownames(cor_mat)
  
  # Create a continuous color scale for correlation from -1 to 1
  # Use viridis or the publication palette
  # Already set scale -1 to 1:
  color_breaks <- seq(-1, 1, length.out = 100)
  correlation_palette <- colorRampPalette(correlation_palette)(100)
  
  #Substitute the domains for a "no feature" annotation instead of NA
  domain_matrix<- domain_matrix %>%
    #mutate(across(everything(), ~ ifelse(!is.na(.), "Feature", .))) %>%
    mutate(across(everything(), ~ ifelse(is.na(.), "", .))) 
  
  # Define color palette for domains
  domain_matrix_vec <- as.vector(as.matrix(domain_matrix))
  all_categories <- factor(unique(as.vector(domain_matrix_vec)))
  #cat_color_mapping <- setNames(interpro_colors[seq_along(all_categories)], all_categories)
  domain_colors <- colorRampPalette(colors = brewer.pal(9, 'PuBuGn')[2:8])(length(all_categories))
  
  cat_color_mapping <- vector("list", length(all_categories))
  names(cat_color_mapping) <- all_categories
  for (i in seq_along(all_categories)) {
    cat_color_mapping[[i]] <- domain_colors[i]
    names(cat_color_mapping[[i]]) <- all_categories[i]
  }
  
  
  # Total plot height
  domain_height <- (ncol(domain_matrix) + 1) * 30
  peptide_height <- (ncol(cor_mat) + 1) * 40
  total_height <- domain_height + peptide_height
  # Proportional heights for subplots
  domain_height <- domain_height / total_height
  peptide_height <- peptide_height / total_height
  
  # Cap height at 600 px
  total_height <- min(total_height, 600)
  
  # Min height of 300
  total_height <- max(total_height, 400)
  
  heatmap <- heatmaply(
    t(data.matrix(cor_mat)),
    #main = paste("InterPro Domain & Correlation -", ioi),
    scale = "none",
    limits = c(-1,1),
    showticklabels = c(TRUE, FALSE),
    colors = correlation_palette,
    ColSideColors = domain_matrix, 
    colv = F,
    Rowv = F,
    plot_method ='plotly',
    col_side_palette = colorRampPalette(brewer.pal(9, 'PuBuGn')[2:8]), # not working at the moment, works if you don't use "plotly" but then rest is broken
    hoverinfo = "text",
    subplot_heights = c(domain_height, peptide_height),
    height = total_height,
    custom_hovertext = t(hover_text), # for hover
    dendrogram = "none"
  ) %>%
    colorbar(
      title = "MS-Olink correlation", 
      titlefont = list(size = 10), 
      tickfont = list(size = 8),
      which = 2,
      x = 1,
      y = 0,
      len = 100,
      lenmode = 'pixels',
      yanchor = 'bottom',
      thickness = 15,
      yref = 'paper'
    )
  
  # Remove colorbar and legend for domain plot
  heatmap$x$data[[1]]$showscale <- FALSE
  
  # Help if the isoform is too short
  if (nrow(cor_mat) < 100) {
    tickvals <- integer(0)
  } else if (nrow(cor_mat) > 2000) {
    tickvals <- seq(1000, nrow(cor_mat), by = 1000)
  }else {
    tickvals <- seq(100, nrow(cor_mat), by = 100)
  }
  
  heatmap <- heatmap %>% 
    layout(
      xaxis = list(
        title = 'Amino Acid Position',
        titlefont = list(size = 10),
        tickfont = list(size = 8),
        tickvals = tickvals, #seq(100, nrow(cor_mat), by = 100),  # Set tick positions every 100
        ticktext = tickvals,  #seq(100, nrow(cor_mat), by = 100),   # Set corresponding labels
        tickangle = 0
      ),
      yaxis = list(
        titlefont = list(size = 10),
        tickfont = list(size = 8),
        side = 'right'
      ),
      yaxis2 = list(
        title = 'MS Peptides',
        titlefont = list(size = 10),
        tickfont = list(size = 8)
      )
    )
  
  return(heatmap)
  
}
# 
# interpro_plot <- function(ioi, interpro_file, uniprot_ids, peptide_seq_list_file = NULL, fasta_file = NULL, abundance_file = NULL, abundance_column = "quant",
#                           correlation_palette = c("red", "gray", "blue"), interpro_colors = c("#FF5376", "#72AFD9", "#E3D26F", "#A288E3", "#1B5299", "#68D8D6", "#B78DA3")){
# 
#   interpro_ioi <- readRDS(interpro_file)
#   interpro_ioi <- interpro_ioi[interpro_ioi$hgnc_symbol == ioi & !is.null(interpro_ioi$interpro_description) & !is.na(interpro_ioi$interpro_description),]
# 
#   ioi_alphafold <- tryCatch({get_alphafold_file(ioi, uniprot_ids)},
#                             error = function(e){0})
# 
#   if(!is.null(fasta_file)){fasta_list <- readRDS(fasta_file)
# 
#   fasta_list <- fasta_list[[which(names(fasta_list) %in% c(ioi, unique(uniprot_ids[uniprot_ids$hgnc_symbol == ioi,]$uniprotswissprot)))]][[1]]}else{
#     #read a FASTA from uniprot
#     #refer to alpha fold so FASTA is consistent
#     fasta_list <- get_FASTA_fromalphafold(ioi_alphafold)
# 
#   }
# 
#   peptide_seq_list <- readRDS(peptide_seq_list_file)
#   peptide_seq_list <- peptide_seq_list[[ioi]]
#   peptide_indices <- tryCatch({get_peptide_and_correlation_numeric(fasta_list, peptide_seq_list)},
#                               error = function(e){0})
#   #retrieve from alphafold so it's consistent with the other visualizations in 3D
# 
#   interpro_options <- unique(interpro_ioi$interpro_description)
# 
# 
#   meta_hmap <- matrix(nrow = nchar(fasta_list[[1]]), ncol = length(c("Amino acid residue", "Olink target correlation"#, interpro_options
#   )))
#   colnames(meta_hmap) <- c("Amino acid residue", "Olink target correlation"#, interpro_options
#   )
#   meta_hmap[,'Amino acid residue'] <- as.numeric(1:nchar(fasta_list[[1]]))
# 
# 
#   # Function to process peptides and handle overlaps
#   process_peptides <- function(col) {
#     overlap_peptides <- list()
# 
#     for (k in seq_along(peptide_indices)) {
#       start <- as.numeric(peptide_indices[[k]][["start_position"]])
#       end <- start + nchar(peptide_indices[[k]][["peptide"]]) - 1
#       added_to_overlap <- FALSE
# 
#       if(sum(!is.na(meta_hmap[start:end, col]))==0){
#         meta_hmap[start:end, col] <- rep(peptide_indices[[k]][["target_correlation"]],times = length(start:end))
#       } else {
#         added_to_overlap = TRUE
#       }
# 
#       if (added_to_overlap) {
#         overlap_peptides[[length(overlap_peptides) + 1]] <- peptide_indices[[k]]
#       }
#     }
# 
#     return(list(overlap_peptides=overlap_peptides, meta_hmap=meta_hmap))
# 
#   }
# 
# 
#   # Main loop to manage and resolve overlaps
# 
#   col = 2
#   repeat {
#     result <- process_peptides(col)
#     current_overlaps <- result$overlap_peptides
#     meta_hmap <- result$meta_hmap
#     if (length(current_overlaps) == 0) {
#       break  # Exit if no overlaps remain
#     }
#     # Prepare for the next round of overlaps
#     peptide_indices <- current_overlaps
#     col <- ncol(meta_hmap) + 1
#     meta_hmap <- cbind(meta_hmap, NA)
#     colnames(meta_hmap)[col] <- paste("Correlation", col)
#   }
# 
#   meta_hmap <- as.data.frame(meta_hmap)
#   colnames_corrrelation <- colnames(meta_hmap)[colnames(meta_hmap) != "Amino acid residue"]
# 
# 
#   for(domain in interpro_ioi$interpro_description){
#     meta_hmap[,domain] <- NA
#   }
#   for(i in as.numeric(meta_hmap[,'Amino acid residue'])){
# 
#     for(j in 1:nrow(interpro_ioi)){
# 
#       if(i >= interpro_ioi$interpro_start[j] & i <= interpro_ioi$interpro_end[j]){
# 
#         meta_hmap[meta_hmap[,'Amino acid residue'] == i,interpro_ioi$interpro_description[j]] <- interpro_ioi$interpro_description[j]
# 
#       }
# 
#     }
#   }
# 
# 
#   #make color for multi mapping
# 
#   # min_correlation <- min(meta_hmap[,'Olink target correlation'], na.rm = TRUE)
#   # max_correlation <- max(meta_hmap[,'Olink target correlation'], na.rm = TRUE)
#   # max_abs_correlation <- max(abs(min_correlation), abs(max_correlation))
# 
#   color_breaks <- seq(-1, 1, length.out = 100)
#   correlation_palette <- colorRampPalette(correlation_palette)(length(color_breaks))
#   color_mapping_function <- colorRamp2(color_breaks, correlation_palette)
# 
# 
#   named_palette_list = setNames(rep(list(color_mapping_function), length(colnames_corrrelation)), colnames_corrrelation)
# 
# 
#   #colors for each interpro option
# 
#   if (length(interpro_colors) < length(interpro_options)) {
#     interpro_colors <- rep(interpro_colors, length.out = length(interpro_options))
#   }
# 
#   domain_color_mapping <- setNames(interpro_colors, interpro_options)
# 
#   for (domain in interpro_options) {
#     named_palette_list[[domain]] <- domain_color_mapping[domain]
#   }
# 
# 
#   #make it horizontal
#   #meta_hmap <- as.data.frame(t(meta_hmap))
# 
#   #set legend things to hide
#   legend_logical <- logical(length = length(colnames(meta_hmap)))
#   legend_logical[1] <- TRUE
#   legend_logical[2] <- TRUE
# 
#   if(length(colnames(meta_hmap)) >= 3){for (i in 3:length(colnames(meta_hmap))){
# 
#     legend_logical[i] <- FALSE
# 
#   }
# 
#     metadata_annotation_obj <- HeatmapAnnotation(df = meta_hmap, na_col = "white", col =   named_palette_list,
#                                                  show_legend = legend_logical,
#                                                  which = "col")
# 
#     plot(metadata_annotation_obj)
#     return(metadata_annotation_obj)
#   }
# 
# }
# 


jenks_density_plot <- function(ioi, peptide_seq_list_file, uniprot_ids, fasta_file = NULL) {
  library(ggplot2)
  library(dplyr)
  
  # Load peptide sequences for the gene of interest
  peptide_seq_list <- readRDS(peptide_seq_list_file)
  if (!ioi %in% names(peptide_seq_list)) {
    stop("IOI not found in peptide_seq_list.")
  }
  
  gene_peptides <- peptide_seq_list[[ioi]]
  
  # Get FASTA sequence
  # if(!is.null(fasta_file)){
  #   fasta_all <- readRDS(fasta_file)
  #   # Extract relevant fasta; uniprot_ids should link gene_symbol to uniprot
  #   relevant_uniprot <- unique(uniprot_ids[uniprot_ids$hgnc_symbol == ioi,]$uniprotswissprot)
  #   relevant_uniprot <- relevant_uniprot[relevant_uniprot %in% names(fasta_all)]
  #   if(length(relevant_uniprot) == 0) stop("No matching uniprot found for ioi in fasta_file.")
  #   
  #   # Use the first matching fasta sequence
  #   fasta_seq <- fasta_all[[relevant_uniprot[1]]]
  #   # Ensure fasta_seq is a character
  #   fasta_seq <- as.character(fasta_seq)
  # } else {
  #   # If no fasta_file is given, try alphafold route (assuming get_alphafold_file and get_FASTA_fromalphafold defined)
  #   ioi_alphafold <- get_alphafold_file(ioi, uniprot_ids)
  #   fasta_seq <- get_FASTA_fromalphafold(ioi_alphafold)
  #   fasta_seq <- as.character(fasta_seq[[1]])
  # }
  # 
  
  # Get FASTA sequence from plasma_data_FASTA_by_gene.RDS (pass in fasta_file)
  fasta_seq <- unique(fasta_file$fasta[fasta_file$Gene.Name == ioi])
  
  jenks_df <- gene_peptides %>% 
    mutate(peptide_id = row_number()) %>% 
    mutate(jenks_class = cor_intervals(as.numeric(target_correlation)))
  
  # jenks_df now has a jenks_class column
  # Now, we need to map each peptide to its amino acid positions
  
  # For each peptide, we find its start position in the sequence
  # The start position is computed by removing the peptide and counting length
  # The first occurrence of the peptide in fasta_seq gives start position
  # But be aware of multiple occurrences. If multiple matches occur, handle accordingly!!!!! RECHECK THIS LOGIC
  
  # Assume the first occurrence is the correct one:
  create_position_data <- function(peptide, jenks_class, fasta_seq, used_positions) {
    # Find the start position of the peptide in fasta_seq
    # We can use gregexpr for this:
    match_pos <- gregexpr(peptide, fasta_seq, fixed = TRUE)[[1]]
    if (match_pos[1] == -1) {
      # Peptide not found, skip
      return(NULL)
    }
    
    # Find the first match that does not overlap with used_positions
    for (pos in match_pos) {
      start_pos <- pos
      end_pos <- start_pos + nchar(peptide) - 1
      if (!any(used_positions >= start_pos & used_positions <= end_pos)) {
        # Mark these positions as used
        used_positions[start_pos:end_pos] <- TRUE
        # Create a data frame of all positions covered by this peptide
        return(data.frame(position = start_pos:end_pos,
                          jenks_class = jenks_class,
                          stringsAsFactors = FALSE))
      }
    }
    # If all matches overlap, skip
    return(NULL)
  }
  
  # Initialize a vector to keep track of used positions
  used_positions <- rep(FALSE, nchar(fasta_seq))
  
  # Apply this to all peptides
  position_data_list <- lapply(seq_len(nrow(jenks_df)), function(i) {
    current_peptide <- jenks_df$peptides[i]
    current_class <- jenks_df$jenks_class[i]
    
    # Check if peptide is non-empty
    if (is.na(current_peptide) || current_peptide == "") return(NULL)
    create_position_data(current_peptide, current_class, fasta_seq, used_positions)
  })
  
  # Combine all
  position_data <- do.call(rbind, position_data_list)
  
  # If no data, return NULL or a simple message
  if (is.null(position_data) || nrow(position_data) == 0) {
    warning("No peptide position data available for plotting.")
    return(NULL)
  }
  
  # Now we have a data frame with columns: position, jenks_class
  # We can make a density plot by jenks_class
  p <- ggplot(position_data, aes(x = position, fill = jenks_class,
                                 text = paste("Class:", jenks_class))) +
    geom_density(alpha = 0.5, lwd = 0.3) +
    # scale_fill_manual(values = c("Negative correlation" = "#4E79A7", 
    #                              "No correlation" = "gray",
    #                              "Weak correlation" = "#F28E4C",  
    #                              "Moderate correlation" = "#E15759", 
    #                              "Strong correlation" = "#F28BBB")) +
    # scale_fill_manual(values = c("No correlation" = "grey",
    #                              "Weak correlation" = "#9ccfe7",  
    #                              "Moderate correlation" = "#629db8", 
    #                              "Strong correlation" = "#206e8c")) +
    scale_fill_manual(values = categorical_colors) +
    theme_minimal() +
    labs(
      title = paste("Peptide Density by Jenks Class -", ioi),
      x = "Amino Acid Position",
      y = "Density",
      fill = "Correlation Category"
    ) +
    theme(
      legend.title = element_blank()
    )
  
  return(p)
}



# Combined Plot Function
combined_interpro_density_plot <- function(ioi, interpro_file, uniprot_ids, 
                                           peptide_seq_list_file = NULL, fasta_file = NULL, 
                                           abundance_file = NULL, abundance_column = "quant",
                                           correlation_palette = c("red", "gray", "blue"), 
                                           interpro_colors = c("#FF5376", "#72AFD9", "#E3D26F", 
                                                               "#A288E3", "#1B5299", "#68D8D6", "#B78DA3")) {
  
  # Generate the Heatmaply Plot
  heatmap_plot <- interpro_plot(
    ioi = ioi,
    interpro_file = interpro_file,
    uniprot_ids = uniprot_ids,
    peptide_seq_list_file = peptide_seq_list_file,
    fasta_file = fasta_file,
    abundance_file = abundance_file,
    abundance_column = abundance_column,
    correlation_palette = correlation_palette,
    interpro_colors = interpro_colors
  )
  
  # Generate the Density Plot using ggplot2
  density_plot_gg <- jenks_density_plot(
    ioi = ioi,
    peptide_seq_list_file = peptide_seq_list_file,
    uniprot_ids = uniprot_ids,
    fasta_file = fasta_file
  )
  
  # Convert ggplot2 Density Plot to Plotly Object
  density_height <- 150
  if(!is.null(density_plot_gg)){
    density_plotly <- ggplotly(density_plot_gg, tooltip = "text", height = density_height) %>%
      layout(showlegend = TRUE,
             legend = list(
               itemwidth = 15,
               tracegroupgap = 1
             )
      )
    
  } else {
    # If density plot is NULL, proceed with heatmap only
    return(heatmap_plot)
  }
  
  # Remove x-axis title and tick labels from density plot to align with heatmap
  density_plotly <- layout(density_plotly, xaxis = list(title = "", showticklabels = FALSE))
  
  # Ensure the heatmap and density plot share the same x-axis range
  heatmap_xrange <- heatmap_plot$x$layout$xaxis$range
  density_plotly <- layout(density_plotly, xaxis = list(range = heatmap_xrange))
  
  # Total height of the heatmap based on the number of rows
  heatmap_height <- heatmap_plot$height
  
  # Total plot height should be 800 px
  total_height <- 800
  
  # Combine Heatmap and Density Plot using subplot
  combined_plot <- subplot(
    density_plotly,
    heatmap_plot,
    nrows = 2, 
    shareX = TRUE, 
    margin = 0,
    heights = c(density_height/total_height, heatmap_height/total_height)
  ) %>%
    layout(
      title = list(
        text = paste("MS Peptide Mapping & Correlation to Olink -", ioi),
        font = list(size = 12)),
      xaxis = list(title = "Amino Acid Position"),
      yaxis = list(showticklabels = FALSE),
      margin = list(t = 30, b = 25, l = 50, r = 0),
      legend = list(font = list(size = 8), title = list(font = list(size = 10))),
      annotations = list(
        list(
          x = -0.01, y = 0.5,
          showarrow = FALSE,
          text = 'Density',
          xref = 'paper', yref = 'y1 domain',
          xanchor = 'right', yanchor = 'middle',
          font = list(size = 10),
          textangle = -90
        ),
        list(
          x = -0.01, y = 0.5,
          showarrow = FALSE,
          text = 'Protein Domains',
          xref = 'paper', yref = 'y3 domain',
          xanchor = 'right', yanchor = 'middle',
          font = list(size = 10),
          textangle = -90
        ),
        list(
          x = -0.01, y = 0.5,
          showarrow = FALSE,
          text = 'MS Peptides',
          xref = 'paper', yref = 'y2 domain',
          xanchor = 'right', yanchor = 'middle',
          font = list(size = 10),
          textangle = -90
        )
      )
    )
  
  return(combined_plot)
}

expand_semicolon_rows <- function(dt, columns_to_expand, sep = ";") {
  # dt:                data.table to expand
  # columns_to_expand: character vector of column names to split by `sep`
  # sep:               the delimiter to split on, default ";"
  
  # Make a list to hold expanded rows
  result_list <- vector("list", nrow(dt))
  
  # Iterate over each row
  for (i in seq_len(nrow(dt))) {
    # Extract the single row as a data.table
    row_data <- dt[i]
    
    # Split each specified column into a character vector
    splitted_cols <- lapply(columns_to_expand, function(col) {
      # Convert to character to avoid factor issues
      strsplit(as.character(row_data[[col]]), sep, fixed = TRUE)[[1]]
    })
    
    # How many times do we need to replicate this row?
    lengths_vec <- sapply(splitted_cols, length)
    n_dups <- max(lengths_vec)
    
    # Create a copy of the row repeated n_dups times
    expanded_dt <- row_data[rep(1, n_dups)]
    
    # Fill in each column with the single split values
    for (j in seq_along(columns_to_expand)) {
      colname <- columns_to_expand[j]
      split_vals <- splitted_cols[[j]]
      
      # If this particular column had fewer splits than n_dups, 
      # fill-forward the last entry
      if (length(split_vals) < n_dups) {
        split_vals <- c(split_vals, rep(tail(split_vals, 1), n_dups - length(split_vals)))
      }
      
      expanded_dt[[colname]] <- split_vals
    }
    
    # Store the expanded data.table for this row
    result_list[[i]] <- expanded_dt
  }
  
  # Combine everything into one expanded data.table
  return(data.table::rbindlist(result_list, use.names = TRUE, fill = TRUE))
}



cor_intervals <- function(x) {
  cut(x, breaks = c(-1, 0.3, 0.5, 0.7, 1),
      labels = c('No correlation','Weak correlation', 
                 'Moderate correlation', 'Strong correlation'),
      include.lowest = TRUE, right = FALSE)
}

# Helper: Compute median correlation per residue based on peptide indices
get_residue_correlation_median <- function(fasta_seq, peptide_indices) {
  n <- nchar(fasta_seq)
  # Create a list to hold correlations for each residue
  correlation_list <- vector("list", n)
  
  # For each peptide, add its target correlation to every residue it covers
  for (pep in peptide_indices) {
    start <- as.numeric(pep$start_position)
    pep_len <- nchar(pep$peptide)
    end <- start + pep_len - 1
    if (!is.na(start) && start > 0 && end <= n) {
      for (pos in start:end) {
        correlation_list[[pos]] <- c(correlation_list[[pos]], pep$target_correlation)
      }
    }
  }
  
  # Compute the median for each residue (if no values, return NA)
  median_corr <- sapply(correlation_list, function(x) {
    if (length(x) > 0) median(x, na.rm = TRUE) else NA
  })
  
  return(median_corr)
}

# Helper: Map numeric correlation values to colors using a continuous palette
map_correlation_to_color <- function(corr_values, palette = c("red", "gray", "blue"), n_bins = 100) {
  # Create a vector of colors based on the provided palette
  palette_colors <- colorRampPalette(palette)(n_bins)
  # Define breaks assuming correlation values range from -1 to 1
  breaks <- seq(-1, 1, length.out = n_bins + 1)
  # Bin each correlation value; if NA, it remains NA
  bins <- cut(corr_values, breaks = breaks, include.lowest = TRUE)
  colors <- ifelse(is.na(bins), NA, palette_colors[as.numeric(bins)])
  return(colors)
}

alphafold_plot <- function(ioi, genesymb, uniprot_ids, 
                           peptide_seq_list_file = NULL, 
                           fasta_file = NULL, 
                           correlation_palette = c("red", "gray", "blue")) {
  # Load peptide sequence list for the isoform
  if (!is.null(peptide_seq_list_file)) {
    peptide_seq_list <- readRDS(peptide_seq_list_file)
    peptide_seq_list <- peptide_seq_list[[genesymb]]
  } else {
    stop("No peptide membership data provided.")
  }
  
  # Get the AlphaFold structure file using the isoform directly
  ioi_alphafold <- tryCatch({
    get_alphafold_file(ioi)
  }, error = function(e) { 0 })
  
  if (!is.list(ioi_alphafold) || is.null(ioi_alphafold[["status_code"]])) {
    stop("No valid AlphaFold structure file returned for this ID of interest.")
  }
  
  if (ioi_alphafold[["status_code"]] == 404) {
    stop("No AlphaFold structure was found for this ID of interest.")
  }
  
  # Obtain FASTA sequence: either from file or via AlphaFold retrieval
  if (!is.null(fasta_file)) {
    fasta_list <- readRDS(fasta_file)
    # Assume fasta_list is a named list with ioi as key
    fasta_seq <- fasta_list[[ioi]]
  } else {
    fasta_seq <- get_FASTA_fromalphafold(ioi_alphafold)
  }
  
  # Process peptides to get indices and correlation values
  peptide_indices <- tryCatch({
    get_peptide_and_correlation_numeric(fasta_seq, peptide_seq_list)
  }, error = function(e) { 0 })
  
  # If no valid peptide mapping, return a default cartoon representation
  if (!is.list(peptide_indices)) {
    p <- NGLVieweR(ioi_alphafold[["content"]]) %>%
      addRepresentation("cartoon", 
                        param = list(name = "cartoon", backgroundColor = "white", 
                                     colorScheme = "uniform", colorValue = "#808080")) %>%
      stageParameters(backgroundColor = "white", colorValue = "#808080") %>%
      setQuality("high") %>%
      setFocus(0) %>%
      setSpin(TRUE)
    return(p)
  }
  
  # Compute median correlation per residue and map to colors
  residue_corr <- get_residue_correlation_median(fasta_seq, peptide_indices)
  residue_colors <- map_correlation_to_color(residue_corr, palette = correlation_palette)
  
  # Start with a base representation of the structure using "tube" style
  p <- NGLVieweR(ioi_alphafold[["content"]]) %>%
    addRepresentation("tube", 
                      param = list(name = "tube", backgroundColor = "white", 
                                   colorScheme = "uniform", colorValue = "#808080")) %>%
    stageParameters(backgroundColor = "white", colorValue = "#808080") %>%
    setQuality("low") %>%
    setFocus(0) %>%
    setSpin(TRUE)
  
  # Instead of grouping residues, add a surface representation for each residue (disjoint)
  seq_length <- nchar(fasta_seq)
  for (pos in 1:seq_length) {
    col_val <- residue_colors[pos]
    if (!is.na(col_val)) {
      # Selection string for a single residue position
      sele_str <- as.character(pos)
      p <- tryCatch({
        p %>% addRepresentation("surface", 
                                param = list(colorScheme = "uniform",
                                             colorValue = col_val,
                                             sele = sele_str,
                                             opacity = 0.25))
      }, error = function(e) { p })
    }
  }
  
  # Clean up temporary AlphaFold file and return the NGL view
  unlink(ioi_alphafold[["content"]])
  return(p)
}

color_scale_plot <- function(correlation_palette, n_bins = 100) {
  # Create a data frame with evenly spaced values from -1 to 1
  df <- data.frame(x = seq(-1, 1, length.out = n_bins), y = 1)
  
  p <- ggplot(df, aes(x = x, y = y, fill = x)) +
    geom_tile() +
    scale_fill_gradientn(colors = colorRampPalette(correlation_palette)(n_bins),
                         limits = c(-1, 1)) +
    scale_x_continuous(breaks = seq(-1, 1, by = 0.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    theme_minimal(base_family = "Open Sans") +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          panel.grid   = element_blank(),
          legend.position = "none",
          plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
    labs(x = "Correlation", fill = NULL)
  
  return(p)
}



#  _______ _    _ _____  _____ 
# |__   __| |  | |_   _|/ ____|
#    | |  | |__| | | | | (___  
#    | |  |  __  | | |  \___ \ 
#    | |  | |  | |_| |_ ____) |
#  __|_|  |_|__|_|_____|_____/ 
# |_   _|/ ____|               
#   | | | (___                 
#   | |  \___ \                
#  _| |_ ____) |               
# |_____|_____/ _ ______       
# |__   __| |  | |  ____|      
#    | |  | |__| | |__         
#    | |  |  __  |  __|        
#    | |  | |  | | |____       
#    |_|  |_|__|_|______|      
#     /\   |  __ \|  __ \      
#    /  \  | |__) | |__) |     
#   / /\ \ |  ___/|  ___/      
#  / ____ \| |    | |          
# /_/    \_\_|    |_|          


ui <- fluidPage(
  theme = my_theme,
  tags$head(
    tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Fascinate+Inline&display=swap"),
    tags$link(rel = "icon", type = "image/png", href = "icon.png")
    ),
  tags$style(HTML("
h1, h2, h3, h4 {
      font-family: 'Fascinate Inline', cursive;
      color: #ff5c8d; /* Vibrant pink-orange color */
}
    .shiny-plot-output, .plotly {
        border-radius: 15px;
        overflow: hidden;
    }
      .nav-tabs + .tab-content {
      background-color: #fff;
      padding: 15px; /* optional, for spacing */
      }
    .nav-tabs .nav-item .nav-link.active {
      background-color: #FEEEEB;
      border-color: #fff;
    }
    .info-icon {
      font-size: 0.9em; 
      color: #666; 
      margin-left: 5px; 
      cursor: help;
    }
    .plot-title {
      font-size: 1rem; 
      font-weight: 600;
      margin-bottom: 0.5em;
    }
    .sidebar-section {
      margin-bottom: 15px;
    }
  ")),
  titlePanel(
    "PeptOlink",
    windowTitle = "PeptOlink"
  ),
  fluidRow(
    column(3,
           div(class = "sidebar-section",
               # Filters
               dropdownButton(
                 inputId = "filters_dropdown_btn",
                 label = "Filters",
                 icon = icon("sliders-h"),
                 status = "primary",
                 circle = FALSE,
                 inline = TRUE,
                 
                 numericInput("mean_corr", "Mean correlation ≤", 
                              value = 1, max = 1, min = -1, step = .01),
                 numericInput("sd_corr", "Correlation SD ≥", 
                              value = 0, max = 1, min = 0, step = .01),
                 numericInput("median_corr", "Median correlation ≤",
                              value = 1, max = 1, min = -1, step = .01),
                 numericInput("iqr_corr", "Correlation IQR ≥",
                              value = 0, max = 1, min = 0, step = .01),
                 numericInput("n_peptides", "Number of peptides ≥", 
                              value = 5, max = 500, min = 1, step = 5),
                 numericInput("n_isoforms", "Number of isoforms ≥", 
                              value = 1, max = 4, min = 1, step = 1),
                 actionBttn("clear_filters", "Clear", icon = icon("times"), style = "unite", size = "xs")
               )
           ),
           div(class = "sidebar-section",
               uiOutput("select_result_ui"),
               uiOutput("isoform_count_ui"),
               uiOutput("isoform_select_ui")
           )
    ),
    column(9,
           tabsetPanel(
             id = "main_tabs",
             
             # Detailed first
             # Replace plotOutput with plotlyOutput for interactive plots
             tabPanel("Detailed",
                      fluidRow(
                        column(12,
                               span("Isoform-specific Data", class = "plot-title"),
                               icon("info-circle", class = "info-icon", id = "detailed_info")
                        )
                      ),
                      br(),
                      fluidRow(
                        column(12,
                               withSpinner(plotlyOutput("detailed_plot", height = "800px"), type = 4) 
                        )
                      )
             ),
             
             tabPanel("Structural",
                      fluidRow(
                        column(12,
                               
                               span("Alphafold Structural Data", class = "plot-title"),
                               icon("info-circle", class = "info-icon", id = "alphafold_info")
                        )
                      ),
                      fluidRow(
                        column(12,
                               div(style = "display: flex; justify-content: center;",
                                   uiOutput("alphafold_warn")  
                               )
                        )
                      ),
                      fluidRow(
                        column(12,
                               withSpinner(NGLVieweROutput("NGL_plot"), type = 4) 
                        )
                      ),

                      br(),
                      fluidRow(
                        column(12,
                               div(style = "display: flex; justify-content: center;",
                                   plotOutput("color_scale", width = "60%", height = "80px")
                               )
                        )
                      )
             ),
             
             
             # Summary second
             tabPanel("Summary",
                      fluidRow(
                        column(12,
                               span("Summary of Filtered Results", class = "plot-title"),
                               icon("info-circle", class = "info-icon", id = "summary_info")
                        )
                      ),
                      br(),
                      fluidRow(
                        column(12, withSpinner(plotOutput("summary_plot", height = "300px"), type = 4))
                      ),
                      
                      br(),
                      
                      fluidRow(
                        column(6, withSpinner(plotlyOutput("volcano_plot", height = "300px"), type = 4)),
                        column(6, withSpinner(plotlyOutput("volcano_plot2", height = "300px"), type = 4))
                      )
             )
           )
    )
  )
)

# Define server logic
server <- function(input, output, session) {

  #anova_results <- readRDS("data/ANOVA_CLUSTER_leiden_1_partition_peptides_bypeptide.RDS")
  #kruskal_results <- readRDS("data/KRUSKAL_CLUSTER_leiden_1_partition_peptides_bypeptide.RDS")
  
  gene_stats <- correlation_long_filt %>%
    group_by(gene_symbol) %>%
    summarise(
      mean_corr = mean(correlation, na.rm = TRUE),
      sd_corr = sd(correlation, na.rm = TRUE),
      median_corr = median(correlation, na.rm = TRUE),
      iqr_corr = IQR(correlation, na.rm = TRUE),
      n_peptides = n(),
      n_isoforms = length(unique(unlist(strsplit(paste(UniProt.MS, collapse = ";"), ";"))))
    )
  
  
  filtered_data <- reactive({
    
    # Ensure the correct inputs are available
    req(
      input$sd_corr,
      input$mean_corr,
      input$median_corr,
      input$iqr_corr,
      input$n_peptides,
      input$n_isoforms
    ) 
    
    # Filter genes based on user inputs
    filtered_genes <- gene_stats %>%
      filter(
        sd_corr >= input$sd_corr,
        mean_corr <= input$mean_corr,
        median_corr <= input$median_corr,
        iqr_corr >= input$iqr_corr,
        n_peptides >= input$n_peptides,
        n_isoforms >= input$n_isoforms
      ) %>%
      pull(gene_symbol)
    
    filtered <- correlation_long_filt %>%
      filter(gene_symbol %in% filtered_genes)
    
    # If no data after filtering, return NULL
    if (nrow(filtered) == 0) {
      return(NULL)
    }
    
    # Apply fixed correlation intervals instead of Jenks
    filtered <- filtered %>%
      mutate(jenks_class = cor_intervals(correlation))
    
    return(filtered)
  })
  
  
  # #Update filtered data when thresholds change
  #Not needed -- REACTIVE
  # observeEvent(input$update, {
  #   anova_filtered <- anova_results %>%
  #     filter(p.value <= input$anova_threshold)
  #   
  #   kruskal_filtered <- kruskal_results %>%
  #     filter(p.value <= input$kruskal_threshold)
  #   
  #   new_data <- correlation_long_filt %>%
  #     filter(gene_symbol %in% anova_filtered$gene_symbol,
  #            gene_symbol %in% kruskal_filtered$gene_symbol)
  #   
  #   filtered_data(new_data)  # Update the reactive value
  #   
  #   # Update UI and indicators
  #   session$sendCustomMessage("toggleFilterIndicator", TRUE)
  #   updateTabsetPanel(session, "main_tabs", selected = "Summary")
  #   session$sendCustomMessage("toggleLoading", FALSE)
  # })
  
  # Clear filters and reset data
  # observeEvent(input$clear_filters, {
  #   filtered_data(correlation_long_filt)  # Reset to the full dataset
  #   
  #   # Reset thresholds
  #   updateNumericInput(session, "anova_threshold", value = 1)
  #   updateNumericInput(session, "kruskal_threshold", value = 1)
  #   
  #   # Update UI and indicators
  #   session$sendCustomMessage("toggleFilterIndicator", FALSE)
  #   updateTabsetPanel(session, "main_tabs", selected = "Summary")
  #   session$sendCustomMessage("toggleLoading", FALSE)
  # })
  observeEvent(input$clear_filters, {
    
    # Reset thresholds
    # Goal: Show everything
    updateNumericInput(session, "sd_corr", value = 0)
    updateNumericInput(session, "mean_corr", value = 1)
    updateNumericInput(session, "median_corr", value = 1)
    updateNumericInput(session, "iqr_corr", value = 0)
    updateNumericInput(session, "n_peptides", value = 1)
    updateNumericInput(session, "n_isoforms", value = 1)
   
  })
  
  # UI for selecting results
  output$select_result_ui <- renderUI({
    fd <- filtered_data()
    if (is.null(fd) || nrow(fd) == 0) {
      return(tags$div("No results match the current filters.", class="text-danger"))
    }
    choices <- unique(fd$gene_symbol)
    selectInput("selected_result", "Select Gene:", choices = choices, selected = choices[1])
  })
  
  observeEvent(input$selected_result, {
    plasma_data <- plasma_dt[.(input$selected_result)] # Fast subset
  })
  
  #Isoform -------
  output$isoform_count_ui <- renderUI({
    req(input$selected_result)
    #Subset to the selected gene
    sub_dt <- plasma_dt[Gene.Name == input$selected_result]
    if (nrow(sub_dt) == 0) return(NULL)
    
    #Collect unique isoforms (can split on ";" for multi matching)
    isoforms <- unique(unlist(strsplit(sub_dt$UniProt.MS, ";")))
    num_iso <- length(isoforms)
    
    # Show the info in a gray box (inline style example)
    div(
      style = "background-color: #f9f9f9; border: 1px solid #ccc; padding: 5px; margin-bottom: 10px;",
      paste0("For ", input$selected_result, ", there are ", num_iso, " isoforms detected.")
    )
  })
  
  output$isoform_select_ui <- renderUI({
    req(input$selected_result)
    sub_dt <- plasma_dt[Gene.Name == input$selected_result]
    if (nrow(sub_dt) == 0) return(NULL)
    
    isoforms <- unique(unlist(strsplit(sub_dt$UniProt.MS, ";")))
    # Safely select first isoform if available
    default_selected <- if (length(isoforms) > 0) isoforms[1] else NULL
    
    selectInput("selected_isoforms", 
                "Select Isoform(s):", 
                choices = isoforms, 
                selected = default_selected, 
                multiple = FALSE)
  })
  
  working_plasma_dt <- reactive({
    req(input$selected_result)  # Must have a gene
    dt <- plasma_dt[Gene.Name == input$selected_result]
    
    # If isoforms are chosen, subset further
    if (!is.null(input$selected_isoforms) && length(input$selected_isoforms) > 0) {
      
      cols_to_expand <- c("UniProt.MS", 
                          "fasta")
      
      
      dt <- expand_semicolon_rows(dt, columns_to_expand = cols_to_expand)
      
      # Commented out - grep strategy to include them without any collapsing
      # pattern <- paste(input$selected_isoforms, collapse="|")
      # dt <- dt[grepl(pattern, UniProt.MS)]
      dt <- dt[UniProt.MS %in% input$selected_isoforms]
    }
    
    # Return a data.table
    dt
  })
  
  #Plots -----
  
  output$summary_plot <- renderPlot({
    fd <- filtered_data()
    req(fd)
    
    ggplot(fd, aes(x = correlation, y = jenks_class, color = jenks_class)) +
      geom_violin(alpha = 0.5) +
      geom_jitter(height = 0.3, alpha = 0.7) +
      scale_color_manual(values = categorical_colors) +
      labs(x = "Peptide-Olink correlation", y = "Correlation category",
           title = "Peptide-Olink correlations by category") +
      theme_classic2() +
      theme(legend.position = 'none',
            text = element_text(size = 14, family = "Open Sans"),    
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 13),
            plot.title = element_text(size = 18)
    )
    
  })

  output$volcano_plot <- renderPlotly({
    fd <-  gene_stats %>% filter(!is.na(sd_corr), gene_symbol %in% filtered_data()$gene_symbol)
    req(fd)
    
    p <- ggplot(fd, aes(
      x = mean_corr,
      y = sd_corr,
      size = n_peptides,
      text = paste("Gene:", gene_symbol,
                   "<br>Mean Corr:", round(mean_corr, 2),
                   "<br>SD:", round(sd_corr, 2),
                   "<br>Peptides:", n_peptides)
    )) +
      geom_point(alpha = 0.8, fill = "#77D9C7", shape = 21, stroke = 0.2, color = 'grey40') +
      labs(x = "Mean correlation", y = "SD", title = "Mean vs. SD") +
      theme_classic2() +
      theme(plot.title = element_text(face = "bold"),
            legend.title = element_blank())
    
    ggplotly(p, tooltip = "text")
    
  })
  
  output$volcano_plot2 <- renderPlotly({
    fd <-  gene_stats %>% filter(!is.na(iqr_corr), gene_symbol %in% filtered_data()$gene_symbol)
    req(fd)
    
    p <- ggplot(fd, aes(
      x = median_corr,
      y = iqr_corr,
      size = n_peptides,
      text = paste("Gene:", gene_symbol,
                   "<br>Median Corr:", round(median_corr, 2),
                   "<br>IQR:", round(iqr_corr, 2),
                   "<br>Peptides:", n_peptides)
    )) +
      geom_point(alpha = 0.8, fill = "#77D9C7", shape = 21, stroke = 0.2, color = 'grey40') +
      labs(x = "Median correlation", y = "IQR", title = "Median vs. IQR") +
      theme_classic2() +
      theme(plot.title = element_text(face = "bold"),
            legend.title = element_blank())
    
    ggplotly(p, tooltip = "text")
    
  })
  

    
  # Check for the plot being possible using a reactive
  ngl_plot_obj <- reactive({
    req(input$selected_result)
    req(input$selected_isoforms)
    req(nrow(working_plasma_dt()) > 0)
    result <- try(
      alphafold_plot(
        ioi = input$selected_isoforms, 
        genesymb = input$selected_result, 
        uniprot_ids = working_plasma_dt(), 
        peptide_seq_list_file = "data/peptide_seq_list.RDS", 
        fasta_file = NULL, 
        correlation_palette = correlation_palette
      ),
      silent = FALSE
    )
    if (inherits(result, "try-error")) return(NULL)
    result
    
  })
  
  # Render the NGL plot output
  output$NGL_plot <- renderNGLVieweR({
    need(ngl_plot_obj(), message = paste("No valid AlphaFold structure file returned for this ID of interest."))
    ngl_plot_obj()
  })
  output$alphafold_warn <- renderUI({
    validate(
      need(ngl_plot_obj(), "No valid AlphaFold structure file returned for this ID of interest.")
    )
    NULL  })
  
  # Render the color scale only if the NGL plot is produced.
  output$color_scale <- renderPlot({
    if (is.null(ngl_plot_obj())) return(NULL)
    color_scale_plot(correlation_palette = sapply(seq(-1, 1, length.out = 100), correlation_palette_function),
                     n_bins = 100)
  })

  
  output$detailed_plot <- renderPlotly({
    req(input$selected_result)
    req(input$selected_isoforms)
    req(nrow(working_plasma_dt()) > 0)
    # # Use a reactive expression for selected_plasma_data
    # selected_plasma <- plasma_dt[.(input$selected_result)]
    # 
    # interpro_plot(
    #   ioi = input$selected_result,
    #   interpro_file = "data/interpro_domains.RDS",
    #   uniprot_ids = selected_plasma,  # Pass the correct subset
    #   peptide_seq_list_file = "data/peptide_seq_list.RDS",
    #   correlation_palette = correlation_palette,
    #   interpro_colors = c("#CCDDAA", "#EEEEBB", "#FFCCCC", "#DDDDDD", "#BBCCEE", "#CCEEFF")
    #)
    
    # Call the combined plot function
    combined_interpro_density_plot(
      ioi = input$selected_result,
      interpro_file = "data/interpro_domains.RDS",
      uniprot_ids = working_plasma_dt(),  
      peptide_seq_list_file = "data/peptide_seq_list.RDS",
      fasta_file = working_plasma_dt(),
      abundance_file = NULL, 
      abundance_column = "quant",
      correlation_palette =  sapply(seq(-1, 1, length.out = 100), correlation_palette_function),
      interpro_colors = c("#CCDDAA", "#EEEEBB", "#FFCCCC", "#DDDDDD", "#BBCCEE", "#CCEEFF")
    ) |> 
      config(
        # add download button
        toImageButtonOptions = list(
          format = "svg",
          filename = paste0('heatmap_', input$selected_result),
          width = 650,
          height = 300
        )
      )
  })
  
  addPopover(session, "summary_info", 
             title = NULL, 
             content = "View overall distributions and volcano plot. Adjust filters in the sidebar to refine.",
             placement = "right", 
             trigger = "click")
  
  addPopover(session, "detailed_info", 
             title = NULL, 
             content = "Explore gene-specific positional data (density & InterPro domains) and classification-based comparisons.",
             placement = "right", 
             trigger = "click")
}


shinyApp(ui = ui, server = server)
