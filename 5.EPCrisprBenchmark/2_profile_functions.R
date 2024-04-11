# Report time execution
report_time_execution <- function(fun){
  #' Report time execution
  #'
  start_time <- Sys.time()
  output <- fun
  print(Sys.time() - start_time)
  return(output)
}

# RSE to Granges
prepRSEtoGranges <- function(rse, assay="counts", coln_assay=1){
  gr = rowRanges(rse)
  gr$score = assay(rse, inputAssay=assay)[,coln_assay]
  return(gr)
}

# Return score if the position is present, 0 otherwise
get_count_vector <- function(pos, df){
  #' Return score if the position is present, 0 otherwise
  #' 
  if (pos %in% df$atac_relative_pos) {
    return(df$score[which(df$atac_relative_pos == pos)])
  } else {return(0)}
}

# Return score vector for each position for both strands
get_count_vector_both_strands <- function(df){
  #' Return score vector for each position for both strands
  #' 
  plus_count_vector <- sapply(1:len_vec, get_count_vector, df = df[df$strand == "+",])
  minus_count_vector <- sapply(1:len_vec, get_count_vector, df = df[df$strand == "-",])
  return(c(plus_count_vector, minus_count_vector))
}

get_chr_windows_profiles_nofilter <- function(cage_granges, atac_granges, chr, lib, setname){
  print(paste0("Performing extraction on ", chr, ".."))                                                                                           
  
  # Select chromosome
  cage_granges <- cage_granges[seqnames(cage_granges) == chr]
  atac_granges <- atac_granges[seqnames(atac_granges) == chr]

  # Add all information into one df
  overlaps <- findOverlaps(cage_granges, atac_granges)                   # Index of overlapping CAGE fragment
  
  # Check if there are ATAC positive windows overlapping CAGE data
  if (length(overlaps) > 0){  
    df <- cage_granges[queryHits(overlaps)]                              # Keep only overlapping CAGE data
    df$index_overlapping_atac <- subjectHits(overlaps)                   # Add index of overlapping ATAC   
    df %>% data.frame() %>%                                              # Add ATAC start site and relative position
      mutate(atac_start = start(atac_granges[subjectHits(overlaps)]),
             atac_relative_pos = pos - atac_start + 1,
             atac_txType = atac_granges$txType[subjectHits(overlaps)]) -> df
    print("#################")

    # Extract profiles of each (CAGE) overlapping ATAC region
    profiles <- by(data = df, 
                   INDICES = df$index_overlapping_atac, 
                   FUN = function(x){get_count_vector_both_strands(x)})
    profiles <- data.frame(do.call(rbind, profiles))
    
    # Add additional profiles which were not overlapped with CAGE at all
    vec <- c(1:length(atac_granges))
    vec1 <- vec[!vec %in% subjectHits(overlaps)] 
    profiles[as.character(vec1),] <- 0
    profiles <- profiles[order(as.numeric(row.names(profiles))), ]
    print(paste("additional profiles", vec1))
    
    # Add column names
    colnames(profiles) <- c(paste("Plus_", 1:len_vec, sep = ""), paste("Minus_", 1:len_vec, sep = ""))
    
    # Add metadata information
    profiles <- profiles %>% mutate(atac_start = start(atac_granges[as.numeric(rownames(profiles))]),
                                    chr = chr,
                                    atac_txType = atac_granges$txType[as.numeric(rownames(profiles))],
                                    lib = lib, 
                                    setname = setname) %>% relocate(c(chr, atac_start, atac_txType), .before = Plus_1) 
    return(profiles)
  } 
  else {
    print(paste(chr, "contains no overlapping windows"))
    return(NULL)
  }
}

# Return the windows profiles for all chromosomes
get_windows_profiles_nofilter <- function(cage_granges, atac_granges, lib_name, set_name){
  #' Return the windows profiles for all chromosomes
  #' 
  chromosomes <- unique(seqnames(cage_granges))
  list_chr_profiles <- lapply(chromosomes, function(x) 
  {report_time_execution(get_chr_windows_profiles_nofilter(cage_granges = cage_granges,                                                           
                                                           atac_granges = atac_granges,
                                                           chr = x, 
                                                           lib = lib_name,
                                                           setname = set_name))})
  output_list <- list()
  print("Concatenating all extracted profiles")                                                                                             
  output_list$profiles <- data.frame(dplyr::bind_rows(list_chr_profiles))                                           
  output_list$metadata <- output_list$profiles %>% dplyr::select(chr, atac_start, atac_txType, lib, setname)
  output_list$profiles <- output_list$profiles %>% dplyr::select(-chr, -atac_start, -atac_txType, -lib, -setname)
  return(output_list)
}

#' # Remove empty profiles
#' filter_empty_profiles <- function(windows){
#'   #' Remove empty profiles
#'   #' 
#'   non_overlapping_ranges <- apply(windows$profiles, 1, sum) == 0
#'   if (sum(non_overlapping_ranges) > 0){
#'     print(paste("Removing", sum(non_overlapping_ranges), "non-overlapping ranges"))
#'     windows$profiles <- windows$profiles[!non_overlapping_ranges,]
#'     windows$metadata <- windows$metadata[!non_overlapping_ranges,]}
#'   return(windows)
#' }

#' # Filter the profiles by CAGE signal 
#' windows_profiles_filter <- function(windows,
#'                                     threshold = 1, fun = max){
#'   #' Filter the profiles by CAGE signal 
#'   #' 
#'   filter <- apply(windows$profiles, 1, fun) > threshold
#'   windows$profiles <-  windows$profiles[filter,]
#'   windows$metadata <- windows$metadata[filter,]
#'   return(windows)
#' }

# Apply forward minus reverse subtraction 
strands_norm_subtraction <- function(vec){
  #' Apply forward minus reverse subtraction 
  #' 
  # Subtract the forward and reverse signal
  p <- as.numeric(vec[1:len_vec])
  m <- as.numeric(vec[(len_vec+1):(len_vec*2)])
  # Normalized strand subtraction
  return ((p - m) / max(abs(c(p, m))))
}

strands_norm_subtraction_all_windows <- function(windows){
  #' Apply normalized forward and reverse subtraction to all windows
  #' 
  p_min_m_norm_df <- as.data.frame(t(apply(windows, 1, strands_norm_subtraction)))
  pos <- seq(1, ncol(p_min_m_norm_df))
  colnames(p_min_m_norm_df) <- paste0("Pos", pos - ATAC_BP_EXT - 1)
  return(p_min_m_norm_df)
}

profile_nofilter = function(cage_rse, atac_gr, col=1, setDir="te_pos", setName="te.pos"){
  
  cage_gr = prepRSEtoGranges(cage_rse, assay="counts", coln_assay=col)
  
  writeLines("\nExtracting profiles:")
  windows <- report_time_execution(get_windows_profiles_nofilter(cage_gr, atac_gr, as.character(col), setName))
  print(paste("Extracted profiles:", nrow(windows$metadata)))
  
  # Filter by minimal CAGE requirement (TSS with atleast 2 [> threshold])
  #windows <- windows_profiles_filter(windows, threshold = 1)
  #print(paste("Extracted profiles after CAGE filtering:", nrow(windows$metadata)))
  
  ### Process the profiles
  writeLines("\nProcessing extracted profiles..")
  windows_subt <- list()
  windows_subt$profiles <- report_time_execution(strands_norm_subtraction_all_windows(windows$profiles))
  windows_subt$metadata <- windows$metadata
  
  writeLines("Exporting extracted profiles..")
  write_csv(windows$profiles, paste0(profiles_path, "/", setDir, "/profiles_", setName, ".", as.character(col), ".csv"))
  write_csv(windows$metadata, paste0(metadata_path, "/", setDir, "/metadata_", setName, ".", as.character(col), ".csv"))
  
  writeLines("Exporting processed extracted profiles..")
  write_csv(windows_subt$profiles, paste0(profiles_subtnorm_path, "/", setDir, "/profiles_", setName, ".", as.character(col), ".csv"))
  write_csv(windows_subt$metadata, paste0(metadata_subtnorm_path, "/", setDir, "/metadata_", setName, ".", as.character(col), ".csv"))
  
}


