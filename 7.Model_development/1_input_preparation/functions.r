#' This function calculates the midpoint positions, referred to as "thick" 
#' positions, for each row in a DataFrame containing "start" and "end" columns. 
calculate_thick_df <- function(df) {
  # Check if the start and end columns exist and are numeric
  assertthat::assert_that(exists("start", where = df) && exists("end", where = df), # nolint: line_length_linter.
                          msg = "The 'start' and 'end' columns must exist.")
  assertthat::assert_that(is.numeric(df$start) && is.numeric(df$end),
                          msg = "The 'start' and 'end' columns must be numeric.") # nolint: line_length_linter.

  # Calculate the thick position
  thick_df <- df %>%
    dplyr::mutate(thick = .data$start + round((.data$end - .data$start) / 2))

  return(thick_df)
}



#' Extends the genomic intervals from their center points, known as "thick"
extend_from_center_thick_df <- function(df,
                                        dis,
                                        keep_same_length = TRUE,
                                        group_exact_positions = TRUE) {
  # Check if the start and end columns exist and are numeric
  assertthat::assert_that(exists("start", where = df) && exists("end", where = df), # nolint: line_length_linter.
                          msg = "The 'start' and 'end' columns must exist.")
  assertthat::assert_that(is.numeric(df$start) && is.numeric(df$end),
                          msg = "The 'start' and 'end' columns must be numeric.") # nolint: line_length_linter.

  # Check if the thick column exists
  if ("thick" %in% colnames(df)) {
    # Check if the thick column is numeric
    if (!is.numeric(df$thick)) {
      df <- calculate_thick_df(df)
    }
  } else {
    df <- calculate_thick_df(df)
  }

  if (keep_same_length) {
    ext_df <- df %>%
      dplyr::mutate(start = .data$thick - dis,
                    end = .data$thick + dis,
                    start = dplyr::case_when(start < 1 ~ 1,
                                             TRUE ~ start),
                    end = dplyr::case_when(start == 1 ~ (dis * 2 + 1),
                                           TRUE ~ end),
                    length = end - start + 1)
    # Assert that all intervals have the same length
    assertthat::assert_that(length(unique(ext_df$length)) == 1,
                            msg = "The length of the intervals is not the same.") # nolint: line_length_linter.
    # Recalculate the thick position to match the adjusted end position
    ext_df <- calculate_thick_df(ext_df)
  } else {
    ext_df <- df %>%
      dplyr::mutate(start = .data$thick - dis,
                    end = .data$thick + dis,
                    start = dplyr::case_when(.data$start < 1 ~ 1,
                                             TRUE ~ start),
                    length = end - start + 1)
    print("Length was recalculated to match the new start/end")
  }

  if (group_exact_positions) {
    ext_df <- group_exact_positions_df(ext_df)
    print("Grouped by exact positions")
  }

  return(ext_df)
}



#' Convert Row Name Strands to No-strand
convert_rowname_to_nostrand <- function(strand_str) {
  strand_str <- gsub("[+-]$", "*", strand_str)
  return(strand_str)
}



#' Modify Profile Row Names
modify_profile_rownames <- function(profiles, count_profiles) {
  plus_converted <- convert_rowname_to_nostrand(rownames(count_profiles$`*`$`+`)) # nolint: line_length_linter.
  minus_converted <- convert_rowname_to_nostrand(rownames(count_profiles$`*`$`-`)) # nolint: line_length_linter.
  if (identical(sort(plus_converted), sort(minus_converted))) {
    message("Both lists are the same after strand conversion.")
    rownames(profiles) <- plus_converted
  } else {
    message("The lists are not the same after strand conversion.")
  }
  return(profiles)
}



#' Combine Plus and Minus Strand Profiles
combine_plus_minus_profiles <- function(count_profiles, len_vec) {
  combined_profiles <- cbind(data.frame(count_profiles$`*`$`+`),
                             data.frame(count_profiles$`*`$`-`))
  colnames(combined_profiles) <- c(paste("Plus_", 1:len_vec, sep = ""),
                                   paste("Minus_", 1:len_vec, sep = ""))
  combined_profiles <- modify_profile_rownames(combined_profiles,
                                               count_profiles)
  return(combined_profiles)
}



#' Extract Components from Row Names
extract_rowname_components <- function(row_names_cpn) {
  # Split by ':' to separate chromosome and rest
  parts <- strsplit(row_names_cpn, ":")[[1]]
  chr <- parts[1]
  # Split the rest by '-' and ';' to get start, end, and strand
  rest <- strsplit(parts[2], "[-;]")[[1]]
  start <- as.numeric(rest[1])
  end <- as.numeric(rest[2])
  strand <- rest[3]
  return(c(chr, start, end, strand))
}



#' Create GRanges Object from Row Names
create_granges_from_rownames <- function(row_names_str) {
  # Apply the function to each row name
  components <- t(sapply(row_names_str, extract_rowname_components))
  # Create GRanges object
  gr <- GenomicRanges::GRanges(seqnames = components[, 1],
                               ranges = IRanges::IRanges(start = as.numeric(components[, 2]), # nolint: line_length_linter.
                                                         end = as.numeric(components[, 3])),  # nolint: line_length_linter.
                               strand = components[, 4])
  return(gr)
}



#' Measure and Report Execution Time of a Function
report_time_execution <- function(fun) {
  start_time <- Sys.time()          # Capture the start time
  output <- fun                     # Execute the function and store the output
  print(Sys.time() - start_time)    # Calculate and print the elapsed time
  return(output)                    # Return the function's output
}



#' Convert SummarizedExperiment to GRanges with Assay Data
cast_rse_to_granges <- function(rse,
                                assay = "counts",
                                coln_assay = 1,
                                colname = "score") {
  # Extract the row ranges from the SummarizedExperiment object
  gr <- SummarizedExperiment::rowRanges(rse) %>% GenomicRanges::GRanges() # nolint: pipe_operator_linter

  # Extract the assay data
  assay_data <- SummarizedExperiment::assay(rse, assay)

  # Assign the assay data to the specified column name in the GRanges object
  GenomicRanges::mcols(gr)[[colname]] <- assay_data[, coln_assay]

  return(gr)
}
