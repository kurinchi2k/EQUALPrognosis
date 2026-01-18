process_data <- function(data_file_path, metadata_file_path) {
  get_date_format <- function(vector_of_dates) {
    vector_of_dates <- vector_of_dates[! is.na(vector_of_dates)]
    if (length(vector_of_dates) == 0) {
      output <- NA
    } else {
      separator <- ifelse(str_detect(vector_of_dates[1], "/"), "/", "-")
      vector_of_dates_split <- t(data.frame(
        str_split(vector_of_dates, separator)
      ))
      date_split <- lapply(1:3, function(x) {
        unique(vector_of_dates_split[,x])
      })
      year_position <- which(c(max(nchar(date_split[[1]])),
                               max(nchar(date_split[[2]])),
                               max(nchar(date_split[[3]]))
      ) == 4)
      remaining_fields <- setdiff(1:3,year_position)
      month_position <- remaining_fields[
        which(c(max(as.numeric(date_split[[remaining_fields[1]]])),
                max(as.numeric(date_split[[remaining_fields[2]]]))) <= 12)]
      if (length(month_position) == 2) {month_position <-
        if (year_position == 1) {
          2
        } else {
          if (str_detect(str_to_lower(Sys.getlocale("LC_TIME")), "_united states") |
              str_detect(str_to_lower(Sys.getlocale("LC_TIME")), "\\.united states") |
              str_detect(str_to_lower(Sys.getlocale("LC_TIME")), "_us\\.") |
              str_detect(str_to_lower(Sys.getlocale("LC_TIME")), "\\.us\\.") |
              str_detect(str_to_lower(Sys.getlocale("LC_TIME")), "\\.us$")) {
            1
          } else {
            2
          }
        }
      }
      date_position <- setdiff(1:3, c(year_position, month_position))
      format_text <- cbind.data.frame(
        text = c("%d", "%m", "%Y"),
        position = c(date_position, month_position, year_position)
      )
      format_text <- format_text[order(format_text$position),]
      format_text <- paste0(format_text$text, collapse = separator)
      output <- format_text
    }
    return(output)
  }
  mandatory_metadata_columns <- c("variable", "data_type")
  additional_metadata_columns <- c("reference", "ordinal_levels")
  accepted_data_types <- c("numerical", "count", "binary", "nominal", "ordinal", "date", "time")
  df <- try(read.csv(data_file_path, header = TRUE, check.names = FALSE, na.strings = c("NA", "", " ", "  ")), silent = TRUE)
  metadata <- try(read.csv(metadata_file_path, header = TRUE, check.names = FALSE, na.strings = c("NA", "", " ", "  ")), silent = TRUE)
  if (TRUE %in% c(class(df) == "try-error", class(metadata) == "try-error")) {
    output <- list(
      outcome = "Unsuccessful",
      message = "Either the data file or the metadata file has no columns",
      data_processed = NA,
      any_type = NA,
      quantitative = NA,
      numerical = NA,
      count = NA,
      binary = NA,
      nominal = NA,
      ordinal = NA,
      date = NA,
      time = NA
    )
  } else if (TRUE %in% c(nrow(df) == 0, nrow(metadata) == 0)) {
    output <- list(
      outcome = "Unsuccessful",
      message = "Either the data file or the metadata file has no data",
      data_processed = NA,
      any_type = NA,
      quantitative = NA,
      numerical = NA,
      count = NA,
      binary = NA,
      nominal = NA,
      ordinal = NA,
      date = NA,
      time = NA
    )
  } else {
    colnames(metadata) <- gsub("[^a-zA-Z0-9_ ()]","_",trimws(colnames(metadata)))
    if (TRUE %in% is.na(match(mandatory_metadata_columns, colnames(metadata)))) {
      output <- list(outcome = "Unsuccessful",
                     message = paste0("The metadata file does not contain the mandatory fields (columns). The metadata must contain the following columns: ",
                                      paste0("'", mandatory_metadata_columns, "'", collapse = ", "),
                                      ". Depending on the data_types, additional columns: ",
                                      paste0("'", additional_metadata_columns, "'", collapse = ", "),
                                      " may also be required. Please check and revise any column names that you may have misspelt."),
                     data_processed = NA,
                     any_type = NA,
                     quantitative = NA,
                     numerical = NA,
                     count = NA,
                     binary = NA,
                     nominal = NA,
                     ordinal = NA,
                     date = NA,
                     time = NA
      )
    } else {
      metadata <- metadata[(! is.na(metadata$variable)) & (metadata$variable %in% colnames(df)), ]
      if (TRUE %in% is.na(match(colnames(df), metadata$variable))) {
        missing_column_names = setdiff(colnames(df), metadata$variable)
        output <- list(outcome = "Unsuccessful",
                       message = paste0("The metadata file does not include the following columns: ",
                                        paste0("'", missing_column_names, "'", collapse = ", "),
                                        ". Please check whether you have missed or misspelt column names in the metadata file."),
                       data_processed = NA,
                       any_type = NA,
                       quantitative = NA,
                       numerical = NA,
                       count = NA,
                       binary = NA,
                       nominal = NA,
                       ordinal = NA,
                       date = NA,
                       time = NA
        )
      } else if (TRUE %in% is.na(metadata$data_type)) {
        missing_data_types = metadata$variable[is.na(metadata$data_type)]
        output <- list(outcome = "Unsuccessful",
                       message = paste0("The metadata file does not include the datatypes for the following columns: ",
                                        paste0("'", missing_data_types, "'", collapse = ", "),
                                        ". Please complete the data types for all columns."),
                       data_processed = NA,
                       any_type = NA,
                       quantitative = NA,
                       numerical = NA,
                       count = NA,
                       binary = NA,
                       nominal = NA,
                       ordinal = NA,
                       date = NA,
                       time = NA
        )
      } else if (TRUE %in% is.na(match(metadata$data_type, accepted_data_types))) {
        unrecognised_data_types <- setdiff(metadata$data_type, accepted_data_types)
        output <- list(outcome = "Unsuccessful",
                       message = paste0("The metadata file includes unrecognised data_types: ",
                                        paste0("'", unrecognised_data_types, "'", collapse = ", "),
                                        ". The recognised data_types are: ",
                                        paste0("'", accepted_data_types, "'", collapse = ", "),
                                        ". Please check whether you have misspelt data_types in the metadata file."),
                       data_processed = NA,
                       any_type = NA,
                       quantitative = NA,
                       numerical = NA,
                       count = NA,
                       binary = NA,
                       nominal = NA,
                       ordinal = NA,
                       date = NA,
                       time = NA
        )
      } else {
        metadata_outcome <- cbind.data.frame(variable = metadata$variable, outcome = NA)
        # First numerical
        numerical_columns <- metadata$variable[metadata$data_type %in% c("numerical","count")]
        if (length(numerical_columns) >0) {
          df[,numerical_columns] <- suppressWarnings(try(lapply(df[,numerical_columns], as.numeric), silent = TRUE))
          metadata_outcome$outcome[metadata_outcome$variable %in% numerical_columns] <- paste0("Successfully converted to ",  metadata$data_type[metadata_outcome$variable %in% numerical_columns], " data.")
        }
        # All non-numerical should be as characters
        df[,setdiff(colnames(df), numerical_columns)] <- sapply(df[,setdiff(colnames(df), numerical_columns)], as.character)
        # Next nominal
        nominal_columns <- metadata$variable[metadata$data_type %in% c("binary", "nominal")]
        # If there is no column called 'reference', then use the default method
        if (length(nominal_columns) > 0) {
          if (is.na(match("reference", colnames(metadata)))) {
            df[,nominal_columns] <- lapply(df[,nominal_columns], as.factor)
            metadata_outcome$outcome[metadata_outcome$variable %in% nominal_columns] <- paste0("Successfully converted to ",  metadata$data_type[metadata_outcome$variable %in% nominal_columns], " data. Since no reference was provided, the default reference by alphabetical order was used. If you think that you provided a column called 'reference', please check whether you have spelt this correctly")
          } else {
            metadata$reference <- as.character(metadata$reference)
            nominal_columns_no_reference <- metadata$variable[(metadata$data_type %in% c("binary", "nominal")) & (is.na(metadata$reference))]
            if (length(nominal_columns_no_reference) > 0) {
              df[,nominal_columns_no_reference] <- lapply(df[,nominal_columns_no_reference], as.factor)
              metadata_outcome$outcome[metadata_outcome$variable %in% nominal_columns_no_reference] <- paste0("Successfully converted to ",  metadata$data_type[metadata_outcome$variable %in% nominal_columns_no_reference], " data. Since no reference was provided, the default reference by alphabetical order was used. If you want a different reference category, please provide this in the correct row.")
            }
            nominal_columns_with_reference <- metadata$variable[(metadata$data_type %in% c("binary", "nominal")) & (! is.na(metadata$reference))]
            if (length(nominal_columns_with_reference) > 0) {
              nominal_columns_with_incorrect_reference <- nominal_columns_with_reference[unlist(lapply(nominal_columns_with_reference, function(x) {
                ! (metadata$reference[metadata$variable == x] %in% df[,x])}))]
              if (length(nominal_columns_with_incorrect_reference) > 0) {
                df[,nominal_columns_with_incorrect_reference] <- lapply(df[,nominal_columns_with_incorrect_reference], as.factor)
                metadata_outcome$outcome[metadata_outcome$variable %in% nominal_columns_with_incorrect_reference] <- paste0("Successfully converted to ",  metadata$data_type[metadata_outcome$variable %in% nominal_columns_with_incorrect_reference], " data. Since incorrect reference was provided, the default reference by alphabetical order was used. If you wanted a different reference category, please provide the correct reference category present in the data.")
              }
              nominal_columns_with_correct_reference <- setdiff(nominal_columns_with_reference, nominal_columns_with_incorrect_reference)
              if (length(nominal_columns_with_correct_reference) > 0) {
                df[,nominal_columns_with_correct_reference] <- lapply(nominal_columns_with_correct_reference, function(x) {
                  reference_category <- metadata$reference[metadata$variable == x]
                  category_levels <- c(reference_category, setdiff(unique(df[,x]), reference_category))
                  category_levels <- category_levels[! is.na(category_levels)]
                  output <- factor(df[,x], levels = category_levels, ordered = FALSE)
                })
                metadata_outcome$outcome[metadata_outcome$variable %in% nominal_columns_with_correct_reference] <- paste0("Successfully converted to ",  metadata$data_type[metadata_outcome$variable %in% nominal_columns_with_correct_reference], " data using the reference category provided.")
              }
            }
          }
        }
        # Next ordinal
        ordinal_columns <- metadata$variable[metadata$data_type == "ordinal"]
        if (length(ordinal_columns) > 0) {
          if (is.na(match("ordinal_levels", colnames(metadata)))) {
            df[,ordinal_columns] <- lapply(df[,ordinal_columns], function(x) {factor(x, ordered = TRUE)})
            metadata_outcome$outcome[metadata_outcome$variable %in% ordinal_columns] <- paste0("Successfully converted to ordinal data. Since no ordinal_levels were provided, the default order by alphabetical order was used. If you think that you provided a column called 'ordinal_levels', please check whether you have spelt this correctly")
          } else {
            metadata$ordinal_levels <- as.character(metadata$ordinal_levels)
            ordinal_columns_no_order <- metadata$variable[(metadata$data_type == "ordinal") & (is.na(metadata$ordinal_levels))]
            if (length(ordinal_columns_no_order) > 0) {
              df[,ordinal_columns_no_order] <- lapply(df[,ordinal_columns_no_order], function(x) {factor(x, ordered = TRUE)})
              metadata_outcome$outcome[metadata_outcome$variable %in% ordinal_columns_no_order] <- paste0("Successfully converted to ordinal data. Since no order was provided, the default order by alphabetical order was used. If you want a different order, please provide this in the correct row.")
            }
            ordinal_columns_with_order <- metadata$variable[(metadata$data_type == "ordinal") & (! is.na(metadata$ordinal_levels))]
            if (length(ordinal_columns_with_order) > 0) {
              ordinal_columns_with_incorrect_order <- ordinal_columns_with_order[unlist(lapply(ordinal_columns_with_order, function(x) {
                TRUE %in% (is.na(match(unique(df[,x])[! is.na(unique(df[,x]))], trimws(unlist(strsplit(metadata$ordinal_levels[metadata$variable == x], ";"))))))}))]
              if (length(ordinal_columns_with_incorrect_order) > 0) {
                df[,ordinal_columns_with_incorrect_order] <- lapply(df[,ordinal_columns_with_incorrect_order], function(x) {factor(x, ordered = TRUE)})
                metadata_outcome$outcome[metadata_outcome$variable %in% ordinal_columns_with_incorrect_order] <- paste0("Successfully converted to ordinal data. Since one or more levels were missing in the order provided, the default order by alphabetical order was used. If you wanted a different order, please provide the correct order present in the data.")
              }
              ordinal_columns_with_correct_order <- setdiff(ordinal_columns_with_order, ordinal_columns_with_incorrect_order)
              if (length(ordinal_columns_with_correct_order) > 0) {
                df[,ordinal_columns_with_correct_order] <- lapply(ordinal_columns_with_correct_order, function(x) {
                  category_levels <- trimws(unlist(strsplit(metadata$ordinal_levels[metadata$variable == x], ";")))
                  category_levels <- category_levels[! is.na(category_levels)]
                  output <- factor(df[,x], levels = category_levels, ordered = TRUE)
                })
                metadata_outcome$outcome[metadata_outcome$variable %in% ordinal_columns_with_correct_order] <- paste0("Successfully converted to ordinal data using the order provided.")
              }
            }
          }
        }
        output <- list(outcome = "Successful",
                       message = paste0("The data has been processed as follows.", "<br>",
                                        paste0("<b>", metadata_outcome$variable, ": ", "<\b>", metadata_outcome$outcome,  "<br>", collapse = "\n")),
                       data_processed = df,
                       any_type = metadata$variable,
                       quantitative = metadata$variable[metadata$data_type %in% c("numerical", "count")],
                       numerical = metadata$variable[metadata$data_type == "numerical"],
                       count = metadata$variable[metadata$data_type == "count"],
                       binary = metadata$variable[metadata$data_type == "binary"],
                       nominal = metadata$variable[metadata$data_type == "nominal"],
                       ordinal = metadata$variable[metadata$data_type == "ordinal"],
                       date = metadata$variable[metadata$data_type == "date"],
                       time = metadata$variable[metadata$data_type == "time"]
        )
      }
    }
  }
  return(output)
}
