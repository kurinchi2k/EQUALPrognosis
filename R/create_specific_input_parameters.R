create_specific_input_parameters <- function(generic_input_parameters, analysis_details_path, verbose = TRUE) {
  if (verbose == TRUE) {cat("\nGeneric input parameters\n")}
  # Check validity of generic_input_parameters
  generic_input_parameters_names <- c('general_title', 'simulations', 'simulations_per_file', 'seed', 'df', 'outcome_name', 'outcome_type', 'outcome_time', 'outcome_count')
  if (TRUE %in% (is.na(match(generic_input_parameters_names, names(generic_input_parameters))))) {
    outcome <- "Unsuccessful"
    message <- "The object provided for the input parameter 'generic_input_parameters' did not have the correct elements. Please create the generic_input_parameters using the R documentation for that function."
    specific_input_parameters <- NULL
  } else if (is.null(
    create_generic_input_parameters(
      general_title = generic_input_parameters$general_title, simulations = generic_input_parameters$simulations, simulations_per_file = generic_input_parameters$simulations_per_file, seed = generic_input_parameters$seed,
      df = generic_input_parameters$df, outcome_name = generic_input_parameters$outcome_name, outcome_type = generic_input_parameters$outcome_type, outcome_time = generic_input_parameters$outcome_time, outcome_count = generic_input_parameters$outcome_count, verbose = FALSE)$generic_input_parameters)) {
    outcome <- "Unsuccessful"
    message <- "The object provided for the input parameter 'generic_input_parameters' had errors. Please create the generic_input_parameters using the R documentation for that function."
    specific_input_parameters <- NULL
  } else {
    default_column_names <- c('name', 'analysis_title', 'develop_model', 'predetermined_model_text', 'mandatory_predictors', 'optional_predictors', 'mandatory_interactions', 'optional_interactions', 'model_threshold_method', 'scoring_system', 'predetermined_threshold', 'higher_values_event')
    analysis_details <- try(read.csv(analysis_details_path, header = TRUE, check.names = FALSE, na.strings = c("NA", "", " ", "  ")), silent = TRUE)
    if (verbose == TRUE) {cat("\n\nSpecific input parameters\n")}
    if (TRUE %in% (class(analysis_details) == "try-error")) {
      outcome <- "Unsuccessful"
      message <- "The outcome was unsuccessful. The analysis details file has no columns. Please upload a file with data."
      specific_input_parameters <- NULL
    } else if (TRUE %in% (nrow(analysis_details) == 0)) {
      outcome <- "Unsuccessful"
      message <- "The outcome was unsuccessful. The analysis details file has no rows. Please upload a file with data."
      specific_input_parameters <- NULL
    } else if (TRUE %in% (is.na(match(default_column_names, colnames(analysis_details))))) {
      outcome <- "Unsuccessful"
      message <- paste0("The outcome was unsuccessful. The analysis details file does not have the mandatory columns. The mandatory columns are: ",
                        paste0("'", default_column_names, "'", collapse = ", "), ". Please check whether you have missed or misspelt some column names.")
      specific_input_parameters <- NULL
    } else {
      checks <- data.frame(matrix(nrow = nrow(analysis_details), ncol = 3))
      colnames(checks) <- c("row_number", "comments", "fatal_errors")
      checks$row_number <- 1:nrow(analysis_details)
      checks$comments <- ""
      checks$fatal_errors <- FALSE
      # Create a new file with processed details
      processed_analysis_details <- analysis_details
      # Fatal errors for all scenarios #####
      missing_names <- ((checks$fatal_errors == FALSE) & (is.na(processed_analysis_details$name)))
      checks$comments[missing_names] <- paste0(checks$comments[missing_names], "\n",
                                               "The name is missing. Please provide a short name.")
      checks$fatal_errors[missing_names] <- TRUE
      duplicate_names <- ((checks$fatal_errors == FALSE) & (duplicated(processed_analysis_details$name)))
      checks$comments[duplicate_names] <- paste0(checks$comments[duplicate_names], "\n",
                                                 "There is a previous row with the same name. Please provide unique names.")
      checks$fatal_errors[duplicate_names] <- TRUE
      processed_analysis_details$develop_model <- suppressWarnings(as.logical(toupper(processed_analysis_details$develop_model)))
      invalid_develop_model <- ((checks$fatal_errors == FALSE) & (is.na(processed_analysis_details$develop_model)))
      checks$comments[invalid_develop_model] <- paste0(checks$comments[invalid_develop_model], "\n",
                                                       "The field develop_model is either empty or contains unrecognised text. The recognised texts are TRUE or FALSE.")
      checks$fatal_errors[invalid_develop_model] <- TRUE
      # Scoring systems ####
      missing_scoring_systems <- ((checks$fatal_errors == FALSE) & (processed_analysis_details$develop_model == FALSE) & (is.na(processed_analysis_details$scoring_system)))
      checks$comments[missing_scoring_systems] <- paste0(checks$comments[missing_scoring_systems], "\n",
                                                         "When develop_model is FALSE, the name of the scoring system is mandatory.")
      checks$fatal_errors[missing_scoring_systems] <- TRUE
      invalid_scoring_system_name <- ((checks$fatal_errors == FALSE) & (processed_analysis_details$develop_model == FALSE) & (is.na(match(processed_analysis_details$scoring_system, colnames(generic_input_parameters$df)))))
      checks$comments[invalid_scoring_system_name] <- paste0(checks$comments[invalid_scoring_system_name], "\n",
                                                             "The name of the scoring system does not match with any column name in the dataset.")
      checks$fatal_errors[invalid_scoring_system_name] <- TRUE
      scoring_system_not_numeric <- unlist(lapply(1:nrow(processed_analysis_details), function(x) {
        if ((checks$fatal_errors[x] == FALSE) & (processed_analysis_details$develop_model[x] == FALSE)) {
          (! is.numeric(generic_input_parameters$df[,processed_analysis_details$scoring_system[x]]))
        } else {
          FALSE
        }
      }))
      checks$comments[scoring_system_not_numeric] <- paste0(checks$comments[scoring_system_not_numeric], "\n",
                                                            "The scoring system must refer to a numerical column.")
      checks$fatal_errors[scoring_system_not_numeric] <- TRUE
      # If the outcome is quantitative, there is no threshold to use
      if (generic_input_parameters$outcome_type != "quantitative") {
        # Scoring system threshold
        processed_analysis_details$predetermined_threshold <- suppressWarnings(as.numeric(processed_analysis_details$predetermined_threshold))
        scoring_system_threshold_missing <- ((checks$fatal_errors == FALSE) & (processed_analysis_details$develop_model == FALSE) & (is.na(processed_analysis_details$predetermined_threshold)))
        checks$comments[scoring_system_threshold_missing] <- paste0(checks$comments[scoring_system_threshold_missing], "\n",
                                                                    "The predetermined threshold for the scoring system is either missing or not a number.")
        checks$fatal_errors[scoring_system_threshold_missing] <- TRUE
        invalid_thresholds <- unlist(lapply(1:nrow(processed_analysis_details), function(x) {
          if ((checks$fatal_errors[x] == FALSE) & (processed_analysis_details$develop_model[x] == FALSE)) {
            ((processed_analysis_details$predetermined_threshold[x]) < min(generic_input_parameters$df[,processed_analysis_details$scoring_system[x]])) |
              ((processed_analysis_details$predetermined_threshold[x]) > max(generic_input_parameters$df[,processed_analysis_details$scoring_system[x]]))
          } else {
            FALSE
          }
        }))
        checks$comments[invalid_thresholds] <- paste0(checks$comments[invalid_thresholds], "\n",
                                                      "The predetermined threshold is less the minimum value of the data or more than the maximum value of the data.")
        checks$fatal_errors[invalid_thresholds] <- TRUE
        # Scoring system higher values = events
        processed_analysis_details$higher_values_event <- suppressWarnings(as.logical(toupper(processed_analysis_details$higher_values_event)))
        higher_values_event_missing <- ((checks$fatal_errors == FALSE) & (processed_analysis_details$develop_model == FALSE) & (is.na(processed_analysis_details$higher_values_event)))
        checks$comments[higher_values_event_missing] <- paste0(checks$comments[higher_values_event_missing], "\n",
                                                               "The higher_values_event for the scoring system is either missing or in the appropriate format. The accepted values are TRUE or FALSE.")
        checks$fatal_errors[higher_values_event_missing] <- TRUE
      }
      # Model development ####
      # Model predictors
      no_predictors <- (checks$fatal_errors == FALSE) & (processed_analysis_details$develop_model == TRUE) &
        (is.na(processed_analysis_details$mandatory_predictors)) & (is.na(processed_analysis_details$optional_predictors))
      checks$comments[no_predictors] <- paste0(checks$comments[no_predictors], "\n",
                                               "When develop_model is TRUE, at least one mandatory or optional predictor is required.")
      checks$fatal_errors[no_predictors] <- TRUE
      invalid_mandatory_predictors <- unlist(lapply(1:nrow(processed_analysis_details), function(x) {
        if ((checks$fatal_errors[x] == FALSE) & (processed_analysis_details$develop_model[x] == TRUE)) {
          if (! is.na(processed_analysis_details$mandatory_predictors[x])) {
            TRUE %in% (is.na(match(trimws(unlist(strsplit(processed_analysis_details$mandatory_predictors[x], ";"))), colnames(generic_input_parameters$df))))
          } else {
            FALSE
          }
        } else {
          FALSE
        }
      }))
      checks$comments[invalid_mandatory_predictors] <- paste0(checks$comments[invalid_mandatory_predictors], "\n",
                                                              "One or more mandatory predictors are invalid. Multiple predictors must be separated by a semicolon. Please check whether there is a typographical error in entering the variable names or while combining them. Please also check whether the analysis details file corresponds to the database provided in the generic_input_parameters.")
      checks$fatal_errors[invalid_mandatory_predictors] <- TRUE
      invalid_optional_predictors <- unlist(lapply(1:nrow(processed_analysis_details), function(x) {
        if ((checks$fatal_errors[x] == FALSE) & (processed_analysis_details$develop_model[x] == TRUE)) {
          if (! is.na(processed_analysis_details$optional_predictors[x])) {
            TRUE %in% (is.na(match(trimws(unlist(strsplit(processed_analysis_details$optional_predictors[x], ";"))), colnames(generic_input_parameters$df))))
          } else {
            FALSE
          }
        } else {
          FALSE
        }
      }))
      checks$comments[invalid_optional_predictors] <- paste0(checks$comments[invalid_optional_predictors], "\n",
                                                             "One or more optional predictors are invalid. Multiple predictors must be separated by a semicolon. Please check whether there is a typographical error in entering the variable names or while combining them. Please also check whether the analysis details file corresponds to the database provided in the generic_input_parameters.")
      checks$fatal_errors[invalid_optional_predictors] <- TRUE
      invalid_mandatory_interactions <- unlist(lapply(1:nrow(processed_analysis_details), function(x) {
        if ((checks$fatal_errors[x] == FALSE) & (processed_analysis_details$develop_model[x] == TRUE)) {
          if (! is.na(processed_analysis_details$mandatory_interactions[x])) {
            TRUE %in% (is.na(match(trimws(unlist(strsplit(processed_analysis_details$mandatory_interactions[x], ";"))), colnames(generic_input_parameters$df))))
          } else {
            FALSE
          }
        } else {
          FALSE
        }
      }))
      checks$comments[invalid_mandatory_interactions] <- paste0(checks$comments[invalid_mandatory_interactions], "\n",
                                                                "One or more mandatory interactions are invalid. Multiple predictors must be separated by a semicolon. Please check whether there is a typographical error in entering the variable names or while combining them. Please also check whether the analysis details file corresponds to the database provided in the generic_input_parameters.")
      checks$fatal_errors[invalid_mandatory_interactions] <- TRUE
      invalid_optional_interactions <- unlist(lapply(1:nrow(processed_analysis_details), function(x) {
        if ((checks$fatal_errors[x] == FALSE) & (processed_analysis_details$develop_model[x] == TRUE)) {
          if (! is.na(processed_analysis_details$optional_interactions[x])) {
            TRUE %in% (is.na(match(trimws(unlist(strsplit(processed_analysis_details$optional_interactions[x], ";"))), colnames(generic_input_parameters$df))))
          } else {
            FALSE
          }
        } else {
          FALSE
        }
      }))
      checks$comments[invalid_optional_interactions] <- paste0(checks$comments[invalid_optional_interactions], "\n",
                                                               "One or more optional interactions are invalid. Multiple predictors must be separated by a semicolon. Please check whether there is a typographical error in entering the variable names or while combining them. Please also check whether the analysis details file corresponds to the database provided in the generic_input_parameters.")
      checks$fatal_errors[invalid_optional_interactions] <- TRUE
      # Invalid threshold
      invalid_threshold_method <- ((checks$fatal_errors == FALSE) & (processed_analysis_details$develop_model == TRUE) &
                                     ((is.na(processed_analysis_details$model_threshold_method)) |
                                        (! processed_analysis_details$model_threshold_method %in% c("youden", "topleft", "heuristic"))))
      checks$comments[invalid_threshold_method] <- paste0(checks$comments[invalid_threshold_method], "\n",
                                                          "The threshold method is either missing or is invalid. Valid threshold methods are 'youden', 'topleft', 'heuristic'.")
      checks$fatal_errors[invalid_threshold_method] <- TRUE
      # Non fatal errors ####
      # If analysis_title is missing, use the name
      missing_analysis_title <- ((checks$fatal_errors == FALSE) & (is.na(processed_analysis_details$analysis_title)))
      checks$comments[missing_analysis_title] <- paste0(checks$comments[missing_analysis_title], "\n",
                                                        "Missing title has been obtained from the name.")
      processed_analysis_details$analysis_title[missing_analysis_title] <- processed_analysis_details$name[missing_analysis_title]
      # Correct invalid names
      processed_analysis_details$name <- make.names(analysis_details$name)
      checks$comments[(checks$fatal_error == FALSE) & (processed_analysis_details$name != analysis_details$name)] <-
        paste0(checks$comments[(checks$fatal_error == FALSE) & (processed_analysis_details$name != analysis_details$name)], "\n",
               "Invalid names of analyses have been revised.")
      # Remove optional predictors that are also mandatory predictors
      optional_predictors <- do.call(rbind.data.frame, lapply(1:nrow(processed_analysis_details), function(x) {
        if (checks$fatal_errors[x] == TRUE) {
          output <- cbind.data.frame(optional_predictors = processed_analysis_details$optional_predictors[x],
                                     comments = NA)
        } else if ((processed_analysis_details$develop_model[x] == FALSE) | is.na(processed_analysis_details$optional_predictors[x])) {
          output <- cbind.data.frame(optional_predictors = processed_analysis_details$optional_predictors[x],
                                     comments = NA)
        } else {
          mandatory_predictors <- trimws(unlist(strsplit(processed_analysis_details$mandatory_predictors[x], ";")))
          optional_predictors <- trimws(unlist(strsplit(processed_analysis_details$optional_predictors[x], ";")))
          # Mandatory predictors cannot also be in optional predictors
          if (!TRUE %in% is.na(mandatory_predictors)) {
            updated_optional_predictors <- optional_predictors[! optional_predictors %in% mandatory_predictors]
          } else {
            updated_optional_predictors <- optional_predictors
          }
          if (TRUE %in% (is.na(updated_optional_predictors))) {
            output <- cbind.data.frame(optional_predictors = NA,
                                       comments = "Optional predictors present in mandatory predictors were removed. After the revisions, there were no optional predictors.")
          } else {
            output <- cbind.data.frame(optional_predictors = paste0(updated_optional_predictors, collapse = ";"),
                                       comments = if (length(updated_optional_predictors) == length(optional_predictors)) {
                                         NA
                                       } else {
                                         comments = "Some optional predictors were removed because of overlap with mandatory predictors."
                                       }
            )
          }
        }
        return(output)
      }))
      processed_analysis_details$optional_predictors <- optional_predictors$optional_predictors
      checks$comments[! is.na(optional_predictors$comments)] <- paste0(checks$comments[! is.na(optional_predictors$comments)], "\n",
                                                                       optional_predictors$comments[! is.na(optional_predictors$comments)])
      # Check that the mandatory interactions and optional interactions are valid
      mandatory_interactions <- do.call(rbind.data.frame, lapply(1:nrow(processed_analysis_details), function(x) {
        if (checks$fatal_errors[x] == TRUE) {
          output <- cbind.data.frame(mandatory_interactions = processed_analysis_details$mandatory_interactions[x],
                                     comments = NA)
        } else if ((processed_analysis_details$develop_model[x] == FALSE) | is.na(processed_analysis_details$mandatory_predictors[x]) | is.na(processed_analysis_details$mandatory_interactions[x])) {
          output <- cbind.data.frame(mandatory_interactions = processed_analysis_details$mandatory_interactions[x],
                                     comments = NA)
        } else {
          mandatory_predictors <- trimws(unlist(strsplit(processed_analysis_details$mandatory_predictors[x], ";")))
          mandatory_interactions <- trimws(unlist(strsplit(processed_analysis_details$mandatory_interactions[x], ";")))
          # Mandatory_interactions can include only mandatory predictors
          updated_mandatory_interactions <- mandatory_interactions[mandatory_interactions %in% mandatory_predictors]
          # There must be minimum of two predictors for mandatory interactions
          updated_mandatory_interactions <- if (length(updated_mandatory_interactions) <= 1){NA} else {updated_mandatory_interactions}
          if (TRUE %in% (is.na(updated_mandatory_interactions))) {
            output <- cbind.data.frame(mandatory_interactions = NA,
                                       comments = "Mandatory interactions must be included in mandatory predictors and must have more than one mandatory predictors. Invalid mandatory interactions have been removed.")
          } else {
            output <- cbind.data.frame(mandatory_interactions = paste0(updated_mandatory_interactions, collapse = ";"),
                                       comments = if (length(updated_mandatory_interactions) == length(mandatory_interactions)) {
                                         NA
                                       } else {
                                         "Mandatory interactions must be included in mandatory predictors. Invalid mandatory interactions have been removed."
                                       }
            )
          }
        }
        return(output)
      }))
      processed_analysis_details$mandatory_interactions <- mandatory_interactions$mandatory_interactions
      checks$comments[! is.na(mandatory_interactions$comments)] <- paste0(checks$comments[! is.na(mandatory_interactions$comments)], "\n",
                                                                          mandatory_interactions$comments[! is.na(mandatory_interactions$comments)])
      optional_interactions <- do.call(rbind.data.frame, lapply(1:nrow(processed_analysis_details), function(x) {
        if (checks$fatal_errors[x] == TRUE) {
          output <- cbind.data.frame(optional_interactions = processed_analysis_details$optional_interactions[x],
                                     comments = NA)
        } else if ((processed_analysis_details$develop_model[x] == FALSE) | is.na(processed_analysis_details$optional_interactions[x])) {
          output <- cbind.data.frame(optional_interactions = processed_analysis_details$optional_interactions[x],
                                     comments = NA)
        } else {
          mandatory_predictors <- if (is.na(processed_analysis_details$mandatory_predictors[x])) {NA} else {trimws(unlist(strsplit(processed_analysis_details$mandatory_predictors[x], ";")))}
          mandatory_interactions <- if (is.na(processed_analysis_details$mandatory_interactions[x])) {NA} else {trimws(unlist(strsplit(processed_analysis_details$mandatory_interactions[x], ";")))}
          optional_predictors <- if(is.na(processed_analysis_details$optional_predictors[x])) {NA} else {trimws(unlist(strsplit(processed_analysis_details$optional_predictors[x], ";")))}
          optional_interactions <- trimws(unlist(strsplit(processed_analysis_details$optional_interactions[x], ";")))
          # Mandatory interactions cannot also be in optional interactions
          if (!TRUE %in% is.na(mandatory_interactions)) {
            updated_optional_interactions <- optional_interactions[! optional_interactions %in% mandatory_interactions]
          } else {
            updated_optional_interactions <- optional_interactions
          }
          # Optional_interactions can include either mandatory or optional predictors
          updated_optional_interactions <- updated_optional_interactions[updated_optional_interactions %in% c(mandatory_predictors, optional_predictors)]
          # There must be minimum of two predictors for optional interactions
          updated_optional_interactions <- if (length(updated_optional_interactions) <= 1){NA} else {updated_optional_interactions}
          if (TRUE %in% (is.na(updated_optional_interactions))) {
            output <- cbind.data.frame(optional_interactions = NA,
                                       comments = "Some optional interactions were removed because they cannot be part of mandatory interactions. There must be more than one predictor in the optional interaction. After the revisions, there were no optional interactions.")
          } else {
            output <- cbind.data.frame(optional_interactions = paste0(updated_optional_interactions, collapse = ";"),
                                       comments = if (length(updated_optional_interactions) == length(optional_interactions)) {
                                         NA
                                       } else {
                                         comments = "Some optional interactions were removed because they cannot be part of mandatory interactions."
                                       }
            )
          }
        }
        return(output)
      }))
      processed_analysis_details$optional_interactions <- optional_interactions$optional_interactions
      checks$comments[! is.na(optional_interactions$comments)] <- paste0(checks$comments[! is.na(optional_interactions$comments)], "\n",
                                                                         optional_interactions$comments[! is.na(optional_interactions$comments)])
      # Invalid rows ####
      invalid_rows <- processed_analysis_details[checks$fatal_errors == TRUE,]
      checks_invalid_rows <- checks[checks$fatal_errors == TRUE,]
      if (nrow(invalid_rows) > 0) {
        invalid_rows_comments <- paste0(paste0("Row ", checks_invalid_rows$row_number, ": "), checks_invalid_rows$comments, collapse = "\n")
      } else {
        invalid_rows_comments <- ""
      }
      # Valid rows ####
      valid_rows <- processed_analysis_details[checks$fatal_errors == FALSE,]
      checks_valid_rows <- checks[checks$fatal_errors == FALSE,]
      checks_valid_rows <- checks_valid_rows[checks_valid_rows$comments != "",]
      if (nrow(valid_rows) > 0) {
        if (nrow(checks_valid_rows) > 0) {
          valid_rows_comments <- paste0(paste0("Row ", checks_valid_rows$row_number, ": "), checks_valid_rows$comments, collapse = "\n")
        } else {
          valid_rows_comments <- ""
        }
        valid_rows[valid_rows$develop_model == FALSE, c('mandatory_predictors', 'optional_predictors', 'mandatory_interactions', 'optional_interactions', 'model_threshold_method')] <- NA
        if (generic_input_parameters$outcome_type == "quantitative") {
          valid_rows[valid_rows$develop_model == FALSE, c('predetermined_threshold', 'higher_values_event')] <- NA
        }
        valid_rows[valid_rows$develop_model == TRUE, c('scoring_system', 'predetermined_threshold', 'higher_values_event')] <- NA
        valid_rows$predetermined_model_text[valid_rows$predetermined_model_text == ""] <- NA
        valid_rows$mandatory_predictors[valid_rows$mandatory_predictors == ""] <- NA
        valid_rows$optional_predictors[valid_rows$optional_predictors == ""] <- NA
        valid_rows$mandatory_interactions[valid_rows$mandatory_interactions == ""] <- NA
        valid_rows$optional_interactions[valid_rows$optional_interactions == ""] <- NA
        specific_input_parameters <- lapply(1:nrow(valid_rows), function(x){
          list(
            analysis_title = valid_rows$analysis_title[x],
            develop_model = valid_rows$develop_model[x],
            predetermined_model_text = valid_rows$predetermined_model_text[x],
            mandatory_predictors = valid_rows$mandatory_predictors[x],
            optional_predictors = valid_rows$optional_predictors[x],
            mandatory_interactions = valid_rows$mandatory_interactions[x],
            optional_interactions = valid_rows$optional_interactions[x],
            model_threshold_method = valid_rows$model_threshold_method[x],
            scoring_system = valid_rows$scoring_system[x],
            predetermined_threshold = valid_rows$predetermined_threshold[x],
            higher_values_event = valid_rows$higher_values_event[x]
          )
        })
        names(specific_input_parameters) <- valid_rows$name
      } else {
        valid_rows_comments <- ""
        specific_input_parameters <- NULL
      }
      # Message ####
      message <- {paste0("Number of valid rows: ", nrow(valid_rows), "\n",
                         "Number of invalid rows: ", nrow(invalid_rows), "\n",
                         if(invalid_rows_comments != "") {
                           paste0(
                             "Reasons for rows being invalid: ", "\n",
                             invalid_rows_comments, "\n"
                           )
                         },
                         if (nrow(valid_rows) > 0){
                           if(valid_rows_comments != "") {
                             paste0(
                               "Modifications made to valid rows: ", "\n",
                               valid_rows_comments, "\n"
                             )
                           } else {
                             "No notable modifications were made to the valid rows."
                           }
                         }
      )}
    }
  }
  # Output ####
  if (verbose == TRUE) {cat(message)}
  output <- list(outcome = message, specific_input_parameters = specific_input_parameters)
  return(output)
}
