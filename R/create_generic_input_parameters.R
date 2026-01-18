create_generic_input_parameters <- function(general_title = "Enter the title here", simulations = 2000, simulations_per_file = 100, seed = 1, df, outcome_name, outcome_type, outcome_time, outcome_count) {
  # Create a data frame for each aspect
  checks <- {cbind.data.frame(
    aspect = c("general_title", "simulations", "simulations_per_file", "seed", "df", "outcome_name", "outcome_type", "outcome_time", "outcome_count"),
    check_passed = NA,
    feedback = NA
  )}
  # Non fatal errors ####
  # Some checks are always passed with comments
  checks$check_passed[
    checks$aspect %in% c("general_title", "simulations", "simulations_per_file", "seed")] <- TRUE
  # Add comments for these items which are passed with corrections
  if ((is.na(general_title)) | (general_title == "") | (! is.character(general_title))) {
    general_title <- "Prediction"
    checks$feedback[checks$aspect == "general_title"] <- "The title was either not supplied or provided in the incorrect format. Therefore, a default title of 'Prediction' was assigned."
  }
  if ((is.na(simulations)) | (! is.numeric(simulations)) | (simulations <= 1)) {
    simulations <- 2000
    checks$feedback[checks$aspect == "simulations"] <- "Simulations were not provided or provided in the incorrect format. Therefore, the default simulations of 2000 was assigned."
  }
  if ((is.na(simulations_per_file)) | (! is.numeric(simulations_per_file)) | (simulations_per_file <= 1)) {
    simulations_per_file <- 100
    checks$feedback[checks$aspect == "simulations_per_file"] <- "Simulations per file were not provided or provided in the incorrect format. Therefore, the default simulations per file of 100 was assigned."
  }
  if ((is.na(seed)) | (! is.numeric(seed)) | (seed < 1)) {
    seed <- 1
    checks$feedback[checks$aspect == "seed"] <- "Seed was not provided or or provided in the incorrect format. Therefore, the default seed of 1 was assigned."
  }
  # Fatal errors ####
  checks$check_passed[checks$aspect == "df"] <- (is.data.frame(df))
  if (checks$check_passed[checks$aspect == "df"] == FALSE) {
    checks$feedback[checks$aspect == "df"] <- "The object provided for the input parameter 'df' was not a data frame. Please check whether you have included any quotes. Please process the data using the process_data function and provide the processed_data$df as input for this function. Please examples in R documentation for this function and process_data function."
  } else {
    # Data frames stored as tables cause errors - So convert to data frames
    df <- data.frame(df)
    checks$check_passed[checks$aspect == "outcome_name"] <- ((is.character(outcome_name)) & (! is.na(outcome_name)) & (outcome_name %in% colnames(df)))
    if (checks$check_passed[checks$aspect == "outcome_name"] == FALSE) {
      checks$feedback[checks$aspect == "outcome_name"] <- "The outcome name was missing or if provided was not a character or was not a column name in the data frame provided. Please provide the correct outcome name."
    } else {
      known_outcome_type <- ((is.character(outcome_type)) & (! is.na(outcome_type)) & (outcome_type %in% c("binary", "time-to-event", "quantitative")))
      if (known_outcome_type == FALSE) {
        checks$feedback[checks$aspect == "outcome_type"] <- "The outcome type was missing or if provided was not a character or was not a recognised outcome type. The recognised outcome types are 'binary', 'time-to-event', or 'quantitative'. Please examples in R documentation for this function and process_data function."
      } else {
        match_outcome_type_field_type <- ((outcome_type == "quantitative") & (is.numeric(df[,outcome_name]))) | ((outcome_type != "quantitative") & (is.factor(df[,outcome_name])))
        if (match_outcome_type_field_type == FALSE) {
          checks$check_passed[checks$aspect == "outcome_type"] <- FALSE
          checks$feedback[checks$aspect == "outcome_type"] <- "There was a mismatch between the outcome type and the data. Only outcomes formatted as numbers are accepted for quantitative outcomes and those formatted as factors are accepted for binary and time-to-event outcomes. Please see examples in R documentation for this function and process_data function."
        } else {
          checks$check_passed[checks$aspect == "outcome_type"] <- TRUE
          if (outcome_type != "time-to-event") {
            checks$check_passed[checks$aspect == "outcome_time"] <- TRUE
            if (! is.na(outcome_time)) {
              checks$feedback[checks$aspect == "outcome_time"] <- "Since this was not a time-to-event outcome, the value has been changed to NA."
              outcome_time <- NA
            }
          } else {
            outcome_time_present <- ((is.character(outcome_time)) & (! is.na(outcome_time)) & (outcome_time %in% colnames(df)))
            if (outcome_time_present == FALSE) {
              checks$check_passed[checks$aspect == "outcome_time"] <- FALSE
              checks$feedback[checks$aspect == "outcome_time"] <- "The outcome time was missing or if provided was not a character or was not a column name in the data frame provided. Please provide the correct outcome name."
            } else {
              outcome_time_numeric <- (is.numeric(df[,outcome_time]))
              if (outcome_time_numeric == FALSE) {
                checks$check_passed[checks$aspect == "outcome_time"] <- FALSE
                checks$feedback[checks$aspect == "outcome_time"] <- "The outcome time was missing or if provided was not a character or was not a column name in the data frame provided. Please provide the correct outcome name."
              } else {
                checks$check_passed[checks$aspect == "outcome_time"] <- TRUE
              }
            }
          }
          if (outcome_type != "quantitative") {
            checks$check_passed[checks$aspect == "outcome_count"] <- TRUE
            if ((is.na(outcome_count)) | (! is.logical(outcome_count)) | (outcome_count != FALSE)) {
              outcome_count <- FALSE
              checks$feedback[checks$aspect == "outcome_count"] <- "The outcome_count was missing or incorrect information was provided. This was changed to FALSE as the outcome is not a quantitative outcome."
            }
          } else {
            if ((is.na(outcome_count)) | (! is.logical(outcome_count)) | (! outcome_count %in% c(TRUE, FALSE))) {
              checks$check_passed[checks$aspect == "outcome_count"] <- FALSE
              checks$feedback[checks$aspect == "outcome_count"] <- "The outcome_count was missing or was provided in the incorrect format. Please provide either TRUE or FALSE without quotes."
            } else {
              checks$check_passed[checks$aspect == "outcome_count"] <- TRUE
            }
          }
        }
      }
    }
  }
  fatal_errors_present <- (FALSE %in% checks$check_passed)
  if (fatal_errors_present) {
    outcome <- "Unsuccessful"
    fatal_errors <- checks[(! is.na(checks$check_passed)) & (checks$check_passed == FALSE),]
    message <- paste0("The outcome was not successful. The reasons for failure are: \n",
                      paste0(fatal_errors$aspect, ": ", fatal_errors$feedback, collapse = "\n"))
    generic_input_parameters <- NULL
  } else {
    outcome <- "Succesful"
    if (FALSE %in% (is.na(checks$feedback))) {
      non_fatal_errors <- checks[! is.na(checks$feedback),]
      message <- paste0("The input parameters have been checked. Some corrections were made. You can use the corrected generic_input_parameters. The corrections are shown below: \n",
                        paste0(non_fatal_errors$aspect, ": ", non_fatal_errors$feedback, collapse = "\n"))
    } else {
      message <- "The input parameters have been checked and are correct for the data frame provided."
    }
    generic_input_parameters <- list(
      general_title = general_title,
      simulations = simulations,
      simulations_per_file = simulations_per_file,
      seed = seed,
      df = df,
      outcome_name = outcome_name,
      outcome_type = outcome_type,
      outcome_time = outcome_time,
      outcome_count = outcome_count
    )
  }
  cat(message)
  # Output ####
  return(generic_input_parameters)
}
