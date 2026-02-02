calculate_actual_predicted <- function(prepared_datasets, outcome_name, outcome_type, outcome_time, outcome_count, develop_model = TRUE, predetermined_model_text, mandatory_predictors, optional_predictors, mandatory_interactions = NULL, optional_interactions = NULL, model_threshold_method = "youden", scoring_system = NA, predetermined_threshold = NA, higher_values_event = NA, each_simulation = 1, bootstrap_sample = TRUE, verbose = TRUE) {
  # Functions ####
  calculate_heuristic_threshold <- function(lp, actual) {
    mean_no_event <- try(mean(lp[actual == levels(as.factor(actual))[[1]]], na.rm = TRUE), silent = TRUE)
    mean_event <- try(mean(lp[actual == levels(as.factor(actual))[[2]]], na.rm = TRUE), silent = TRUE)
    if (TRUE %in% c(class(mean_event) == "try-error", class(mean_no_event) == "try-error")) {
      output <- list(threshold = NA, direction = NA)
      return(output)
    }
    direction <- ifelse(mean_no_event < mean_event, "<", ">")
    # For unweighted - the cut-off is guided by simply prevalence
    event_prevalence <- length(actual[actual == levels(as.factor(actual))[[2]]])/length(actual)
    cut_off_point <- (1 - event_prevalence)
    distance <- (abs(max(lp, na.rm = TRUE) - min(lp, na.rm = TRUE)))*cut_off_point
    # Whether the distance should be min(lp) or max(lp) is determined by the direction
    threshold <- ifelse(direction == "<", min(lp, na.rm = TRUE) + distance, max(lp, na.rm = TRUE) - distance)
    output <- list(threshold = threshold, direction = direction)
    return(output)
  }
  calculate_intercept_slope_adjusted_lp <- function(outcome_type, outcome_count, intercept_slope_adjustment_model, manual_lp) {
    coef_intercept_slope <- cbind.data.frame(
      coef_names = names(coef(intercept_slope_adjustment_model)),
      coef = coef(intercept_slope_adjustment_model)
    )
    coef_table <- manual_lp$coef_table
    new_data_model_matrix <- manual_lp$new_data_model_matrix
    calibration_slope <- coef_intercept_slope$coef[match('lp_training', coef_intercept_slope$coef_names)]
    # Intercept is added for all outcomes other than time-to-event outcomes
    if (outcome_type != "time-to-event") {
      coef_table$coef[match('(Intercept)', coef_table$names)] <- coef_table$coef[match('(Intercept)', coef_table$names)] +
        coef_intercept_slope$coef[match('(Intercept)', coef_table$names)]
    }
    # Multiply the coefficients with calibration slope
    coef_table$coef <- coef_table$coef * calibration_slope
    new_data <- manual_lp$new_data
    new_data_model_matrix <- manual_lp$new_data_model_matrix
    training_model_matrix <- manual_lp$training_model_matrix
    if (outcome_type == "time-to-event") {
      base_hazard <- manual_lp$base_hazard
      closest_time_row <- manual_lp$closest_time_row
      # Cumulative hazard of the first subject at the time point closest to the subject
      htz <- base_hazard[closest_time_row, 1]
      variable_values_first_subject <- cbind.data.frame(t(training_model_matrix[1,coef_table$names]))
      colnames(variable_values_first_subject) <- coef_table$names
      variable_values_first_subject[1:nrow(new_data),] <- variable_values_first_subject
      variable_values_new_data_minus_first_subject <- data.frame(new_data_model_matrix[,coef_table$names] -
                                                                   variable_values_first_subject[,coef_table$names])
      colnames(variable_values_new_data_minus_first_subject) <- coef_table$names
      each_term <- lapply(1:nrow(coef_table), function(x) {
        output <- coef_table$coef[x] * variable_values_new_data_minus_first_subject[,coef_table$names[x]]
        output[is.na(output)] <- 0
        return(output)
      })
      names(each_term) <- coef_table$names
      each_term <- do.call(cbind.data.frame, each_term)
      rxz <- exp(rowSums(each_term))
      htx <- htz * rxz
      lp_new_data_untransformed <- htx
      lp_new_data <- (1-exp(-htx))
    } else {
      base_hazard <- NULL
      closest_time_row <- NULL
      each_term <- lapply(1:nrow(coef_table), function(x) {
        output <- coef_table$coef[x] * new_data_model_matrix[,coef_table$names[x]]
        output[is.na(output)] <- 0
        return(output)
      })
      names(each_term) <- coef_table$names
      each_term <- do.call(cbind.data.frame, each_term)
      lp_new_data_untransformed <- rowSums(each_term)
      if (outcome_type == "binary") {
        lp_new_data <- plogis(lp_new_data_untransformed)
      } else if (outcome_count == TRUE) {
        lp_new_data <- exp(lp_new_data_untransformed)
      } else {
        lp_new_data <- lp_new_data_untransformed
      }
    }
    output <- list(new_data = new_data,
                   training_model_matrix = training_model_matrix,
                   new_data_model_matrix = new_data_model_matrix,
                   base_hazard = base_hazard, closest_time_row = closest_time_row,
                   coef_table = coef_table,
                   lp_new_data_untransformed = lp_new_data_untransformed,
                   lp_new_data = lp_new_data)
  }
  calculate_manual_lp <- function(outcome_type, outcome_count, new_data, regression_model, variables_in_model_summary, df_training_complete) {
    coef_table <- {cbind.data.frame(
      names = names(coef(regression_model)),
      coef = coef(regression_model)
    )}
    coef_table <- coef_table[! is.na(coef_table$coef),]
    coef_table <- coef_table[coef_table$names %in% c(variables_in_model_summary$names_in_coefficient_table, "(Intercept)"),]
    coef_table$interaction_term <- (nchar(coef_table$names) != nchar(gsub(":", "", coef_table$names)))
    if (nrow(coef_table) > 0) {
      training_model_frame <- model.frame(regression_model)
      training_model_matrix <- model.matrix(regression_model)
      variables_in_model_frame <- unique(variables_in_model_summary$name[match(coef_table$names[coef_table$interaction_term == FALSE],
                                                                               variables_in_model_summary$names_in_coefficient_table)])
      variables_in_model_frame <- variables_in_model_frame[! is.na(variables_in_model_frame)]
      variables_in_model_matrix <- coef_table$names[coef_table$interaction_term == FALSE]
      new_data_model_matrix <- do.call(cbind.data.frame, lapply(1:length(variables_in_model_frame), function(x) {
        unprocessed_data <- new_data[,variables_in_model_frame[x]]
        variable_type <- variables_in_model_summary$type[match(variables_in_model_frame[x], variables_in_model_summary$name)]
        if (variable_type == "numeric") {
          output <- cbind.data.frame(unprocessed_data)
          colnames(output) <- variables_in_model_summary$revised_name[match(variables_in_model_frame[x], variables_in_model_summary$name)]
        } else {
          # If there are unrecognised levels in the data compared to that used in training, this should be changed to NA as the prediction cannot be made
          part_processed_data <- factor(as.character(unprocessed_data), levels = levels(training_model_frame[,variables_in_model_frame[x]]), ordered = (variable_type == "ordinal"))
          look_up_table <- variables_in_model_summary[variables_in_model_summary$name %in% variables_in_model_frame[x],]
          matching_columns <- look_up_table$names_in_coefficient_table[look_up_table$names_in_coefficient_table %in% coef_table$names]
          matching_rows <- match(look_up_table$levels, training_model_frame[,variables_in_model_frame[x]])
          look_up_table <- cbind.data.frame(
            look_up_table,
            training_model_matrix[matching_rows, matching_columns]
          )
          colnames(look_up_table)[6:ncol(look_up_table)] <- matching_columns
          output <- data.frame(matrix(nrow = length(part_processed_data), ncol = length(matching_columns)))
          colnames(output) <- matching_columns
          for (i in 1:nrow(look_up_table)) {
            output[(!is.na(part_processed_data)) & (part_processed_data == look_up_table$levels[i]),] <- look_up_table[i,matching_columns]
          }
        }
        return(output)
      }))
      if (outcome_type != "time-to-event") {
        new_data_model_matrix <- cbind.data.frame(
          `(Intercept)` = 1,
          new_data_model_matrix
        )
      }
      # Include interactions if any
      new_data_model_matrix <- do.call(cbind.data.frame, lapply(1:nrow(coef_table), function(x) {
        interaction_terms <- unlist(strsplit(coef_table$names[x], ":"))
        output <- eval(parse(text = paste0("new_data_model_matrix[,'", interaction_terms,"']", collapse = " * ")))
      }))
      colnames(new_data_model_matrix) <- coef_table$names
      if (outcome_type == "time-to-event") {
        base_hazard <- basehaz(regression_model, newdata = df_training_complete)
        closest_time_row <- unlist(lapply(new_data[,outcome_time], function(x) {
          output <- which.max(base_hazard$time[base_hazard$time <= x])
          if (length(output) < 1) {output <- 1}
          return(output)
        }))
        # Cumulative hazard of the first subject at the time point closest to the subject
        htz <- base_hazard[closest_time_row, 1]
        variable_values_first_subject <- cbind.data.frame(t(training_model_matrix[1,coef_table$names]))
        colnames(variable_values_first_subject) <- coef_table$names
        variable_values_first_subject[1:nrow(new_data),] <- variable_values_first_subject
        variable_values_new_data_minus_first_subject <- data.frame(new_data_model_matrix[,coef_table$names] -
                                                                     variable_values_first_subject[,coef_table$names])
        colnames(variable_values_new_data_minus_first_subject) <- coef_table$names
        each_term <- lapply(1:nrow(coef_table), function(x) {
          output <- coef_table$coef[x] * variable_values_new_data_minus_first_subject[,coef_table$names[x]]
          output[is.na(output)] <- 0
          return(output)
        })
        names(each_term) <- coef_table$names
        each_term <- do.call(cbind.data.frame, each_term)
        rxz <- exp(rowSums(each_term))
        htx <- htz * rxz
        lp_new_data_untransformed <- htx
        lp_new_data <- (1-exp(-htx))
      } else {
        base_hazard <- NULL
        closest_time_row <- NULL
        each_term <- lapply(1:nrow(coef_table), function(x) {
          output <- coef_table$coef[x] * new_data_model_matrix[,coef_table$names[x]]
          output[is.na(output)] <- 0
          return(output)
        })
        names(each_term) <- coef_table$names
        each_term <- do.call(cbind.data.frame, each_term)
        lp_new_data_untransformed <- rowSums(each_term)
        if (outcome_type == "binary") {
          lp_new_data <- plogis(lp_new_data_untransformed)
        } else if (outcome_count == TRUE) {
          lp_new_data <- exp(lp_new_data_untransformed)
        } else {
          lp_new_data <- lp_new_data_untransformed
        }
      }
      new_data = data.frame(new_data[,variables_in_model_frame])
      colnames(new_data) <- variables_in_model_frame
      output <- list(new_data = new_data,
                     training_model_matrix = training_model_matrix,
                     new_data_model_matrix = new_data_model_matrix,
                     base_hazard = base_hazard, closest_time_row = closest_time_row,
                     coef_table = coef_table,
                     lp_new_data_untransformed = lp_new_data_untransformed,
                     lp_new_data = lp_new_data)
    } else {
      output <- list(new_data = new_data,
                     training_model_matrix = NA,
                     new_data_model_matrix = NA,
                     base_hazard = NA, closest_time_row = NA,
                     coef_table = coef_table,
                     lp_new_data_untransformed = NA,
                     lp_new_data = NA)
    }
  }
  calculate_pROC_threshold <- function(lp, actual) {
    AUROC <- try(pROC::roc(actual,lp,ci = TRUE, ci.alpha = 0.95, print.auc=TRUE, quiet = TRUE), silent = TRUE)
    if (TRUE %in% (class(AUROC) == "try-error")) {
      output <- list(threshold_youden = NA, threshold_topleft = NA, direction = NA)
    } else {
      direction <- AUROC$direction
      threshold_youden <- try(pROC::coords(AUROC, x= "best", best.method = "youden", ret = "threshold"), silent = TRUE)
      threshold_topleft <- try(pROC::coords(AUROC, x= "best", best.method = "closest.topleft", ret = "threshold"), silent = TRUE)
      output <- list(threshold_youden = if(length(unlist(threshold_youden)) > 1) {NA} else {threshold_youden} , threshold_topleft = if(length(unlist(threshold_topleft)) > 1) {NA} else {threshold_topleft}, direction = direction)
    }
    return(output)
  }
  create_intercept_slope_adjustment_text <- function(outcome_name, outcome_type, outcome_time, outcome_count) {
    if (outcome_type %in% c("binary", "time-to-event", "quantitative")) {
      if ((outcome_type == "time-to-event") & (is.na(outcome_time))) {
        output <- "outcome_time is a mandatory variable for time-to-event outcome_type."
      } else {
        output <- {paste0(
          ifelse(outcome_type == 'time-to-event',
                 'coxph',
                 'glm'
          ),
          '(',
          ifelse(outcome_type == 'time-to-event',
                 paste0('Surv(', '`', outcome_time, '`', ', ', '`', outcome_name, '`', ')'),
                 paste0('`', outcome_name, '`')
          ),
          ' ~ back_transformed_lp_training' ,
          ', data = df_training_complete_copy',
          if (outcome_type != 'time-to-event') {
            paste0(
              ', family = ',
              ifelse(outcome_type == 'binary',
                     'binomial',
                     ifelse(outcome_count == TRUE,
                            'poisson',
                            'gaussian'
                     )
              )
            )
          },
          ')'
        )}
      }
    } else {
      output <- "Unrecognised outcome_type. Accepted outcome_types are: 'binary', 'time-to-event', 'quantitative'"
    }
    return(output)
  }
  create_model_text <- function(outcome_name, outcome_type, outcome_time, outcome_count, mandatory_predictors, optional_predictors, mandatory_interactions = NULL, optional_interactions = NULL) {
    if (outcome_type %in% c("binary", "time-to-event", "quantitative")) {
      if ((outcome_type == "time-to-event") & (is.na(outcome_time))) {
        output <- "outcome_time is a mandatory variable for time-to-event outcome_type."
      } else {
        if (! is.null(mandatory_interactions)) {
          # Mandatory_interactions can include only mandatory predictors
          mandatory_interactions <- mandatory_interactions[mandatory_interactions %in% mandatory_predictors]
          # There must be minimum of two predictors for mandatory interactions
          mandatory_interactions <- if (length(mandatory_interactions) <= 1){NULL} else {mandatory_interactions}
        }
        if (! is.null(optional_interactions)) {
          # Mandatory interactions cannot also be in optional interactions
          if (! is.null(mandatory_interactions)) {
            optional_interactions <- optional_interactions[! optional_interactions %in% mandatory_interactions]
          }
          # Optional interactions can include either mandatory or optional predictors
          optional_interactions <- optional_interactions[optional_interactions %in% c(mandatory_predictors, optional_predictors)]
          # There must be minimum of two predictors for optional interactions
          optional_interactions <- if (length(optional_interactions) <= 1) {NULL} else {optional_interactions}
        }
        include_interactions <- ((! is.null(mandatory_interactions)) | (! is.null(optional_interactions)))
        if ((is.null(optional_predictors)) &
            (
              (length(mandatory_predictors) == 1) |
              (include_interactions == FALSE) |
              (is.null(optional_interactions))
            )
        ) {
          output <- {paste0(
            ifelse(outcome_type == 'time-to-event',
                   'coxph',
                   'glm'
            ),
            '(',
            ifelse(outcome_type == 'time-to-event',
                   paste0('Surv(', '`', outcome_time, '`', ', ', '`', outcome_name, '`', ')'),
                   paste0('`', outcome_name, '`')
            ),
            ' ~ ' ,
            paste0(paste0('`', mandatory_predictors, '`'), collapse = ' + '),
            if (! is.null(mandatory_interactions)) {
              paste0(' + ', ' (',  paste0(paste0('`', mandatory_interactions, '`'), collapse = ' * '), ')')
            },
            ', data = df_training_complete',
            if (outcome_type != 'time-to-event') {
              paste0(
                ', family = ',
                ifelse(outcome_type == 'binary',
                       'binomial',
                       ifelse(outcome_count == TRUE,
                              'poisson',
                              'gaussian'
                       )
                )
              )
            } else {
              ', x = TRUE, y = TRUE'
            },
            ')'
          )}
        } else {
          output <- {paste0(
            'step(',
            ifelse(outcome_type == 'time-to-event',
                   'coxph',
                   'glm'
            ),
            '(',
            ifelse(outcome_type == 'time-to-event',
                   paste0('Surv(', '`', outcome_time, '`', ', ', '`', outcome_name, '`', ')'),
                   paste0('`', outcome_name, '`')
            ),
            ' ~ ' ,
            paste0(paste0('`', c(mandatory_predictors, optional_predictors), '`'), collapse = ' + '),
            if (! is.null(mandatory_interactions)) {
              paste0(' + ', ' (',  paste0(paste0('`', mandatory_interactions, '`'), collapse = ' * '), ')')
            },
            if (! is.null(optional_interactions)) {
              paste0(' + ', ' (',  paste0(paste0('`', optional_interactions, '`'), collapse = ' * '), ')')
            },
            ', data = df_training_complete',
            if (outcome_type != 'time-to-event') {
              paste0(
                ', family = ',
                ifelse(outcome_type == 'binary',
                       'binomial',
                       ifelse(outcome_count == TRUE,
                              'poisson',
                              'gaussian'
                       )
                )
              )
            } else {
              ', x = TRUE, y = TRUE'
            },
            ')',
            ', scope=list(upper = ~ ',
            paste0(paste0('`', c(mandatory_predictors, optional_predictors), '`'), collapse = ' + '),
            if (! is.null(mandatory_interactions)) {
              paste0(' + ', ' (',  paste0(paste0('`', mandatory_interactions, '`'), collapse = ' * '), ')')
            },
            if (! is.null(optional_interactions)) {
              paste0(' + ', ' (',  paste0(paste0('`', optional_interactions, '`'), collapse = ' * '), ')')
            },
            ', lower = ~ ',
            ifelse(is.null(mandatory_predictors), 1,
                   paste0(
                     paste0(paste0('`', mandatory_predictors, '`'), collapse = ' + '),
                     if (! is.null(mandatory_interactions)) {
                       paste0(' + ', ' (',  paste0(paste0('`', mandatory_interactions, '`'), collapse = ' * '), ')')
                     }
                   )
            ),
            '), trace = FALSE',
            ')'
          )}
        }
      }
    } else {
      output <- "Unrecognised outcome_type. Accepted outcome_types are: 'binary', 'time-to-event', 'quantitative'"
    }
    return(output)
  }
  r_object_to_html_code <- function(r_object) {
    r_object_content <- capture.output(r_object)
    r_object_content <- gsub("&", "&amp;", r_object_content)
    r_object_content <- gsub("<", "&lt;", r_object_content)
    r_object_content <- gsub(">", "&gt;", r_object_content)
    r_object_content <- gsub("'", "&apos;", r_object_content)
    r_object_content <- gsub('"', "&quot;", r_object_content)
    r_object_content <- gsub("\\t", "&emsp;", r_object_content)
    r_object_content <- gsub(" ", "&nbsp;", r_object_content)
    output <- paste0('<p>', paste0(r_object_content, collapse = "<br>"), '</p>')
  }
  remove_dropped_predictors_from_models <- function(model_text, dropped_moderators) {
    updated_model_text <- model_text
    for (i in 1:length(dropped_moderators)) {
      updated_model_text <- gsub(dropped_moderators[i], 0, updated_model_text)
    }
    return(updated_model_text)
  }
  # Start analysis ####
  html_file <- list()
  if (
    (TRUE %in% c(is.na(outcome_name), is.na(outcome_type))) |
    ((develop_model == TRUE) & (is.na(mandatory_predictors) & is.na(optional_predictors))) |
    ((develop_model == FALSE) & ((outcome_type != "quantitative") & ((is.na(scoring_system)) | is.na(predetermined_threshold) | is.na(higher_values_event)))) |
    ((develop_model == FALSE) & ((outcome_type == "quantitative") & ((is.na(scoring_system)))))
  ) {
    html_file$outcome <- "One or more mandatory fields are missing. For model development, the mandatory fields are model_text and additional_fields_to_include. For preexisting scoring systems, the mandatory fields are scoring_sytem, predetermined_threshold or predetermined_threshold_direction for binary and time-to-event outcomes, and scoring_sytem and quantitative variable formula for quantitative outcomes."
    output <- {list(
      actual_training = actual_training,
      predicted_training = NA,
      predicted_training_calibration_adjusted = NA,
      predicted_training_adjusted_mandatory_predictors_only = NA,
      actual_only_validation = actual_only_validation,
      predicted_only_validation = NA,
      predicted_only_validation_calibration_adjusted = NA,
      predicted_only_validation_adjusted_mandatory_predictors_only = NA,
      actual_all_subjects = actual_all_subjects,
      predicted_all_subjects = NA,
      predicted_all_subjects_calibration_adjusted = NA,
      predicted_all_subjects_adjusted_mandatory_predictors_only = NA,
      lp_training = NA,
      lp_only_validation = NA,
      lp_all_subjects = NA,
      lp_training_calibration_adjusted = NA,
      lp_only_validation_calibration_adjusted = NA,
      lp_all_subjects_calibration_adjusted = NA,
      lp_training_adjusted_mandatory_predictors_only = NA,
      lp_only_validation_adjusted_mandatory_predictors_only = NA,
      lp_all_subjects_adjusted_mandatory_predictors_only = NA,
      time_training = NA,
      time_only_validation = NA,
      time_all_subjects = NA,
      regression_model = NA,
      df_all_subjects = df_all_subjects,
      html_file = html_file,
      outcome = "One or more mandatory fields are missing. For model development, the mandatory fields are model_text and additional_fields_to_include. For preexisting scoring systems, the mandatory fields are scoring_sytem, predetermined_threshold or predetermined_threshold_direction for binary and time-to-event outcomes, and scoring_sytem and quantitative variable formula for quantitative outcomes."
    )}
  } else {
    # Set-up ####
    if (bootstrap_sample == TRUE) {
      if(verbose == TRUE) {cat(paste0(each_simulation, "..."))}
      df_training <- prepared_datasets$df_training_list[[each_simulation]]
    } else {
      if(verbose == TRUE) {cat(paste0("Apparent performance..."))}
      # Instead of the training data set provide the complete data set as the training set.
      df_training <- prepared_datasets$df_all_subjects_list[[each_simulation]]
    }
    df_training_complete <- df_training # For later use
    df_training_complete_copy <- df_training # For later use
    df_only_validation <- prepared_datasets$df_only_validation_list[[each_simulation]]
    df_all_subjects <- prepared_datasets$df_all_subjects_list[[each_simulation]]
    actual_training <- df_training_complete_copy[,outcome_name]
    actual_only_validation <- df_only_validation[,outcome_name]
    actual_all_subjects <- df_all_subjects[,outcome_name]
    predicted_training <- rep(NA, nrow(df_training_complete_copy))
    predicted_only_validation <- rep(NA, nrow(df_only_validation))
    predicted_all_subjects <- rep(NA, nrow(df_all_subjects))
    # If this is a time-to-event outcome, the outcome must be converted to numeric
    if (outcome_type == "time-to-event") {
      df_training[,outcome_name] <- as.numeric(df_training[,outcome_name])
      df_training_complete[,outcome_name] <- as.numeric(df_training_complete[,outcome_name])
      df_training_complete_copy[,outcome_name] <- as.numeric(df_training_complete_copy[,outcome_name])
      df_only_validation[,outcome_name] <- as.numeric(df_only_validation[,outcome_name])
      df_all_subjects[,outcome_name] <- as.numeric(df_all_subjects[,outcome_name])
    }
    # For scoring systems ####
    if (develop_model == FALSE) {
      # For existing scoring systems, for binary or time-to-event outcomes, the classification is done based on a predetermined threshold and direction of the predetermined threshold; for quantitative outcomes, the predictions are simply the value of the scoring sytem
      if (outcome_type == "quantitative") {
        predicted_training <- df_training_complete_copy[,scoring_system]
        predicted_only_validation <- df_only_validation[,scoring_system]
        predicted_all_subjects <- df_all_subjects[,scoring_system]
      } else {
        if (higher_values_event == TRUE) {
          predicted_training[df_training_complete_copy[,scoring_system] <= predetermined_threshold] <- levels(actual_training)[1]
          predicted_training[df_training_complete_copy[,scoring_system] > predetermined_threshold] <- levels(actual_training)[2]
          predicted_only_validation[df_only_validation[,scoring_system] <= predetermined_threshold] <- levels(actual_training)[1]
          predicted_only_validation[df_only_validation[,scoring_system] > predetermined_threshold] <- levels(actual_training)[2]
          predicted_all_subjects[df_all_subjects[,scoring_system] <= predetermined_threshold] <- levels(actual_all_subjects)[1]
          predicted_all_subjects[df_all_subjects[,scoring_system] > predetermined_threshold] <- levels(actual_all_subjects)[2]
        } else {
          predicted_training[df_training_complete_copy[,scoring_system] <= predetermined_threshold] <- levels(actual_training)[2]
          predicted_training[df_training_complete_copy[,scoring_system] > predetermined_threshold] <- levels(actual_training)[1]
          predicted_only_validation[df_only_validation[,scoring_system] <= predetermined_threshold] <- levels(actual_training)[2]
          predicted_only_validation[df_only_validation[,scoring_system] > predetermined_threshold] <- levels(actual_training)[1]
          predicted_all_subjects[df_all_subjects[,scoring_system] <= predetermined_threshold] <- levels(actual_all_subjects)[2]
          predicted_all_subjects[df_all_subjects[,scoring_system] > predetermined_threshold] <- levels(actual_all_subjects)[1]
        }
        predicted_training <- factor(predicted_training, levels = levels(actual_training))
        predicted_only_validation <- factor(predicted_only_validation, levels = levels(actual_only_validation))
        predicted_all_subjects <- factor(predicted_all_subjects, levels = levels(actual_all_subjects))
      }
      html_file$regression_model <- {paste0(
        '<h3>Regression model</h3>', '\n',
        'Not applicable as a preexisting scoring system and cut-off was used.'
      )}
      output <- {list(
        actual_training = actual_training,
        predicted_training = predicted_training,
        predicted_training_calibration_adjusted = NA,
        predicted_training_adjusted_mandatory_predictors_only = NA,
        actual_only_validation = actual_only_validation,
        predicted_only_validation = predicted_only_validation,
        predicted_only_validation_calibration_adjusted = NA,
        predicted_only_validation_adjusted_mandatory_predictors_only = NA,
        actual_all_subjects = actual_all_subjects,
        predicted_all_subjects = predicted_all_subjects,
        predicted_all_subjects_calibration_adjusted = NA,
        predicted_all_subjects_adjusted_mandatory_predictors_only = NA,
        lp_training = NA,
        lp_only_validation = NA,
        lp_all_subjects = NA,
        lp_training_calibration_adjusted = NA,
        lp_only_validation_calibration_adjusted = NA,
        lp_all_subjects_calibration_adjusted = NA,
        lp_training_adjusted_mandatory_predictors_only = NA,
        lp_only_validation_adjusted_mandatory_predictors_only = NA,
        lp_all_subjects_adjusted_mandatory_predictors_only = NA,
        time_training = NA,
        time_only_validation = NA,
        time_all_subjects = NA,
        regression_model = NA,
        df_all_subjects = df_all_subjects,
        html_file = html_file,
        outcome = "Successful"
      )}
    }
    # For model development ####
    if (develop_model == TRUE) {
      # Set-up ####
      if (is.na(mandatory_predictors)) {
        mandatory_predictors <- NULL
      } else {
        mandatory_predictors <- trimws(unlist(strsplit(mandatory_predictors, ";")))
        mandatory_predictors <- mandatory_predictors[mandatory_predictors %in% colnames(df_training_complete_copy)]
      }
      if (is.na(optional_predictors)) {
        optional_predictors <- NULL
      } else {
        optional_predictors <- trimws(unlist(strsplit(optional_predictors, ";")))
        optional_predictors <- optional_predictors[optional_predictors %in% colnames(df_training_complete_copy)]
      }
      if (is.na(mandatory_interactions)) {
        mandatory_interactions <- NULL
      } else {
        mandatory_interactions <- trimws(unlist(strsplit(mandatory_interactions, ";")))
        mandatory_interactions <- mandatory_interactions[mandatory_interactions %in% colnames(df_training_complete_copy)]
      }
      if (is.na(optional_interactions)) {
        optional_interactions <- NULL
      } else {
        optional_interactions <- trimws(unlist(strsplit(optional_interactions, ";")))
        optional_interactions <- optional_interactions[optional_interactions %in% colnames(df_training_complete_copy)]
      }
      df_training_complete <- df_training_complete_copy[,c(outcome_name, if (outcome_type == "time-to-event") {outcome_time}, mandatory_predictors, optional_predictors)]
      df_training_complete <- na.omit(df_training_complete)
      # After removing some rows, some categorical predictors may have a single level and this can cause error
      levels_with_data <- sapply(1:ncol(df_training_complete), function(y) {length(table(df_training_complete[,y])[table(df_training_complete[,y]) > 0])})
      # Update the models (after removing these moderators with only one level)
      updated_moderators <- setdiff(colnames(df_training_complete), colnames(df_training_complete)[which(levels_with_data < 2)])
      dropped_moderators <- colnames(df_training_complete)[which(levels_with_data < 2)]
      mandatory_predictors <- mandatory_predictors[mandatory_predictors %in% updated_moderators]
      optional_predictors <- optional_predictors[optional_predictors %in% updated_moderators]
      # Regression model ####
      if ((! is.na(predetermined_model_text)) & (predetermined_model_text != "")) {
        model_text <- predetermined_model_text
        if(length(dropped_moderators) > 0) {model_text <- remove_dropped_predictors_from_models(model_text = model_text, dropped_moderators = dropped_moderators)}
      } else {
        model_text <- {create_model_text(
          outcome_name = outcome_name, outcome_type = outcome_type,
          outcome_time = outcome_time, outcome_count = outcome_count,
          mandatory_predictors = mandatory_predictors,
          optional_predictors = optional_predictors,
          mandatory_interactions = mandatory_interactions,
          optional_interactions = optional_interactions)}
      }
      regression_model <- suppressWarnings(try(eval(parse(text = model_text)), silent = TRUE))
      failed_simpler_model <- FALSE
      if (TRUE %in% (class(regression_model) == "try-error")) {
        # Try a simpler model
        model_text <- {create_model_text(
          outcome_name = outcome_name, outcome_type = outcome_type,
          outcome_time = outcome_time, outcome_count = outcome_count,
          mandatory_predictors = mandatory_predictors,
          optional_predictors = NULL,
          mandatory_interactions = NULL,
          optional_interactions = NULL)}
        regression_model <- suppressWarnings(try(eval(parse(text = model_text)), silent = TRUE))
        if (TRUE %in% (class(regression_model) == "try-error")) {
          failed_simpler_model <- TRUE
          output <- {list(
            actual_training = actual_training,
            predicted_training = NA,
            predicted_training_calibration_adjusted = NA,
            predicted_training_adjusted_mandatory_predictors_only = NA,
            actual_only_validation = actual_only_validation,
            predicted_only_validation = NA,
            predicted_only_validation_calibration_adjusted = NA,
            predicted_only_validation_adjusted_mandatory_predictors_only = NA,
            actual_all_subjects = actual_all_subjects,
            predicted_all_subjects = NA,
            predicted_all_subjects_calibration_adjusted = NA,
            predicted_all_subjects_adjusted_mandatory_predictors_only = NA,
            lp_training = NA,
            lp_only_validation = NA,
            lp_all_subjects = NA,
            lp_training_calibration_adjusted = NA,
            lp_only_validation_calibration_adjusted = NA,
            lp_all_subjects_calibration_adjusted = NA,
            lp_training_adjusted_mandatory_predictors_only = NA,
            lp_only_validation_adjusted_mandatory_predictors_only = NA,
            lp_all_subjects_adjusted_mandatory_predictors_only = NA,
            time_training = NA,
            time_only_validation = NA,
            time_all_subjects = NA,
            regression_model = regression_model,
            df_all_subjects = df_all_subjects,
            html_file = html_file,
            outcome = "There was an error in the model: try a simpler model."
          )}
        }
      }
      if (! TRUE %in% (class(regression_model) == "try-error")) {
        # Manual predictions ####
        variables_in_model <- setdiff(updated_moderators, c(outcome_name, if(outcome_type == "time-to-event"){outcome_time}))
        revised_names <- unlist(lapply(variables_in_model, function(x) {
          if (x == make.names(x)) {
            x
          } else {
            paste0("`", x, "`")
          }
        }))
        variables_in_model_summary <- do.call(rbind.data.frame, lapply(1:length(variables_in_model), function(x) {
          type <- class(df_training_complete_copy[,variables_in_model[x]])
          if (("numeric" %in% type) | ("integer" %in% type)) {
            output <- cbind.data.frame(
              name = variables_in_model[x],
              revised_name = revised_names[x],
              type = "numeric",
              levels = NA,
              names_in_coefficient_table = revised_names[x]
            )
          } else if ("ordered" %in% type) {
            categories <- levels(df_training_complete_copy[,variables_in_model[x]])
            output <- cbind.data.frame(
              name = variables_in_model[x],
              revised_name = revised_names[x],
              type = "ordinal",
              levels = categories,
              names_in_coefficient_table = c(NA, {paste0(
                revised_names[x],
                paste0(
                  c(
                    ".L",
                    if (length(categories) > 2) {".Q"},
                    if (length(categories) > 3) {".C"},
                    if (length(categories) > 4) {paste0("^", 4:(length(categories)-1))}
                  )
                )
              )})
            )
          } else {
            categories <- levels(df_training_complete_copy[,variables_in_model[x]])
            output <- cbind.data.frame(
              name = variables_in_model[x],
              revised_name = revised_names[x],
              type = "nominal",
              levels = categories,
              names_in_coefficient_table = {paste0(
                revised_names[x],
                categories
              )}
            )
          }
          return(output)
        }))
        manual_lp_training <- try(calculate_manual_lp(outcome_type = outcome_type, outcome_count = outcome_count,
                                                      new_data = df_training_complete_copy, regression_model = regression_model,
                                                      variables_in_model_summary = variables_in_model_summary, df_training_complete = df_training_complete), silent = TRUE)
        manual_lp_only_validation <- try(calculate_manual_lp(outcome_type = outcome_type, outcome_count = outcome_count,
                                                             new_data = df_only_validation, regression_model = regression_model,
                                                             variables_in_model_summary = variables_in_model_summary, df_training_complete = df_training_complete), silent = TRUE)
        manual_lp_all_subjects <- try(calculate_manual_lp(outcome_type = outcome_type, outcome_count = outcome_count,
                                                          new_data = df_all_subjects, regression_model = regression_model,
                                                          variables_in_model_summary = variables_in_model_summary, df_training_complete = df_training_complete), silent = TRUE)
        # Linear predictors directly from the model ####
        {
          if (outcome_type != "time-to-event") {
            lp_training_with_se <- try(predict(regression_model, newdata = df_training_complete_copy, type = "response", se.fit = TRUE), silent = TRUE)
            lp_only_validation_with_se <- try(predict(regression_model, newdata = df_only_validation, type = "response", se.fit = TRUE), silent = TRUE)
            lp_all_subjects_with_se <- try(predict(regression_model, newdata = df_all_subjects, type = "response", se.fit = TRUE), silent = TRUE)
          } else {
            lp_training_with_se <- try(predict(regression_model, newdata = df_training_complete_copy, type = "expected", se.fit = TRUE), silent = TRUE)
            lp_only_validation_with_se <- try(predict(regression_model, newdata = df_only_validation, type = "expected", se.fit = TRUE), silent = TRUE)
            lp_all_subjects_with_se <- try(predict(regression_model, newdata = df_all_subjects, type = "expected", se.fit = TRUE), silent = TRUE)
          }
        }
        # Error in linear prediction ####
        if ((TRUE %in% (class(manual_lp_training) == "try-error")) & (TRUE %in% (class(lp_training_with_se) == "try-error"))) {
          # Try a simpler model
          model_text <- {create_model_text(
            outcome_name = outcome_name, outcome_type = outcome_type,
            outcome_time = outcome_time, outcome_count = outcome_count,
            mandatory_predictors = mandatory_predictors,
            optional_predictors = NULL,
            mandatory_interactions = NULL,
            optional_interactions = NULL)}
          regression_model <- suppressWarnings(try(eval(parse(text = model_text)), silent = TRUE))
          if (TRUE %in% (class(regression_model) == "try-error")) {
            failed_simpler_model <- TRUE
            output <- {list(
              actual_training = actual_training,
              predicted_training = NA,
              predicted_training_calibration_adjusted = NA,
              predicted_training_adjusted_mandatory_predictors_only = NA,
              actual_only_validation = actual_only_validation,
              predicted_only_validation = NA,
              predicted_only_validation_calibration_adjusted = NA,
              predicted_only_validation_adjusted_mandatory_predictors_only = NA,
              actual_all_subjects = actual_all_subjects,
              predicted_all_subjects = NA,
              predicted_all_subjects_calibration_adjusted = NA,
              predicted_all_subjects_adjusted_mandatory_predictors_only = NA,
              lp_training = NA,
              lp_only_validation = NA,
              lp_all_subjects = NA,
              lp_training_calibration_adjusted = NA,
              lp_only_validation_calibration_adjusted = NA,
              lp_all_subjects_calibration_adjusted = NA,
              lp_training_adjusted_mandatory_predictors_only = NA,
              lp_only_validation_adjusted_mandatory_predictors_only = NA,
              lp_all_subjects_adjusted_mandatory_predictors_only = NA,
              time_training = NA,
              time_only_validation = NA,
              time_all_subjects = NA,
              regression_model = regression_model,
              df_all_subjects = df_all_subjects,
              html_file = html_file,
              outcome = "There was an error in the model: try a simpler model."
            )}
          } else {
            # Manual predictions ####
            variables_in_model <- setdiff(updated_moderators, c(outcome_name, if(outcome_type == "time-to-event"){outcome_time}))
            revised_names <- unlist(lapply(variables_in_model, function(x) {
              if (x == make.names(x)) {
                x
              } else {
                paste0("`", x, "`")
              }
            }))
            variables_in_model_summary <- do.call(rbind.data.frame, lapply(1:length(variables_in_model), function(x) {
              type <- class(df_training_complete_copy[,variables_in_model[x]])
              if (("numeric" %in% type) | ("integer" %in% type)) {
                output <- cbind.data.frame(
                  name = variables_in_model[x],
                  revised_name = revised_names[x],
                  type = "numeric",
                  levels = NA,
                  names_in_coefficient_table = revised_names[x]
                )
              } else if ("ordered" %in% type) {
                categories <- levels(df_training_complete_copy[,variables_in_model[x]])
                output <- cbind.data.frame(
                  name = variables_in_model[x],
                  revised_name = revised_names[x],
                  type = "ordinal",
                  levels = categories,
                  names_in_coefficient_table = c(NA, {paste0(
                    revised_names[x],
                    paste0(
                      c(
                        ".L",
                        if (length(categories) > 2) {".Q"},
                        if (length(categories) > 3) {".C"},
                        if (length(categories) > 4) {paste0("^", 4:(length(categories)-1))}
                      )
                    )
                  )})
                )
              } else {
                categories <- levels(df_training_complete_copy[,variables_in_model[x]])
                output <- cbind.data.frame(
                  name = variables_in_model[x],
                  revised_name = revised_names[x],
                  type = "nominal",
                  levels = categories,
                  names_in_coefficient_table = {paste0(
                    revised_names[x],
                    categories
                  )}
                )
              }
              return(output)
            }))
            manual_lp_training <- try(calculate_manual_lp(outcome_type = outcome_type, outcome_count = outcome_count,
                                                          new_data = df_training_complete_copy, regression_model = regression_model,
                                                          variables_in_model_summary = variables_in_model_summary, df_training_complete = df_training_complete), silent = TRUE)
            manual_lp_only_validation <- try(calculate_manual_lp(outcome_type = outcome_type, outcome_count = outcome_count,
                                                                 new_data = df_only_validation, regression_model = regression_model,
                                                                 variables_in_model_summary = variables_in_model_summary, df_training_complete = df_training_complete), silent = TRUE)
            manual_lp_all_subjects <- try(calculate_manual_lp(outcome_type = outcome_type, outcome_count = outcome_count,
                                                              new_data = df_all_subjects, regression_model = regression_model,
                                                              variables_in_model_summary = variables_in_model_summary, df_training_complete = df_training_complete), silent = TRUE)
            # Linear predictors directly from the model ####
            {
              if (outcome_type != "time-to-event") {
                lp_training_with_se <- try(predict(regression_model, newdata = df_training_complete_copy, type = "response", se.fit = TRUE), silent = TRUE)
                lp_only_validation_with_se <- try(predict(regression_model, newdata = df_only_validation, type = "response", se.fit = TRUE), silent = TRUE)
                lp_all_subjects_with_se <- try(predict(regression_model, newdata = df_all_subjects, type = "response", se.fit = TRUE), silent = TRUE)
              } else {
                lp_training_with_se <- try(predict(regression_model, newdata = df_training_complete_copy, type = "expected", se.fit = TRUE), silent = TRUE)
                lp_only_validation_with_se <- try(predict(regression_model, newdata = df_only_validation, type = "expected", se.fit = TRUE), silent = TRUE)
                lp_all_subjects_with_se <- try(predict(regression_model, newdata = df_all_subjects, type = "expected", se.fit = TRUE), silent = TRUE)
              }
            }

            # Error in linear prediction ####
            if ((TRUE %in% (class(manual_lp_training) == "try-error")) & (TRUE %in% (class(lp_training_with_se) == "try-error"))) {
              failed_simpler_model <- TRUE
              output <- {list(
                actual_training = actual_training,
                predicted_training = NA,
                predicted_training_calibration_adjusted = NA,
                predicted_training_adjusted_mandatory_predictors_only = NA,
                actual_only_validation = actual_only_validation,
                predicted_only_validation = NA,
                predicted_only_validation_calibration_adjusted = NA,
                predicted_only_validation_adjusted_mandatory_predictors_only = NA,
                actual_all_subjects = actual_all_subjects,
                predicted_all_subjects = NA,
                predicted_all_subjects_calibration_adjusted = NA,
                predicted_all_subjects_adjusted_mandatory_predictors_only = NA,
                lp_training = NA,
                lp_only_validation = NA,
                lp_all_subjects = NA,
                lp_training_calibration_adjusted = NA,
                lp_only_validation_calibration_adjusted = NA,
                lp_all_subjects_calibration_adjusted = NA,
                lp_training_adjusted_mandatory_predictors_only = NA,
                lp_only_validation_adjusted_mandatory_predictors_only = NA,
                lp_all_subjects_adjusted_mandatory_predictors_only = NA,
                time_training = NA,
                time_only_validation = NA,
                time_all_subjects = NA,
                regression_model = regression_model,
                df_all_subjects = df_all_subjects,
                html_file = html_file,
                outcome = "There was an error in linear prediction using the model: try a simpler model."
              )}
            }
          }
        }
      }
      if (failed_simpler_model == FALSE) {
        # Demonstrate the regression model (by getting to html_file) if this was apparent performance ####
        if (bootstrap_sample == FALSE){
          html_file$regression_model <- {paste0(
            '<h3>Regression model</h3>', '\n',
            r_object_to_html_code(r_object = regression_model)
          )}
        }
        # Further processing ####
        if (TRUE %in% c(class(lp_training_with_se) == "try-error", class(lp_only_validation_with_se) == "try-error", class(lp_all_subjects_with_se) == "try-error")) {
          if (! TRUE %in% (class(manual_lp_training) == "try-error")) {lp_training <- manual_lp_training$lp_new_data}
          if (! TRUE %in% (class(manual_lp_only_validation) == "try-error")) {lp_only_validation <- manual_lp_only_validation$lp_new_data}
          if (! TRUE %in% (class(manual_lp_all_subjects) == "try-error")) {lp_all_subjects <- manual_lp_all_subjects$lp_new_data}
        } else {
          # Obtain the transformed linear predictors ####
          if (outcome_type != "time-to-event") {
            lp_training <- lp_training_with_se$fit
            lp_only_validation <- lp_only_validation_with_se$fit
            lp_all_subjects <- lp_all_subjects_with_se$fit
          } else {
            lp_training <- 1 - exp(-lp_training_with_se$fit)
            lp_only_validation <- 1 - exp(-lp_only_validation_with_se$fit)
            lp_all_subjects <- 1 - exp(-lp_all_subjects_with_se$fit)
          }
          # Impute missing linear predictors from manual predictions ####
          {
            if (! TRUE %in% (class(manual_lp_training) == "try-error")) {lp_training[is.na(lp_training)] <- manual_lp_training$lp_new_data[is.na(lp_training)]}
            if (! TRUE %in% (class(manual_lp_only_validation) == "try-error")) {lp_only_validation[is.na(lp_only_validation)] <- manual_lp_only_validation$lp_new_data[is.na(lp_only_validation)]}
            if (! TRUE %in% (class(manual_lp_all_subjects) == "try-error")) {lp_all_subjects[is.na(lp_all_subjects)] <- manual_lp_all_subjects$lp_new_data[is.na(lp_all_subjects)]}
          }
        }
        # Now perform intercept-slope adjustment ####
        if (outcome_type == 'time-to-event') {
          # Reverse the transformation 1 - (exp(-lp))
          # transformed_lp = 1 - (exp(-untransformed_lp))
          # exp(-untransformed_lp) = 1 - transformed_lp
          # -untransformed_lp = log(1 - transformed_lp)
          # untransformed_lp = -log(1-transformed_lp)
          # To avoid error due to lp being exactly 1
          lp_training[lp_training == 1] <- 1 - 0.000001
          # When the results are NaN, then assign a probability of event of 0
          lp_training[is.nan(lp_training)] <- 0
          back_transformed_lp_training <- (-log(1-lp_training))
        } else if (outcome_type == 'binary') {
          # Reverse the transformation plogis(lp)
          # transformed_lp = plogis(untransformed_lp)
          # untransformed_lp = qlogis(transformed_lp)
          # To avoid error due to lp being exactly 0 or 1
          lp_training[lp_training == 0] <- 0.000001
          lp_training[lp_training == 1] <- 1- 0.000001
          back_transformed_lp_training <- qlogis(lp_training)
        } else if (outcome_count == TRUE) {
          # Reverse the transformation exp(lp)
          # transformed_lp = exp(untransformed_lp)
          # untransformed_lp = log(transformed_lp)
          # To avoid error due to lp being exactly 1
          lp_training[lp_training == 1] <- 1 - 0.000001
          back_transformed_lp_training <- log(lp_training)
        } else {
          # For continuous outcomes, no transformation was performed; therefore, none is required
          back_transformed_lp_training <- lp_training
        }
        # If all the linear predictors are of the same value, regression cannot be achieved. Therefore, add a very small value (0.000001) to first sample
        if (length(unique(round(back_transformed_lp_training,6))) == 1) {
          back_transformed_lp_training[1] <- back_transformed_lp_training[1] + 0.000001
        }
        intercept_slope_adjustment_text <- create_intercept_slope_adjustment_text(outcome_name, outcome_type, outcome_time, outcome_count)
        intercept_slope_adjustment_model <- try(eval(parse(text = intercept_slope_adjustment_text)), silent = TRUE)
        if (TRUE %in% (class(intercept_slope_adjustment_model == "try-error"))) {
          predicted_training_calibration_adjusted = NA
          predicted_only_validation_calibration_adjusted = NA
          predicted_all_subjects_calibration_adjusted = NA
          lp_training_calibration_adjusted = NA
          lp_only_validation_calibration_adjusted = NA
          lp_all_subjects_calibration_adjusted = NA
        } else {
          if (! TRUE %in% (class(manual_lp_training) == "try-error")) {
            manual_lp_training_calibration_adjusted <- calculate_intercept_slope_adjusted_lp(
              outcome_type = outcome_type, outcome_count = outcome_count,
              intercept_slope_adjustment_model = intercept_slope_adjustment_model, manual_lp = manual_lp_training
            )
            lp_training_calibration_adjusted = manual_lp_all_subjects_calibration_adjusted$lp_new_data
          } else {
            lp_training_calibration_adjusted = NA
          }
          if (! TRUE %in% (class(manual_lp_only_validation) == "try-error")) {
            manual_lp_only_validation_calibration_adjusted <- calculate_intercept_slope_adjusted_lp(
              outcome_type = outcome_type, outcome_count = outcome_count,
              intercept_slope_adjustment_model = intercept_slope_adjustment_model, manual_lp = manual_lp_only_validation
            )
            lp_only_validation_calibration_adjusted = manual_lp_only_validation_calibration_adjusted$lp_new_data
          } else {
            lp_only_validation_calibration_adjusted = NA
          }
          if (! TRUE %in% (class(manual_lp_all_subjects) == "try-error")) {
            manual_lp_all_subjects_calibration_adjusted <- calculate_intercept_slope_adjusted_lp(
              outcome_type = outcome_type, outcome_count = outcome_count,
              intercept_slope_adjustment_model = intercept_slope_adjustment_model, manual_lp = manual_lp_all_subjects
            )
            lp_all_subjects_calibration_adjusted = manual_lp_all_subjects_calibration_adjusted$lp_new_data
          } else {
            lp_all_subjects_calibration_adjusted = NA
          }
        }
        # Adjusted lp (include coefficient of mandatory predictors only) ####
        if (length(mandatory_predictors) > 0) {
          lp_training_adjusted_mandatory_predictors_only <- try(calculate_manual_lp(outcome_type = outcome_type, outcome_count = outcome_count,
                                                                                new_data = df_training_complete_copy, regression_model = regression_model,
                                                                                variables_in_model_summary = variables_in_model_summary[variables_in_model_summary$name %in% mandatory_predictors,], df_training_complete = df_training_complete), silent = TRUE)
          if (! TRUE %in% (class(lp_training_adjusted_mandatory_predictors_only) == "try-error")) {
            lp_training_adjusted_mandatory_predictors_only <- lp_training_adjusted_mandatory_predictors_only$lp_new_data
          } else {
            lp_training_adjusted_mandatory_predictors_only <- NA
          }
          lp_only_validation_adjusted_mandatory_predictors_only <- try(calculate_manual_lp(outcome_type = outcome_type, outcome_count = outcome_count,
                                                                                       new_data = df_only_validation, regression_model = regression_model,
                                                                                       variables_in_model_summary = variables_in_model_summary[variables_in_model_summary$name %in% mandatory_predictors,], df_training_complete = df_training_complete), silent = TRUE)
          if (! TRUE %in% (class(lp_only_validation_adjusted_mandatory_predictors_only) == "try-error")) {
            lp_only_validation_adjusted_mandatory_predictors_only <- lp_only_validation_adjusted_mandatory_predictors_only$lp_new_data
          } else {
            lp_only_validation_adjusted_mandatory_predictors_only <- NA
          }
          lp_all_subjects_adjusted_mandatory_predictors_only <- try(calculate_manual_lp(outcome_type = outcome_type, outcome_count = outcome_count,
                                                                                    new_data = df_all_subjects, regression_model = regression_model,
                                                                                    variables_in_model_summary = variables_in_model_summary[variables_in_model_summary$name %in% mandatory_predictors,], df_training_complete = df_training_complete), silent = TRUE)
          if (! TRUE %in% (class(lp_all_subjects_adjusted_mandatory_predictors_only) == "try-error")) {
            lp_all_subjects_adjusted_mandatory_predictors_only <- lp_all_subjects_adjusted_mandatory_predictors_only$lp_new_data
          } else {
            lp_all_subjects_adjusted_mandatory_predictors_only <- NA
          }
        } else {
          lp_training_adjusted_mandatory_predictors_only <- NA
          lp_only_validation_adjusted_mandatory_predictors_only <- NA
          lp_all_subjects_adjusted_mandatory_predictors_only <- NA
        }
        # Predictions ####
        # If the outcome type is quantitative, the lp gives the predictions
        if (outcome_type == "quantitative") {
          predicted_training <- lp_training
          predicted_only_validation <- lp_only_validation
          predicted_all_subjects <- lp_all_subjects
          predicted_training_calibration_adjusted <- lp_training_calibration_adjusted
          predicted_only_validation_calibration_adjusted <- lp_only_validation_calibration_adjusted
          predicted_all_subjects_calibration_adjusted <- lp_all_subjects_calibration_adjusted
          predicted_training_adjusted_mandatory_predictors_only <- lp_training_adjusted_mandatory_predictors_only
          predicted_only_validation_adjusted_mandatory_predictors_only <- lp_only_validation_adjusted_mandatory_predictors_only
          predicted_all_subjects_adjusted_mandatory_predictors_only <- lp_all_subjects_adjusted_mandatory_predictors_only
        } else {
          predicted_training <- NA
          predicted_only_validation <- NA
          predicted_all_subjects <- NA
          predicted_training_calibration_adjusted <- NA
          predicted_only_validation_calibration_adjusted <- NA
          predicted_all_subjects_calibration_adjusted <- NA
          predicted_training_adjusted_mandatory_predictors_only <- NA
          predicted_only_validation_adjusted_mandatory_predictors_only <- NA
          predicted_all_subjects_adjusted_mandatory_predictors_only <- NA
          # First standard lp ####
          {
            optimal_threshold <- calculate_pROC_threshold(lp_training, actual_training)
            # If the pROC threshold fails or if the heuristic threshold is the preferred method, try the heuristic threshold based on prevalence
            if (is.na(optimal_threshold$direction) |
                ((is.na(optimal_threshold$threshold_youden)) & (is.na(optimal_threshold$threshold_topleft))) |
                (model_threshold_method == "heuristic")) {
              optimal_threshold <- calculate_heuristic_threshold(lp_training, actual_training)
            }
            # If there is no direction for the threshold calculations, no classification is possible
            if (is.na(optimal_threshold$direction)) {
              predicted_training = NA
              predicted_only_validation = NA
              predicted_all_subjects = NA
            } else {
              # If the youden or topleft methods could not be calculated or if the heuristic method is the preferred method, then the threshold is the heuristic threshold
              threshold <- if (is.null(optimal_threshold$threshold_youden)) {
                optimal_threshold$threshold
              } else {
                # If youden is the preferred method, use the youden threshold, but if this is absent, use the topleft threshold
                # On the other hand, if topleft is the preferred method, use the topleft threshold, but if this is absent, use the youden threshold
                if (model_threshold_method == "youden") {
                  ifelse(is.na(optimal_threshold$threshold_youden), optimal_threshold$threshold_topleft, optimal_threshold$threshold_youden)
                } else {
                  ifelse(is.na(optimal_threshold$threshold_topleft), optimal_threshold$threshold_youden, optimal_threshold$threshold_topleft)
                }
              }
              if (optimal_threshold$direction == "<") {
                predicted_training[lp_training <= threshold] <- levels(actual_training)[1]
                predicted_training[lp_training > threshold] <- levels(actual_training)[2]
                predicted_only_validation[lp_only_validation <= threshold] <- levels(actual_all_subjects)[1]
                predicted_only_validation[lp_only_validation > threshold] <- levels(actual_all_subjects)[2]
                predicted_all_subjects[lp_all_subjects <= threshold] <- levels(actual_all_subjects)[1]
                predicted_all_subjects[lp_all_subjects > threshold] <- levels(actual_all_subjects)[2]
              } else {
                predicted_training[lp_training <= threshold] <- levels(actual_training)[2]
                predicted_training[lp_training > threshold] <- levels(actual_training)[1]
                predicted_only_validation[lp_only_validation <= threshold] <- levels(actual_all_subjects)[2]
                predicted_only_validation[lp_only_validation > threshold] <- levels(actual_all_subjects)[1]
                predicted_all_subjects[lp_all_subjects <= threshold] <- levels(actual_all_subjects)[2]
                predicted_all_subjects[lp_all_subjects > threshold] <- levels(actual_all_subjects)[1]
              }
              predicted_training <- factor(predicted_training, levels = levels(actual_training))
              predicted_only_validation <- factor(predicted_only_validation, levels = levels(actual_only_validation))
              predicted_all_subjects <- factor(predicted_all_subjects, levels = levels(actual_all_subjects))
            }
          }
          # Next calibration-adjusted lp ####
          {
            if (length(lp_training_calibration_adjusted) > 1) {
              optimal_threshold_calibration_adjusted <- calculate_pROC_threshold(lp_training_calibration_adjusted, actual_training)
            # If the pROC threshold fails or if the heuristic threshold is the preferred method, try the heuristic threshold based on prevalence
            if (is.na(optimal_threshold_calibration_adjusted$direction) |
                ((is.na(optimal_threshold_calibration_adjusted$threshold_youden)) & (is.na(optimal_threshold_calibration_adjusted$threshold_topleft))) |
                (model_threshold_method == "heuristic")) {
              optimal_threshold_calibration_adjusted <- calculate_heuristic_threshold(lp_training_calibration_adjusted, actual_training)
            }
            # If there is no direction for the threshold calculations, no classification is possible
            if (is.na(optimal_threshold_calibration_adjusted$direction)) {
              predicted_training_calibration_adjusted = NA
              predicted_only_validation_calibration_adjusted = NA
              predicted_all_subjects_calibration_adjusted = NA
            } else {
              # If the youden or topleft methods could not be calculated or if the heuristic method is the preferred method, then the threshold is the heuristic threshold
              threshold_calibration_adjusted <- if (is.null(optimal_threshold_calibration_adjusted$threshold_youden)) {
                optimal_threshold_calibration_adjusted$threshold
              } else {
                # If youden is the preferred method, use the youden threshold, but if this is absent, use the topleft threshold
                # On the other hand, if topleft is the preferred method, use the topleft threshold, but if this is absent, use the youden threshold
                if (model_threshold_method == "youden") {
                  ifelse(is.na(optimal_threshold_calibration_adjusted$threshold_youden), optimal_threshold_calibration_adjusted$threshold_topleft, optimal_threshold_calibration_adjusted$threshold_youden)
                } else {
                  ifelse(is.na(optimal_threshold_calibration_adjusted$threshold_topleft), optimal_threshold_calibration_adjusted$threshold_youden, optimal_threshold_calibration_adjusted$threshold_topleft)
                }
              }
              if (optimal_threshold_calibration_adjusted$direction == "<") {
                predicted_training_calibration_adjusted[lp_training_calibration_adjusted <= threshold_calibration_adjusted] <- levels(actual_training)[1]
                predicted_training_calibration_adjusted[lp_training_calibration_adjusted > threshold_calibration_adjusted] <- levels(actual_training)[2]
                predicted_only_validation_calibration_adjusted[lp_only_validation_calibration_adjusted <= threshold_calibration_adjusted] <- levels(actual_only_validation)[1]
                predicted_only_validation_calibration_adjusted[lp_only_validation_calibration_adjusted > threshold_calibration_adjusted] <- levels(actual_only_validation)[2]
                predicted_all_subjects_calibration_adjusted[lp_all_subjects_calibration_adjusted <= threshold_calibration_adjusted] <- levels(actual_all_subjects)[1]
                predicted_all_subjects_calibration_adjusted[lp_all_subjects_calibration_adjusted > threshold_calibration_adjusted] <- levels(actual_all_subjects)[2]
              } else {
                predicted_training_calibration_adjusted[lp_training_calibration_adjusted <= threshold_calibration_adjusted] <- levels(actual_training)[2]
                predicted_training_calibration_adjusted[lp_training_calibration_adjusted > threshold_calibration_adjusted] <- levels(actual_training)[1]
                predicted_only_validation_calibration_adjusted[lp_only_validation_calibration_adjusted <= threshold_calibration_adjusted] <- levels(actual_only_validation)[2]
                predicted_only_validation_calibration_adjusted[lp_only_validation_calibration_adjusted > threshold_calibration_adjusted] <- levels(actual_only_validation)[1]
                predicted_all_subjects_calibration_adjusted[lp_all_subjects_calibration_adjusted <= threshold_calibration_adjusted] <- levels(actual_all_subjects)[2]
                predicted_all_subjects_calibration_adjusted[lp_all_subjects_calibration_adjusted > threshold_calibration_adjusted] <- levels(actual_all_subjects)[1]
              }
              predicted_training_calibration_adjusted <- factor(predicted_training_calibration_adjusted, levels = levels(actual_training))
              predicted_only_validation_calibration_adjusted <- factor(predicted_only_validation_calibration_adjusted, levels = levels(actual_only_validation))
              predicted_all_subjects_calibration_adjusted <- factor(predicted_all_subjects_calibration_adjusted, levels = levels(actual_all_subjects))
            }
            }
          }
          # Next adjusted mandatory predictors only lp ####
          {
            if (length(lp_training_adjusted_mandatory_predictors_only) > 1) {
            optimal_threshold_adjusted_mandatory_predictors_only <- calculate_pROC_threshold(lp_training_adjusted_mandatory_predictors_only, actual_training)
            # If the pROC threshold fails or if the heuristic threshold is the preferred method, try the heuristic threshold based on prevalence
            if (is.na(optimal_threshold_adjusted_mandatory_predictors_only$direction) |
                ((is.na(optimal_threshold_adjusted_mandatory_predictors_only$threshold_youden)) & (is.na(optimal_threshold_adjusted_mandatory_predictors_only$threshold_topleft))) |
                (model_threshold_method == "heuristic")) {
              optimal_threshold_adjusted_mandatory_predictors_only <- calculate_heuristic_threshold(lp_training_adjusted_mandatory_predictors_only, actual_training)
            }
            # If there is no direction for the threshold calculations, no classification is possible
            if (is.na(optimal_threshold_adjusted_mandatory_predictors_only$direction)) {
              predicted_training_adjusted_mandatory_predictors_only = NA
              predicted_only_validation_adjusted_mandatory_predictors_only = NA
              predicted_all_subjects_adjusted_mandatory_predictors_only = NA
            } else {
              # If the youden or topleft methods could not be calculated or if the heuristic method is the preferred method, then the threshold is the heuristic threshold
              threshold_adjusted_mandatory_predictors_only <- if (is.null(optimal_threshold_adjusted_mandatory_predictors_only$threshold_youden)) {
                optimal_threshold_adjusted_mandatory_predictors_only$threshold
              } else {
                # If youden is the preferred method, use the youden threshold, but if this is absent, use the topleft threshold
                # On the other hand, if topleft is the preferred method, use the topleft threshold, but if this is absent, use the youden threshold
                if (model_threshold_method == "youden") {
                  ifelse(is.na(optimal_threshold_adjusted_mandatory_predictors_only$threshold_youden), optimal_threshold_adjusted_mandatory_predictors_only$threshold_topleft, optimal_threshold_adjusted_mandatory_predictors_only$threshold_youden)
                } else {
                  ifelse(is.na(optimal_threshold_adjusted_mandatory_predictors_only$threshold_topleft), optimal_threshold_adjusted_mandatory_predictors_only$threshold_youden, optimal_threshold_adjusted_mandatory_predictors_only$threshold_topleft)
                }
              }
              if (optimal_threshold_adjusted_mandatory_predictors_only$direction == "<") {
                predicted_training_adjusted_mandatory_predictors_only[lp_training_adjusted_mandatory_predictors_only <= threshold_adjusted_mandatory_predictors_only] <- levels(actual_training)[1]
                predicted_training_adjusted_mandatory_predictors_only[lp_training_adjusted_mandatory_predictors_only > threshold_adjusted_mandatory_predictors_only] <- levels(actual_training)[2]
                predicted_only_validation_adjusted_mandatory_predictors_only[lp_only_validation_adjusted_mandatory_predictors_only <= threshold_adjusted_mandatory_predictors_only] <- levels(actual_only_validation)[1]
                predicted_only_validation_adjusted_mandatory_predictors_only[lp_only_validation_adjusted_mandatory_predictors_only > threshold_adjusted_mandatory_predictors_only] <- levels(actual_only_validation)[2]
                predicted_all_subjects_adjusted_mandatory_predictors_only[lp_all_subjects_adjusted_mandatory_predictors_only <= threshold_adjusted_mandatory_predictors_only] <- levels(actual_all_subjects)[1]
                predicted_all_subjects_adjusted_mandatory_predictors_only[lp_all_subjects_adjusted_mandatory_predictors_only > threshold_adjusted_mandatory_predictors_only] <- levels(actual_all_subjects)[2]
              } else {
                predicted_training_adjusted_mandatory_predictors_only[lp_training_adjusted_mandatory_predictors_only <= threshold_adjusted_mandatory_predictors_only] <- levels(actual_training)[2]
                predicted_training_adjusted_mandatory_predictors_only[lp_training_adjusted_mandatory_predictors_only > threshold_adjusted_mandatory_predictors_only] <- levels(actual_training)[1]
                predicted_only_validation_adjusted_mandatory_predictors_only[lp_only_validation_adjusted_mandatory_predictors_only <= threshold_adjusted_mandatory_predictors_only] <- levels(actual_only_validation)[2]
                predicted_only_validation_adjusted_mandatory_predictors_only[lp_only_validation_adjusted_mandatory_predictors_only > threshold_adjusted_mandatory_predictors_only] <- levels(actual_only_validation)[1]
                predicted_all_subjects_adjusted_mandatory_predictors_only[lp_all_subjects_adjusted_mandatory_predictors_only <= threshold_adjusted_mandatory_predictors_only] <- levels(actual_all_subjects)[2]
                predicted_all_subjects_adjusted_mandatory_predictors_only[lp_all_subjects_adjusted_mandatory_predictors_only > threshold_adjusted_mandatory_predictors_only] <- levels(actual_all_subjects)[1]
              }
              predicted_training_adjusted_mandatory_predictors_only <- factor(predicted_training_adjusted_mandatory_predictors_only, levels = levels(actual_training))
              predicted_only_validation_adjusted_mandatory_predictors_only <- factor(predicted_only_validation_adjusted_mandatory_predictors_only, levels = levels(actual_only_validation))
              predicted_all_subjects_adjusted_mandatory_predictors_only <- factor(predicted_all_subjects_adjusted_mandatory_predictors_only, levels = levels(actual_all_subjects))
            }
            }
          }
        }
        output <- {list(
          actual_training = actual_training,
          predicted_training = predicted_training,
          predicted_training_calibration_adjusted = predicted_training_calibration_adjusted,
          predicted_training_adjusted_mandatory_predictors_only = predicted_training_adjusted_mandatory_predictors_only,
          actual_only_validation = actual_only_validation,
          predicted_only_validation = predicted_only_validation,
          predicted_only_validation_calibration_adjusted = predicted_only_validation_calibration_adjusted,
          predicted_only_validation_adjusted_mandatory_predictors_only = predicted_only_validation_adjusted_mandatory_predictors_only,
          actual_all_subjects = actual_all_subjects,
          predicted_all_subjects = predicted_all_subjects,
          predicted_all_subjects_calibration_adjusted = predicted_all_subjects_calibration_adjusted,
          predicted_all_subjects_adjusted_mandatory_predictors_only = predicted_all_subjects_adjusted_mandatory_predictors_only,
          lp_training = lp_training,
          lp_only_validation = lp_only_validation,
          lp_all_subjects = lp_all_subjects,
          lp_training_calibration_adjusted = lp_training_calibration_adjusted,
          lp_only_validation_calibration_adjusted = lp_only_validation_calibration_adjusted,
          lp_all_subjects_calibration_adjusted = lp_all_subjects_calibration_adjusted,
          lp_training_adjusted_mandatory_predictors_only = lp_training_adjusted_mandatory_predictors_only,
          lp_only_validation_adjusted_mandatory_predictors_only = lp_only_validation_adjusted_mandatory_predictors_only,
          lp_all_subjects_adjusted_mandatory_predictors_only = lp_all_subjects_adjusted_mandatory_predictors_only,
          time_training = if (outcome_type == "time-to-event") {df_training_complete_copy[,outcome_time]} else {NA},
          time_only_validation = if (outcome_type == "time-to-event") {df_only_validation[,outcome_time]} else {NA},
          time_all_subjects = if (outcome_type == "time-to-event") {df_all_subjects[,outcome_time]} else {NA},
          regression_model = regression_model,
          df_all_subjects = df_all_subjects,
          html_file = html_file,
          outcome = "Successful"
        )}
      }
    }
  }
  return(output)
}
