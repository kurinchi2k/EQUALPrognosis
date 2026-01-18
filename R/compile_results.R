compile_results <- function(generic_input_parameters, specific_input_parameters, prepared_datasets) {
  # Functions ####
  figure_to_html_code <- function(caption, figure_path, alt = "see text", width_percent = 100) {
    mime_type <- mime::guess_type(figure_path)
    if (width_percent > 100) {width_percent <- 100}
    if (width_percent < 20) {width_percent <- 20}
    paste0('<figure>','\n',
           '<figcaption>',caption,'</figcaption>','\n',
           '<img src="', base64enc::dataURI(file = figure_path, mime = mime_type), '" alt="',alt,'" style="width:',width_percent,'%">','\n',
           '</figure>','\n'
    )
  }
  table_to_html_code <- function(caption, df, width_percent = 100) {
    if (width_percent > 100) {width_percent <- 100}
    if (width_percent < 20) {width_percent <- 20}
    paste0('<table style="width:',width_percent,'%">','\n',
           '<caption>','\n',
           caption,'\n',
           '</caption>','\n',
           paste0('<tr>', paste0('<th>', colnames(df), '</th>', collapse = '\n'), '</tr>'), '\n',
           paste0(unlist(lapply(1:nrow(df), function(x) {paste0('<tr>', paste0('<td>', df[x,], '</td>', collapse = '\n'), '</tr>')})), collapse = "\n"),'\n',
           '</table>','\n'
    )
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
  baseline_html_code <- function(figure_border = TRUE, table_double_border = FALSE, table_header_background_colour = "darkblue", table_header_font_colour = "white", table_text_center_alignment = FALSE, caption_background_colour = "black", caption_font_colour = "white", caption_center = FALSE, caption_pixel_size = 16, caption_italics = FALSE, caption_bold = TRUE) {
    paste0(
      '<html>','\n',
      '<head>','\n',
      '<style>','\n',
      'figure {','\n',
      if (figure_border == TRUE) {
        paste0(
          'border: 1px #cccccc solid;','\n',
          'padding: 4px;','\n'
        )
      },
      'margin: auto;','\n',
      '}','\n',
      'figcaption {','\n',
      'background-color: ', caption_background_colour, ';','\n',
      'color: ', caption_font_colour, ';','\n',
      'font-size: ', caption_pixel_size, 'px;','\n',
      if (caption_italics == TRUE) {paste0('font-style: italic;', '\n')},
      if (caption_bold == TRUE) {paste0('font-weight: bold;', '\n')},
      'padding: 2px;', '\n',
      if (caption_center == TRUE) {'text-align: center;'} else {'text-align: left;'}, '\n',
      '}','\n',
      'table, th, td {','\n',
      'border: 1px solid black;','\n',
      if (table_double_border == FALSE) {paste0('border-collapse: collapse;','\n')},
      if (table_text_center_alignment == FALSE) {paste0('text-align: left;','\n')} else {paste0('text-align: center;','\n')},
      '}','\n',
      'th {','\n',
      'background-color: ', table_header_background_colour, ';','\n',
      'color: ', table_header_font_colour, ';','\n',
      '}','\n',
      'caption {','\n',
      'background-color: ', caption_background_colour, ';','\n',
      'color: ', caption_font_colour, ';','\n',
      'font-size: ', caption_pixel_size, 'px;','\n',
      if (caption_italics == TRUE) {paste0('font-style: italic;', '\n')},
      if (caption_bold == TRUE) {paste0('font-weight: bold;', '\n')},
      'padding: 2px;', '\n',
      if (caption_center == TRUE) {'text-align: center;'} else {'text-align: left;'}, '\n',
      '}','\n',
      '</style>','\n',
      '</head>','\n',
      '<body>','\n'
    )
  }
  estimate_with_ci_forest <- function(est_ci_table, parameter_name = NA, minimum = NA, maximum = NA, reference_line = NA,
                                      convert_to_html_code = TRUE) {
    est_ci_table$analysis <- factor(est_ci_table$analysis, levels = rev(est_ci_table$analysis))
    plot_text <- paste0(
      c(
        "ggplot(data = est_ci_table, aes(x=analysis, y=est, ymin = lci, ymax = uci, colour = analysis)) +
      geom_pointrange() +
      coord_flip() +
      theme(legend.position = 'none') +
      xlab('Model')",
        if (! is.na(reference_line)) {paste0("\ngeom_hline(yintercept = ", reference_line, ", linetype = 'dotted')")},
        if ((! is.na(minimum)) | (! is.na(maximum))) {paste0("\nscale_y_continuous(limits = c(", minimum, ",", maximum,"))")},
        if (! is.na(parameter_name)) {paste0("\nylab('", parameter_name, "')")}
      ),
      collapse = " + "
    )
    placeholder <- eval(parse(text = plot_text))
    if (convert_to_html_code == TRUE) {
      temp_file <- tempfile(fileext = ".png")
      suppressMessages(suppressWarnings(ggsave(temp_file, width = 7, height = nrow(est_ci_table))))
      output <- figure_to_html_code(caption = ifelse(is.na(parameter_name), "Parameter", parameter_name), figure_path = temp_file, alt = ifelse(is.na(parameter_name), "Parameter", parameter_name), width_percent = 100)
      unlink(temp_file)
    } else {
      output <- placeholder
    }
    return(output)
  }
  backtransform_lp <- function(lp, outcome_type, outcome_count) {
    if (outcome_type == 'time-to-event') {
      # Reverse the transformation 1 - (exp(-lp))
      # transformed_lp = 1 - (exp(-untransformed_lp))
      # exp(-untransformed_lp) = 1 - transformed_lp
      # -untransformed_lp = log(1 - transformed_lp)
      # untransformed_lp = -log(1-transformed_lp)
      # To avoid error due to lp being exactly 1
      lp[lp == 1] <- 1 - 0.000001
      # When the results are NaN, then assign a probability of event of 0
      lp[is.nan(lp)] <- 0
      back_transformed_lp <- (-log(1-lp))
    } else if (outcome_type == 'binary') {
      # Reverse the transformation plogis(lp)
      # transformed_lp = plogis(untransformed_lp)
      # untransformed_lp = qlogis(transformed_lp)
      # To avoid error due to lp being exactly 0 or 1
      lp[lp == 0] <- 0.000001
      lp[lp == 1] <- 1- 0.000001
      back_transformed_lp <- qlogis(lp)
    } else if (outcome_count == TRUE) {
      # Reverse the transformation exp(lp)
      # transformed_lp = exp(untransformed_lp)
      # untransformed_lp = log(transformed_lp)
      # To avoid error due to lp being exactly 1
      lp[lp == 1] <- 1 - 0.000001
      back_transformed_lp <- log(lp)
    } else {
      # For continuous outcomes, no transformation was performed; therefore, none is required
      back_transformed_lp <- lp
    }
    return(back_transformed_lp)
  }
  # Dummy variables to avoid CRAN errors. These variables are related to the names of variables in box plot.
  group <- NA
  optimism <- NA
  # Baseline html file ####
  html_file <- {list(baseline_code = baseline_html_code(),
                     basic_details = {paste0(
                       '<h1>',generic_input_parameters$general_title,'</h1>',
                       '<h2>Database structure</h2>',
                       r_object_to_html_code(str(generic_input_parameters$df)),
                       '<h2>Basic details (common across all analyses)</h2>',
                       '<p>',
                       'Number of simulations: ', '<i>', generic_input_parameters$simulations, '</i>', '<br>',
                       'Outcome name: ', '<i>', generic_input_parameters$outcome_name, '</i>', '<br>',
                       'Outcome type: ', '<i>', generic_input_parameters$outcome_type, '</i>', '<br>',
                       'Outcome time: ', '<i>', ifelse(is.na(generic_input_parameters$outcome_time), "Not applicable", generic_input_parameters$outcome_time), '</i>', '<br>',
                       'Is this a count outcome: ', '<i>', generic_input_parameters$outcome_count, '</i>', '<br>'
                     )}
  )}
  # html code for specific_input_parameters ####
  html_specific_input_parameters <- lapply(specific_input_parameters, function(x) {
    paste0(
      '<h2>', 'Analysis: ', x$analysis_title, '</h2>',
      'Develop model: ', '<i>', x$develop_model, '</i>', '<br>',
      'Mandatory predictors: ', '<i>',
      ifelse(is.null(x$mandatory_predictors), "None",
             ifelse(is.na(x$mandatory_predictors), "Not applicable",
                    paste0(x$mandatory_predictors, collapse = ", "))),
      '</i>', '<br>',
      'Optional predictors: ', '<i>',
      ifelse(is.null(x$optional_predictors), "None",
             ifelse(is.na(x$optional_predictors), "Not applicable",
                    paste0(x$optional_predictors, collapse = ", "))),
      '</i>', '<br>',
      'Mandatory interactions: ', '<i>',
      ifelse(is.null(x$mandatory_interactions), "None",
             ifelse(is.na(x$mandatory_interactions), "Not applicable",
                    paste0(x$mandatory_interactions, collapse = ", "))),
      '</i>', '<br>',
      'Optional interactions: ', '<i>',
      ifelse(is.null(x$optional_interactions), "None",
             ifelse(is.na(x$optional_interactions), "Not applicable",
                    paste0(x$optional_interactions, collapse = ", "))),
      '</i>', '<br>',
      'Method used for threshold development: ', '<i>',
      ifelse(is.na(x$model_threshold_method), "Not applicable",
             x$model_threshold_method), '</i>', '<br>',
      'Scoring system: ', '<i>',
      ifelse(is.na(x$scoring_system), "Not applicable", x$scoring_system),
      '</i>', '<br>',
      'Predetermined threshold: ', '<i>',
      ifelse(is.na(x$predetermined_threshold), "Not applicable",
             x$predetermined_threshold), '</i>', '<br>',
      'Do higher values indicate event: ', '<i>',
      ifelse(is.na(x$higher_values_event), "Not applicable",
             x$higher_values_event), '</i>',
      '</p>'
    )
  })
  names(html_specific_input_parameters) <- names(specific_input_parameters)
  # Each analysis ####
  results <- lapply(1:length(specific_input_parameters), function(x) {
    cat(paste0("\n", names(specific_input_parameters)[x], "...", "\n"))
    output <- try(perform_analysis(
      generic_input_parameters = generic_input_parameters,
      specific_input_parameters_each_analysis = specific_input_parameters[[x]],
      prepared_datasets = prepared_datasets), silent = TRUE)
  })
  names(results) <- names(specific_input_parameters)
  # Develop html code for each analysis and calibration plots ####
  cat("\nSummarising results...")
  html_each_analysis <- unlist(lapply(1:length(specific_input_parameters), function(x) {
    if (TRUE %in% (class(results[[x]]) == "try-error")) {
      apparent_performance <- c(
        '<h2>Apparent performance</h2>',
        'Apparent performance cannot be calculated because of an error in running the model. Please try a simpler model.'
      )
      test_performance <- c(
        '<h2>Performance in simulations</h2>',
        'Performance in simulations cannot be calculated because of an error in running the model. Please try a simpler model.'
      )
    } else {
      # Apparent performance ####
      if (generic_input_parameters$outcome_type != "time-to-event") {
        actual <- results[[x]]$actual_predicted_results_apparent$actual_all_subjects
        lp <- results[[x]]$actual_predicted_results_apparent$lp_all_subjects
        # remove missing lp
        actual <- actual[! is.na(lp)]
        lp <- lp[! is.na(lp)]
        if (length(actual) > 0) {
          if (generic_input_parameters$outcome_type != "quantitative") {
            actual_values <- rep(NA, length(actual))
            actual_values[actual == levels(actual)[1]] <- 0
            actual_values[actual == levels(actual)[2]] <- 1
            calib_plot <- suppressWarnings(try(valProbggplot(lp, actual_values, allowPerfectPredictions = TRUE), silent = TRUE))
            if (TRUE %in% (class(calib_plot) == "try-error")) {
              # If this resulted in error use predtools. However this is based on untransformed lp. The backtransformation is applied at the stage of creating the calibration_plot_data
              calibration_plot_data <- cbind.data.frame(
                actual = actual_values,
                predicted = backtransform_lp(lp = lp, outcome_type = generic_input_parameters$outcome_type, outcome_count = generic_input_parameters$outcome_count)
                )
              calib_plot <- suppressWarnings(try(calibration_plot(data = calibration_plot_data, obs = "actual", pred = "predicted", y_lim = c(0, 1), x_lim=c(0, 1), title = "Calibration plot"), silent = TRUE))
              # If even this resulted in error, indicate this.
              if (TRUE %in% (class(calib_plot) == "try-error")) {
                apparent_performance <- c(
                  '<h2>Apparent performance</h2>',
                  'Apparent performance cannot be calculated because of an error in estimating the calibration plot.'
                )
              } else {
                temp_file <- tempfile(fileext = ".png")
                suppressMessages(suppressWarnings(ggsave(temp_file)))
                apparent_performance <- c(
                  '<h2>Apparent performance</h2>',
                  figure_to_html_code(caption = "Apparent performance", figure_path = temp_file, alt = paste0("Calibration plot"), width_percent = 100)
                )
                unlink(temp_file)
              }
            } else {
              # ggplot version available for binary and time-to-event outcomes
              temp_file <- tempfile(fileext = ".png")
              suppressMessages(suppressWarnings(ggsave(temp_file)))
              apparent_performance <- c(
                '<h2>Apparent performance</h2>',
                figure_to_html_code(caption = "Apparent performance", figure_path = temp_file, alt = paste0("Calibration plot"), width_percent = 100)
              )
              unlink(temp_file)
            }
          } else {
            actual_values <- actual
            calib_plot <- suppressWarnings(try(genCalCurve(actual_values, lp, family = ifelse(generic_input_parameters$outcome_count, "poisson", "gaussian")), silent = TRUE))
            if (TRUE %in% (class(calib_plot) == "try-error")) {
              # If this resulted in error use predtools. However this is based on untransformed lp. The backtransformation is applied at the stage of creating the calibration_plot_data
              calibration_plot_data <- cbind.data.frame(
                actual = actual_values,
                predicted = backtransform_lp(lp = lp, outcome_type = generic_input_parameters$outcome_type, outcome_count = generic_input_parameters$outcome_count)
                )
              calib_plot <- suppressWarnings(try(calibration_plot(data = calibration_plot_data, obs = "actual", pred = "predicted", title = "Calibration plot"), silent = TRUE))
              # If even this resulted in error, indicate this.
              if (TRUE %in% (class(calib_plot) == "try-error")) {
                apparent_performance <- c(
                  '<h2>Apparent performance</h2>',
                  'Apparent performance cannot be calculated because of an error in estimating the calibration plot.'
                )
              } else {
                temp_file <- tempfile(fileext = ".png")
                suppressMessages(suppressWarnings(ggsave(temp_file)))
                apparent_performance <- c(
                  '<h2>Apparent performance</h2>',
                  figure_to_html_code(caption = "Apparent performance", figure_path = temp_file, alt = paste0("Calibration plot"), width_percent = 100)
                )
                unlink(temp_file)
              }
            } else {
              # ggplot version not available for GLM
              temp_file <- tempfile(fileext = ".png")
              png(temp_file)
              genCalCurve(actual_values, lp, family = ifelse(generic_input_parameters$outcome_count, "poisson", "gaussian"))
              recordPlot()
              dev.off()
              apparent_performance <- c(
                '<h2>Apparent performance</h2>',
                figure_to_html_code(caption = "Apparent performance", figure_path = temp_file, alt = paste0("Calibration plot"), width_percent = 100)
              )
              unlink(temp_file)
            }
          }
        } else {
          apparent_performance <- c(
            '<h2>Apparent performance</h2>',
            'Apparent performance cannot be calculated because linear predictor calculation is not applicable.'
          )
        }
      } else {
        calib_plot <- suppressWarnings(try(valProbSurvival(fit = results[[x]]$actual_predicted_results_apparent$regression_model,
                                                           valdata = results[[x]]$actual_predicted_results_apparent$df_all_subjects, plotCal = "ggplot"), silent = TRUE))
        if (TRUE %in% (class(calib_plot) == "try-error")) {
          # If this resulted in error use predtools. However this is based on untransformed lp. The backtransformation is applied at the stage of creating the calibration_plot_data
          actual <- results[[x]]$actual_predicted_results_apparent$actual_all_subjects
          lp <- results[[x]]$actual_predicted_results_apparent$lp_all_subjects
          follow_up <- results[[x]]$actual_predicted_results_apparent$time_all_subjects
          # remove missing lp
          actual <- actual[! is.na(lp)]
          lp <- lp[! is.na(lp)]
          follow_up <- follow_up[! is.na(lp)]
          if (length(actual) > 0) {
            actual_values <- rep(NA, length(actual))
            actual_values[actual == levels(actual)[1]] <- 0
            actual_values[actual == levels(actual)[2]] <- 1
            calibration_plot_data <- cbind.data.frame(
              actual = actual_values,
              predicted = backtransform_lp(lp = lp, outcome_type = generic_input_parameters$outcome_type, outcome_count = generic_input_parameters$outcome_count),
              follow_up = follow_up)
            calib_plot <- suppressWarnings(try(calibration_plot(data = calibration_plot_data, obs = "actual", follow_up = "follow_up", pred = "predicted", x_lim=c(0, 1), title = "Calibration plot"), silent = TRUE))
            # If even this resulted in error, ignore follow-up time
            if (TRUE %in% (class(calib_plot) == "try-error")) {
              calib_plot <- suppressWarnings(try(calibration_plot(data = calibration_plot_data, obs = "actual", pred = "predicted", y_lim = c(0, 1), x_lim=c(0, 1), title = "Calibration plot"), silent = TRUE))
              # If even this caused an error indicate so
              if (TRUE %in% (class(calib_plot) == "try-error")) {
              apparent_performance <- c(
                '<h2>Apparent performance</h2>',
                'Apparent performance cannot be calculated because of an error in estimating the calibration plot.'
              )
              } else {
                # Indicate that the calibration plot based on rate could not be calculated and with a warning that the plot was based ignoring the follow-up time
                temp_file <- tempfile(fileext = ".png")
                suppressMessages(suppressWarnings(ggsave(temp_file)))
                apparent_performance <-
                  c(
                    '<h2>Apparent performance</h2>',
                    'Apparent performance could not be calculated based on rate. Therefore, follow-up time was removed. Therefore, the calibration plot is likely to be unreliable.<br>',
                    figure_to_html_code(caption = "Apparent performance", figure_path = temp_file, alt = paste0("Calibration plot"), width_percent = 100)
                  )
                unlink(temp_file)
              }
            } else {
              temp_file <- tempfile(fileext = ".png")
              suppressMessages(suppressWarnings(ggsave(temp_file)))
              apparent_performance <- c(
                '<h2>Apparent performance</h2>',
                figure_to_html_code(caption = "Apparent performance", figure_path = temp_file, alt = paste0("Calibration plot"), width_percent = 100)
              )
              unlink(temp_file)
            }
          } else {
            apparent_performance <- c(
              '<h2>Apparent performance</h2>',
              'Apparent performance cannot be calculated because linear predictor calculation is not applicable.'
            )
          }
        } else {
          temp_file <- tempfile(fileext = ".png")
          suppressMessages(suppressWarnings(ggsave(temp_file)))
          apparent_performance <- c(
            '<h2>Apparent performance</h2>',
            figure_to_html_code(caption = "Apparent performance", figure_path = temp_file, alt = paste0("Calibration plot"), width_percent = 100)
          )
          unlink(temp_file)
        }
      }
      # Average of simulations ####
      if (generic_input_parameters$outcome_type != "time-to-event") {
        actual <- results[[x]]$actual_predicted_results_apparent$actual_all_subjects
        lp <- results[[x]]$average_lp_all_subjects
        # remove missing lp
        actual <- actual[! is.na(lp)]
        lp <- lp[! is.na(lp)]
        if (length(actual) > 0) {
          if (generic_input_parameters$outcome_type != "quantitative") {
            actual_values <- rep(NA, length(actual))
            actual_values[actual == levels(actual)[1]] <- 0
            actual_values[actual == levels(actual)[2]] <- 1
            calib_plot <- suppressWarnings(try(valProbggplot(lp, actual_values, allowPerfectPredictions = TRUE), silent = TRUE))
            if (TRUE %in% (class(calib_plot) == "try-error")) {
              # If this resulted in error use predtools. However this is based on untransformed lp. The backtransformation is applied at the stage of creating the calibration_plot_data
              calibration_plot_data <- cbind.data.frame(
                actual = actual_values,
                predicted = backtransform_lp(lp = lp, outcome_type = generic_input_parameters$outcome_type, outcome_count = generic_input_parameters$outcome_count)
                )
              calib_plot <- suppressWarnings(try(calibration_plot(data = calibration_plot_data, obs = "actual", pred = "predicted", y_lim = c(0, 1), x_lim=c(0, 1), title = "Calibration plot"), silent = TRUE))
              # If even this resulted in error, indicate this.
              if (TRUE %in% (class(calib_plot) == "try-error")) {
                test_performance <- c(
                  '<h2>Performance in simulations</h2>',
                  'Performance in simulations cannot be calculated because of an error in estimating the calibration curve.'
                )
              } else {
                temp_file <- tempfile(fileext = ".png")
                suppressMessages(suppressWarnings(ggsave(temp_file)))
                test_performance <- c(
                  '<h2>Performance in simulations</h2>',
                  figure_to_html_code(caption = "Performance in simulations", figure_path = temp_file, alt = paste0("Calibration plot"), width_percent = 100)
                )
                unlink(temp_file)
              }
            } else {
              # ggplot version is available
              temp_file <- tempfile(fileext = ".png")
              suppressMessages(suppressWarnings(ggsave(temp_file)))
              test_performance <- c(
                '<h2>Performance in simulations</h2>',
                figure_to_html_code(caption = "Performance in simulations", figure_path = temp_file, alt = paste0("Calibration plot"), width_percent = 100)
              )
              unlink(temp_file)
            }
          } else {
            actual_values <- actual
            calib_plot <- suppressWarnings(try(genCalCurve(actual_values, lp, family = ifelse(generic_input_parameters$outcome_count, "poisson", "gaussian")), silent = TRUE))
            if (TRUE %in% (class(calib_plot) == "try-error")) {
              # If this resulted in error use predtools. However this is based on untransformed lp. The backtransformation is applied at the stage of creating the calibration_plot_data
              calibration_plot_data <- cbind.data.frame(
                actual = actual_values,
                predicted = backtransform_lp(lp = lp, outcome_type = generic_input_parameters$outcome_type, outcome_count = generic_input_parameters$outcome_count)
                )
              calib_plot <- suppressWarnings(try(calibration_plot(data = calibration_plot_data, obs = "actual", pred = "predicted", title = "Calibration plot"), silent = TRUE))
              # If even this resulted in error, indicate this.
              if (TRUE %in% (class(calib_plot) == "try-error")) {
                test_performance <- c(
                  '<h2>Performance in simulations</h2>',
                  'Performance in simulations cannot be calculated because of an error in estimating the calibration curve.'
                )
              } else {
                temp_file <- tempfile(fileext = ".png")
                suppressMessages(suppressWarnings(ggsave(temp_file)))
                test_performance <- c(
                  '<h2>Performance in simulations</h2>',
                  figure_to_html_code(caption = "Performance in simulations", figure_path = temp_file, alt = paste0("Calibration plot"), width_percent = 100)
                )
                unlink(temp_file)
              }
            } else {
              # ggplot version is not available
              temp_file <- tempfile(fileext = ".png")
              png(temp_file)
              genCalCurve(actual_values, lp, family = ifelse(generic_input_parameters$outcome_count, "poisson", "gaussian"))
              recordPlot()
              dev.off()
              test_performance <- c(
                '<h2>Performance in simulations</h2>',
                figure_to_html_code(caption = "Performance in simulations", figure_path = temp_file, alt = paste0("Calibration plot"), width_percent = 100)
              )
              unlink(temp_file)
            }
          }
        } else {
          test_performance <- c(
            '<h2>Performance in simulations</h2>',
            'Performance in simulations cannot be calculated because linear predictor calculation is not applicable.'
          )
        }
      } else {
        # Use predtools. However this is based on untransformed lp. The backtransformation is applied at the stage of creating the calibration_plot_data
        actual <- results[[x]]$actual_predicted_results_apparent$actual_all_subjects
        lp <- results[[x]]$average_lp_all_subjects
        follow_up <- results[[x]]$actual_predicted_results_apparent$time_all_subjects
        # remove missing lp
        actual <- actual[! is.na(lp)]
        lp <- lp[! is.na(lp)]
        follow_up <- follow_up[! is.na(lp)]
        if (length(actual) > 0) {
          actual_values <- rep(NA, length(actual))
          actual_values[actual == levels(actual)[1]] <- 0
          actual_values[actual == levels(actual)[2]] <- 1
          calibration_plot_data <- cbind.data.frame(
            actual = actual_values,
            predicted = backtransform_lp(lp = lp, outcome_type = generic_input_parameters$outcome_type, outcome_count = generic_input_parameters$outcome_count),
            follow_up = follow_up)
          calib_plot <- suppressWarnings(try(calibration_plot(data = calibration_plot_data, obs = "actual", follow_up = "follow_up", pred = "predicted", x_lim=c(0, 1), title = "Calibration plot"), silent = TRUE))
          # If even this resulted in error, ignore follow-up time
          if (TRUE %in% (class(calib_plot) == "try-error")) {
            calib_plot <- suppressWarnings(try(calibration_plot(data = calibration_plot_data, obs = "actual", pred = "predicted", y_lim = c(0, 1), x_lim=c(0, 1), title = "Calibration plot"), silent = TRUE))
            # If even this caused an error indicate so
            if (TRUE %in% (class(calib_plot) == "try-error")) {
              test_performance <- c(
                '<h2>Performance in simulations</h2>',
                'Performance in simulations cannot be calculated because of an error in estimating the calibration curve.'
              )
            } else {
              # Indicate that the calibration plot based on rate could not be calculated and with a warning that the plot was based ignoring the follow-up time
              temp_file <- tempfile(fileext = ".png")
              suppressMessages(suppressWarnings(ggsave(temp_file)))
              test_performance <- c(
                '<h2>Performance in simulations</h2>',
                  'Performance in simulations could not be calculated based on rate. Therefore, follow-up time was removed. Therefore, the calibration plot is likely to be unreliable.<br>',
                figure_to_html_code(caption = "Performance in simulations", figure_path = temp_file, alt = paste0("Calibration plot"), width_percent = 100)
              )
              unlink(temp_file)
            }
          } else {
            temp_file <- tempfile(fileext = ".png")
            suppressMessages(suppressWarnings(ggsave(temp_file)))
            test_performance <- c(
              '<h2>Performance in simulations</h2>',
              figure_to_html_code(caption = "Performance in simulations", figure_path = temp_file, alt = paste0("Calibration plot"), width_percent = 100)
            )
            unlink(temp_file)
          }
        } else {
          test_performance <- c(
            '<h2>Performance in simulations</h2>',
            'Performance in simulations cannot be calculated because linear predictor calculation is not applicable.'
          )
        }
      }
    }
    # Combine the two ####
    html_calibration_plot <- c('<h2>Calibration plot</h2>', apparent_performance, test_performance)
    output <- {c(
      html_specific_input_parameters[[x]],
      if (! TRUE %in% (class(results[[x]]) == "try-error")) {
        results[[x]]$actual_predicted_results_apparent$html_file
      },
      html_calibration_plot
    )}
  }))
  relevant_parameters <- names(results[[1]]$optimism_corrected_performance)
  # Optimism across analyses
  html_optimism_all_analyses <- {c(
    "<h1>Optimism variability plots</h1>",
    "The plots below show the optimism variability for each model or score.",
    unlist(lapply(1:length(relevant_parameters), function(x) {
      box_plot_data <- do.call(rbind.data.frame, lapply(1:length(specific_input_parameters), function(z) {
        if (!TRUE %in% (class(results[[z]]) == "try-error")) {
          cbind.data.frame(
            group = names(specific_input_parameters[z]),
            optimism = results[[z]]$optimism[,relevant_parameters[x]]
          )
        }
      }))
      box_plot_data <- box_plot_data[! is.na(box_plot_data$optimism),]
      if (nrow(box_plot_data) > 0) {
        optimism_variability_plot <- ggplot(box_plot_data, aes(x = factor(group, levels = rev(unique(box_plot_data$group))), y = optimism, fill = factor(group))) +
          geom_boxplot() +
          coord_flip() +
          theme_minimal() +
          labs(
            title = "Variability in optimism",
            x = "Model or score",
            y = "Optimism") +
          theme_minimal() + theme(legend.position="none")
        temp_file <- tempfile(fileext = ".png")
        suppressWarnings(ggsave(temp_file, width = 6, height = length(unique(box_plot_data$group))/2, units = "in"))
        output <- figure_to_html_code(caption = relevant_parameters[x], figure_path = temp_file, alt = paste0("Optimism variability plot"), width_percent = 100)
        unlink(temp_file)
      } else {
        output <- "None of the groups had optimism data for this parameter"
      }
      return(output)
    }))
  )}
  # Performance parameters across all analyses ####
  performance_parameters <- {cbind.data.frame(
    parameters = c("accuracy", "c_statistic", "OE_ratio", "m_OE_ratio",
                   "calibration_intercept", "m_calibration_intercept",
                   "calibration_slope", "m_calibration_slope",
                   "error"),
    minimum = c(0, 0, NA, 0,
                NA, 0,
                NA, 0,
                0),
    maximum = c(1, 1, NA, 1,
                NA, NA,
                NA, NA,
                NA),
    higher_better = c(TRUE, TRUE, NA, FALSE,
                      NA, FALSE,
                      NA, FALSE,
                      FALSE),
    reference_line = c(0.5, 0.5, 1, 1,
                       0, 0,
                       1, 0,
                       0)
  )}
  performance_parameters_across_analyses <- do.call(rbind.data.frame, lapply(1:length(specific_input_parameters), function(x) {
    if (!TRUE %in% (class(results[[x]]) == "try-error")) {
      cbind.data.frame(analysis = names(specific_input_parameters)[x], lapply(results[[x]]$optimism_corrected_performance, function(z) {round(z,4)}))
    }
  }))
  performance_parameters_across_analyses_calibration_adjusted <- do.call(rbind.data.frame, lapply(1:length(specific_input_parameters), function(x) {
    if (!TRUE %in% (class(results[[x]]) == "try-error")) {
      if (specific_input_parameters[[x]]$develop_model == TRUE) {
        cbind.data.frame(analysis = names(specific_input_parameters)[x], lapply(results[[x]]$optimism_corrected_performance_calibration_adjusted, function(z) {round(z,4)}))
      }
    }
  }))
  performance_parameters_across_analyses_out_of_sample <- do.call(rbind.data.frame, lapply(1:length(specific_input_parameters), function(x) {
    if (!TRUE %in% (class(results[[x]]) == "try-error")) {
      cbind.data.frame(analysis = names(specific_input_parameters)[x], lapply(results[[x]]$out_of_sample_performance_summary, function(z) {round(z$est,4)}))
    }
  }))
  performance_est_ci <- lapply(1:length(relevant_parameters), function(x) {
    do.call(rbind.data.frame,
            lapply(1:length(specific_input_parameters), function(z) {
              if (!TRUE %in% (class(results[[z]]) == "try-error")) {
                cbind.data.frame(
                  analysis = names(specific_input_parameters)[z],
                  round(results[[z]]$optimism_corrected_performance_with_CI[[relevant_parameters[x]]],4)
                )
              }
            }))
  })
  names(performance_est_ci) <- relevant_parameters
  performance_est_ci_calibration_adjusted <- lapply(1:length(relevant_parameters), function(x) {
    do.call(rbind.data.frame,
            lapply(1:length(specific_input_parameters), function(z) {
              if (!TRUE %in% (class(results[[z]]) == "try-error")) {
                if (specific_input_parameters[[z]]$develop_model == TRUE) {
                  cbind.data.frame(
                    analysis = names(specific_input_parameters)[z],
                    round(results[[z]]$optimism_corrected_performance_with_CI_calibration_adjusted[[relevant_parameters[x]]],4)
                  )}
              }
            }))
  })
  names(performance_est_ci_calibration_adjusted) <- relevant_parameters
  performance_est_ci_out_of_sample <- lapply(1:length(relevant_parameters), function(x) {
    do.call(rbind.data.frame,
            lapply(1:length(specific_input_parameters), function(z) {
              if (!TRUE %in% (class(results[[z]]) == "try-error")) {
                cbind.data.frame(
                  analysis = names(specific_input_parameters)[z],
                  round(results[[z]]$out_of_sample_performance_summary[[relevant_parameters[x]]],4)
                )
              }
            }))
  })
  names(performance_est_ci_out_of_sample) <- relevant_parameters
  forest_plots <- unlist(lapply(relevant_parameters, function(x) {
    performance_parameters_row <- match(x, performance_parameters$parameters)
    est_ci_table <- performance_est_ci[[x]]
    est_ci_table <- est_ci_table[! is.na(est_ci_table$est),]
    if (nrow(est_ci_table) > 0) {
      estimate_with_ci_forest(
        est_ci_table = est_ci_table,
        parameter_name = x,
        minimum = performance_parameters$minimum[performance_parameters_row],
        maximum = performance_parameters$maximum[performance_parameters_row],
        reference_line = performance_parameters$reference_line[performance_parameters_row],
        convert_to_html_code = TRUE
      )
    }
  }))
  forest_plots_calibration_adjusted <- unlist(lapply(relevant_parameters, function(x) {
    performance_parameters_row <- match(x, performance_parameters$parameters)
    est_ci_table <- performance_est_ci_calibration_adjusted[[x]]
    est_ci_table <- est_ci_table[! is.na(est_ci_table$est),]
    if (nrow(est_ci_table) > 0) {
      estimate_with_ci_forest(
        est_ci_table = est_ci_table,
        parameter_name = x,
        minimum = performance_parameters$minimum[performance_parameters_row],
        maximum = performance_parameters$maximum[performance_parameters_row],
        reference_line = performance_parameters$reference_line[performance_parameters_row],
        convert_to_html_code = TRUE
      )
    }
  }))
  forest_plots_out_of_sample <- unlist(lapply(relevant_parameters, function(x) {
    performance_parameters_row <- match(x, performance_parameters$parameters)
    est_ci_table <- performance_est_ci_out_of_sample[[x]]
    est_ci_table <- est_ci_table[! is.na(est_ci_table$est),]
    if (nrow(est_ci_table) > 0) {
      estimate_with_ci_forest(
        est_ci_table = est_ci_table,
        parameter_name = x,
        minimum = performance_parameters$minimum[performance_parameters_row],
        maximum = performance_parameters$maximum[performance_parameters_row],
        reference_line = performance_parameters$reference_line[performance_parameters_row],
        convert_to_html_code = TRUE
      )
    }
  }))

  summary_results <- {list(
    performance_parameters_across_analyses = performance_parameters_across_analyses,
    performance_parameters_across_analyses_calibration_adjusted = performance_parameters_across_analyses_calibration_adjusted,
    performance_parameters_across_analyses_out_of_sample = performance_parameters_across_analyses_out_of_sample,
    performance_est_ci = performance_est_ci,
    performance_est_ci_calibration_adjusted = performance_est_ci_calibration_adjusted,
    performance_est_ci_out_of_sample = performance_est_ci_out_of_sample
  )}
  html_performance_parameters_across_all_analyses <- {c(
    '<h1>Performance</h1>',
    table_to_html_code(caption = "Performance", df = performance_parameters_across_analyses, width_percent = 100),
    forest_plots,
    '<h1>Calibration-adjusted performance</h1>',
    table_to_html_code(caption = "Performance of calibration-adjusted models", df = performance_parameters_across_analyses_calibration_adjusted, width_percent = 100),
    forest_plots_calibration_adjusted,
    '<h1>Out-of-sample performance</h1>',
    table_to_html_code(caption = "Performance in out of sample subjects", df = performance_parameters_across_analyses_out_of_sample, width_percent = 100),
    forest_plots_out_of_sample
  )}
  # Write to html ####
  html_file_location <- paste0(gsub("\\\\", "/", tempdir()), "/results.html")
  writeLines(paste0(unlist(c(html_file, html_each_analysis, html_optimism_all_analyses, html_performance_parameters_across_all_analyses)), collapse = "\n"),
             html_file_location
  )
  # Output ####
  output <- list(results = results, summary_results = summary_results, html_file_location = html_file_location)
}
