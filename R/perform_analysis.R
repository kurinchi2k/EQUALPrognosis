perform_analysis <- function(generic_input_parameters, specific_input_parameters_each_analysis, prepared_datasets) {
  # Functions ####
  bca <- function(theta){
    if (TRUE %in% (! is.na(theta))) {
      theta <- theta[! is.na(theta)]
      theta <- theta[! is.nan(theta)]
      theta <- theta[! theta %in% c(Inf, -Inf)]
      sims <- length(theta)
      z.inv <- length(theta[theta < mean(theta)])/sims
      z <- qnorm(z.inv)
      U <- (sims - 1) * (mean(theta, na.rm=TRUE) - theta)
      top <- sum(U^3)
      under <- 6 * (sum(U^2))^{3/2}
      a <- top / under
      lower.inv <-  pnorm(z + (z + qnorm(0.025))/(1 - a * (z + qnorm(0.025))))
      lower <- quantile(theta, lower.inv, names=FALSE)
      upper.inv <-  pnorm(z + (z + qnorm(0.975))/(1 - a * (z + qnorm(0.975))))
      upper <- quantile(theta, upper.inv, names=FALSE)
    } else {
      lower = NA
      upper = NA
    }
    return(c(lower, upper))
  }
  calculate_modified_observed_expected_ratio <- function(OE_ratio) {
    if (OE_ratio > 1) {1/OE_ratio} else {OE_ratio}
  }
  calculate_modified_calibration_slope <- function(calibration_slope) {
    abs(1-calibration_slope)
  }
  calculate_modified_calibration_intercept <- function(calibration_intercept) {
    abs(calibration_intercept)
  }
  performance_parameters <- cbind.data.frame(
    parameters = c("accuracy", "c_statistic", "OE_ratio", "calibration_intercept", "calibration_slope", "error",
                   "m_OE_ratio", "m_calibration_intercept", "m_calibration_slope"),
    transformation = c("logit", "logit", "log", NA, NA, NA,
                       "log", NA, NA)
  )
  if(generic_input_parameters$outcome_type == "quantitative") {
    relevant_parameters <- performance_parameters[performance_parameters$parameters %in% c("calibration_intercept", "calibration_slope", "error", "m_calibration_intercept", "m_calibration_slope"),]
  } else if (generic_input_parameters$outcome_type == "binary") {
    relevant_parameters <- performance_parameters[performance_parameters$parameters %in% c("accuracy", "c_statistic", "OE_ratio", "calibration_intercept", "calibration_slope", "m_OE_ratio", "m_calibration_intercept", "m_calibration_slope"),]
  } else if (generic_input_parameters$outcome_type == "time-to-event") {
    relevant_parameters <- performance_parameters[performance_parameters$parameters %in% c("accuracy", "c_statistic", "OE_ratio", "calibration_slope", "m_OE_ratio", "m_calibration_slope"),]
  }
  # Create temporary folder for results storage for managing memory requirements ####
  create_temporary_folder <- function(simulations, simulations_per_file) {
    number_of_files = ceiling(simulations/simulations_per_file)
    simulation_start_end <- do.call(rbind.data.frame, lapply(1:number_of_files, function(x) {
      start_simulation <- (x-1)*simulations_per_file+1
      if (x < number_of_files) {
        end_simulation <- x*simulations_per_file
      } else {
        end_simulation <- simulations
      }
      output <- cbind.data.frame(start_simulation = start_simulation, end_simulation = end_simulation)
    }))
    temp_results_directory <- tempfile(pattern = "results_")
    placeholder <- dir.create(temp_results_directory)
    list(simulation_start_end = simulation_start_end, temp_results_directory = temp_results_directory)
  }
  placeholder <- {create_temporary_folder(
    simulations = generic_input_parameters$simulations,
    simulations_per_file = generic_input_parameters$simulations_per_file)}
  simulation_start_end <- placeholder$simulation_start_end
  temp_results_directory <- placeholder$temp_results_directory
  # Calculate the actual and predicted
  {
    actual_predicted_results_apparent <- {calculate_actual_predicted(
      prepared_datasets = prepared_datasets,
      outcome_name = generic_input_parameters$outcome_name,
      outcome_type = generic_input_parameters$outcome_type,
      outcome_time = generic_input_parameters$outcome_time,
      outcome_count = generic_input_parameters$outcome_count,
      develop_model = specific_input_parameters_each_analysis$develop_model,
      mandatory_predictors = specific_input_parameters_each_analysis$mandatory_predictors,
      optional_predictors = specific_input_parameters_each_analysis$optional_predictors,
      mandatory_interactions = specific_input_parameters_each_analysis$mandatory_interactions,
      optional_interactions = specific_input_parameters_each_analysis$optional_interactions,
      model_threshold_method = specific_input_parameters_each_analysis$model_threshold_method,
      scoring_system = specific_input_parameters_each_analysis$scoring_system,
      predetermined_threshold = specific_input_parameters_each_analysis$predetermined_threshold,
      higher_values_event = specific_input_parameters_each_analysis$higher_values_event,
      each_simulation = 1, bootstrap_sample = FALSE
    )}
    cat("\nBootstrap performance...")
    for (i in 1:nrow(simulation_start_end)) {
      results_each_batch <- lapply(simulation_start_end$start_simulation[i]:simulation_start_end$end_simulation[i], function(each_simulation) {
        calculate_actual_predicted(
          prepared_datasets = prepared_datasets,
          outcome_name = generic_input_parameters$outcome_name,
          outcome_type = generic_input_parameters$outcome_type,
          outcome_time = generic_input_parameters$outcome_time,
          outcome_count = generic_input_parameters$outcome_count,
          develop_model = specific_input_parameters_each_analysis$develop_model,
          mandatory_predictors = specific_input_parameters_each_analysis$mandatory_predictors,
          optional_predictors = specific_input_parameters_each_analysis$optional_predictors,
          mandatory_interactions = specific_input_parameters_each_analysis$mandatory_interactions,
          optional_interactions = specific_input_parameters_each_analysis$optional_interactions,
          model_threshold_method = specific_input_parameters_each_analysis$model_threshold_method,
          scoring_system = specific_input_parameters_each_analysis$scoring_system,
          predetermined_threshold = specific_input_parameters_each_analysis$predetermined_threshold,
          higher_values_event = specific_input_parameters_each_analysis$higher_values_event,
          each_simulation = each_simulation, bootstrap_sample = TRUE
        )
      })
      placeholder <- saveRDS(results_each_batch, paste0(temp_results_directory, "/batch_", i, ".rds"))
    }
    cat("\nCompiling the results from saved files and calculating the performance...")
    actual_predicted_results_bootstrap <- do.call(c, lapply(1:nrow(simulation_start_end), function(x) {
      readRDS(paste0(temp_results_directory, "/batch_", x, ".rds"))
    }))
    unlink(temp_results_directory, recursive = TRUE, force = TRUE)
  }
  # Calculate performance ####
  {
    apparent_performance <- {cbind.data.frame(
      performance = "apparent", simulation = NA,
      calculate_performance(
        outcome_type = generic_input_parameters$outcome_type,
        time = actual_predicted_results_apparent$time_all_subjects,
        outcome_count = generic_input_parameters$outcome_count,
        actual = actual_predicted_results_apparent$actual_all_subjects,
        predicted = actual_predicted_results_apparent$predicted_all_subjects,
        develop_model = specific_input_parameters_each_analysis$develop_model,
        lp = actual_predicted_results_apparent$lp_all_subjects
      ))}
    test_performance <- do.call(rbind.data.frame, lapply(1:length(actual_predicted_results_bootstrap), function(x) {
      {cbind.data.frame(
        performance = "test", simulation = x,
        calculate_performance(
          outcome_type = generic_input_parameters$outcome_type,
          time = actual_predicted_results_bootstrap[[x]]$time_all_subjects,
          outcome_count = generic_input_parameters$outcome_count,
          actual = actual_predicted_results_bootstrap[[x]]$actual_all_subjects,
          predicted = actual_predicted_results_bootstrap[[x]]$predicted_all_subjects,
          develop_model = specific_input_parameters_each_analysis$develop_model,
          lp = actual_predicted_results_bootstrap[[x]]$lp_all_subjects
        ))}
    }))
    bootstrap_performance <- do.call(rbind.data.frame, lapply(1:length(actual_predicted_results_bootstrap), function(x) {
      {cbind.data.frame(
        performance = "bootstrap", simulation = x,
        calculate_performance(
          outcome_type = generic_input_parameters$outcome_type,
          time = actual_predicted_results_bootstrap[[x]]$time_training,
          outcome_count = generic_input_parameters$outcome_count,
          actual = actual_predicted_results_bootstrap[[x]]$actual_training,
          predicted = actual_predicted_results_bootstrap[[x]]$predicted_training,
          develop_model = specific_input_parameters_each_analysis$develop_model,
          lp = actual_predicted_results_bootstrap[[x]]$lp_training
        ))}
    }))
    out_of_sample_performance <- do.call(rbind.data.frame, lapply(1:length(actual_predicted_results_bootstrap), function(x) {
      {cbind.data.frame(
        performance = "out_of_sample", simulation = x,
        calculate_performance(
          outcome_type = generic_input_parameters$outcome_type,
          time = actual_predicted_results_bootstrap[[x]]$time_only_validation,
          outcome_count = generic_input_parameters$outcome_count,
          actual = actual_predicted_results_bootstrap[[x]]$actual_only_validation,
          predicted = actual_predicted_results_bootstrap[[x]]$predicted_only_validation,
          develop_model = specific_input_parameters_each_analysis$develop_model,
          lp = actual_predicted_results_bootstrap[[x]]$lp_only_validation
        ))}
    }))
  }
  # Calculate the average lp of all_subjects from each simulation
  average_lp_all_subjects <- unlist(rowMeans(do.call(cbind.data.frame, lapply(actual_predicted_results_bootstrap, function(x) {x$lp_all_subjects})), na.rm = TRUE))
  # Calculate calibration-adjusted performance ####
  if (specific_input_parameters_each_analysis$develop_model == TRUE) {
    apparent_performance_calibration_adjusted <- {
      cbind.data.frame(
        performance = "apparent_calibration_adjusted", simulation = NA,calculate_performance(
          outcome_type = generic_input_parameters$outcome_type,
          time = actual_predicted_results_apparent$time_all_subjects,
          outcome_count = generic_input_parameters$outcome_count,
          actual = actual_predicted_results_apparent$actual_all_subjects,
          predicted = actual_predicted_results_apparent$predicted_all_subjects_calibration_adjusted,
          develop_model = specific_input_parameters_each_analysis$develop_model,
          lp = actual_predicted_results_apparent$lp_all_subjects_calibration_adjusted
        ))}
    test_performance_calibration_adjusted <- do.call(rbind.data.frame, lapply(1:length(actual_predicted_results_bootstrap), function(x) {
      {cbind.data.frame(
        performance = "test_calibration_adjusted", simulation = x,
        calculate_performance(
          outcome_type = generic_input_parameters$outcome_type,
          time = actual_predicted_results_bootstrap[[x]]$time_all_subjects,
          outcome_count = generic_input_parameters$outcome_count,
          actual = actual_predicted_results_bootstrap[[x]]$actual_all_subjects,
          predicted = actual_predicted_results_bootstrap[[x]]$predicted_all_subjects_calibration_adjusted,
          develop_model = specific_input_parameters_each_analysis$develop_model,
          lp = actual_predicted_results_bootstrap[[x]]$lp_all_subjects_calibration_adjusted
        ))}
    }))
    bootstrap_performance_calibration_adjusted <- do.call(rbind.data.frame, lapply(1:length(actual_predicted_results_bootstrap), function(x) {
      {cbind.data.frame(
        performance = "bootstrap_calibration_adjusted", simulation = x,
        calculate_performance(
          outcome_type = generic_input_parameters$outcome_type,
          time = actual_predicted_results_bootstrap[[x]]$time_training,
          outcome_count = generic_input_parameters$outcome_count,
          actual = actual_predicted_results_bootstrap[[x]]$actual_training,
          predicted = actual_predicted_results_bootstrap[[x]]$predicted_training_calibration_adjusted,
          develop_model = specific_input_parameters_each_analysis$develop_model,
          lp = actual_predicted_results_bootstrap[[x]]$lp_training_calibration_adjusted
        ))}
    }))
    out_of_sample_performance_calibration_adjusted <- do.call(rbind.data.frame, lapply(1:length(actual_predicted_results_bootstrap), function(x) {
      {cbind.data.frame(
        performance = "out_of_sample_calibration_adjusted", simulation = x,
        calculate_performance(
          outcome_type = generic_input_parameters$outcome_type,
          time = actual_predicted_results_bootstrap[[x]]$time_only_validation,
          outcome_count = generic_input_parameters$outcome_count,
          actual = actual_predicted_results_bootstrap[[x]]$actual_only_validation,
          predicted = actual_predicted_results_bootstrap[[x]]$predicted_only_validation_calibration_adjusted,
          develop_model = specific_input_parameters_each_analysis$develop_model,
          lp = actual_predicted_results_bootstrap[[x]]$lp_only_validation_calibration_adjusted
        ))}
    }))
  } else {
    apparent_performance_calibration_adjusted <- NA
    bootstrap_performance_calibration_adjusted <- NA
    test_performance_calibration_adjusted <- NA
    out_of_sample_performance_calibration_adjusted <- NA
  }
  # Calculate optimism ####
  optimism <- cbind.data.frame(
    performance = "optimism",
    simulation = bootstrap_performance$simulation,
    bootstrap_performance[,-match(c("performance", "simulation", "outcome", "message"), colnames(bootstrap_performance))] -
      test_performance[,-match(c("performance", "simulation", "outcome", "message"), colnames(test_performance))]
  )
  # Calculate optimism-corrected performance ####
  # The average optimism is calculated by applying transformation to a linear scale, i.e., logit transformation for probabiities and log transformation for ratios and reversing the transformation.
  average_optimism <- lapply(1:nrow(relevant_parameters), function(x) {
    placeholder <- optimism[,relevant_parameters$parameters[x]]
    placeholder <- placeholder[!is.na(placeholder)]
    if (length(placeholder) > 0) {
      if (! is.na(relevant_parameters$transformation[x])) {
        if (relevant_parameters$transformation[x] == "logit") {
          placeholder[placeholder <= 0] <- 0.000001
          placeholder[placeholder >= 1] <- 1- 0.000001
          placeholder <- qlogis(placeholder)
          optimism_summary <- cbind.data.frame(plogis(mean(placeholder)), t(plogis(bca(placeholder))))
        } else if (relevant_parameters$transformation[x] == "log") {
          placeholder[placeholder <= 0] <- 0.000001
          placeholder <- log(placeholder)
          optimism_summary <- cbind.data.frame(exp(mean(placeholder)), t(exp(bca(placeholder))))
        }
      } else {
        optimism_summary <- cbind.data.frame(mean(placeholder), t(bca(placeholder)))
      }
      colnames(optimism_summary) <- c("est", "lci", "uci")
    } else {
      optimism_summary <- cbind.data.frame(est = NA, lci = NA, uci = NA)
    }
    return(optimism_summary)
  })
  names(average_optimism) <- relevant_parameters$parameters
  optimism_corrected_performance <- lapply(1:nrow(relevant_parameters), function(x) {
    apparent_performance[1,relevant_parameters$parameters[x]] - average_optimism[[relevant_parameters$parameters[x]]]$est
  })
  names(optimism_corrected_performance) <- relevant_parameters$parameters
  optimism_corrected_performance_with_CI <- lapply(1:nrow(relevant_parameters), function(x) {
    placeholder <- apparent_performance[1,relevant_parameters$parameters[x]] - optimism[, relevant_parameters$parameters[x]]
    placeholder <- placeholder[!is.na(placeholder)]
    if (length(placeholder) > 0) {
      if (! is.na(relevant_parameters$transformation[x])) {
        if (relevant_parameters$transformation[x] == "logit") {
          placeholder[placeholder <= 0] <- 0.000001
          placeholder[placeholder >= 1] <- 1- 0.000001
          placeholder <- qlogis(placeholder)
          performance_summary <- cbind.data.frame(plogis(mean(placeholder)), t(plogis(bca(placeholder))))
        } else if (relevant_parameters$transformation[x] == "log") {
          placeholder[placeholder <= 0] <- 0.000001
          placeholder <- log(placeholder)
          performance_summary <- cbind.data.frame(exp(mean(placeholder)), t(exp(bca(placeholder))))
        }
      } else {
        performance_summary <- cbind.data.frame(mean(placeholder), t(bca(placeholder)))
      }
      colnames(performance_summary) <- c("est", "lci", "uci")
    } else {
      performance_summary <- cbind.data.frame(est = NA, lci = NA, uci = NA)
    }
    return(performance_summary)
  })
  names(optimism_corrected_performance_with_CI) <- relevant_parameters$parameters
  # Calculate out-of-sample performance summary ####
  out_of_sample_performance_summary <- lapply(1:nrow(relevant_parameters), function(x) {
    placeholder <- out_of_sample_performance[, relevant_parameters$parameters[x]]
    placeholder <- placeholder[!is.na(placeholder)]
    if (length(placeholder) > 0) {
      if (! is.na(relevant_parameters$transformation[x])) {
        if (relevant_parameters$transformation[x] == "logit") {
          placeholder[placeholder <= 0] <- 0.000001
          placeholder[placeholder >= 1] <- 1- 0.000001
          placeholder <- qlogis(placeholder)
          performance_summary <- cbind.data.frame(plogis(mean(placeholder)), t(plogis(bca(placeholder))))
        } else if (relevant_parameters$transformation[x] == "log") {
          placeholder[placeholder <= 0] <- 0.000001
          placeholder <- log(placeholder)
          performance_summary <- cbind.data.frame(exp(mean(placeholder)), t(exp(bca(placeholder))))
        }
      } else {
        performance_summary <- cbind.data.frame(mean(placeholder), t(bca(placeholder)))
      }
      colnames(performance_summary) <- c("est", "lci", "uci")
    } else {
      performance_summary <- cbind.data.frame(est = NA, lci = NA, uci = NA)
    }
    return(performance_summary)
  })
  names(out_of_sample_performance_summary) <- relevant_parameters$parameters
  # Calculate optimism-corrected calibration-adjusted performance ####
  if (specific_input_parameters_each_analysis$develop_model == TRUE) {
    # Calculate calibration-adjusted optimism ####
    optimism_calibration_adjusted <- cbind.data.frame(
      performance = "optimism_calibration_adjusted",
      simulation = bootstrap_performance_calibration_adjusted$simulation,
      bootstrap_performance_calibration_adjusted[,-match(c("performance", "simulation", "outcome", "message"), colnames(bootstrap_performance_calibration_adjusted))] -
        test_performance_calibration_adjusted[,-match(c("performance", "simulation", "outcome", "message"), colnames(test_performance_calibration_adjusted))]
    )
    # Calculate optimism-corrected calibration-adjusted performance ####
    # The average optimism is calculated by applying transformation to a linear scale, i.e., logit transformation for probabiities and log transformation for ratios and reversing the transformation.
    average_optimism_calibration_adjusted <- lapply(1:nrow(relevant_parameters), function(x) {
      placeholder <- optimism_calibration_adjusted[,relevant_parameters$parameters[x]]
      placeholder <- placeholder[!is.na(placeholder)]
      if (length(placeholder) > 0) {
        if (! is.na(relevant_parameters$transformation[x])) {
          if (relevant_parameters$transformation[x] == "logit") {
            placeholder[placeholder <= 0] <- 0.000001
            placeholder[placeholder >= 1] <- 1- 0.000001
            placeholder <- qlogis(placeholder)
            optimism_summary <- cbind.data.frame(plogis(mean(placeholder)), t(plogis(bca(placeholder))))
          } else if (relevant_parameters$transformation[x] == "log") {
            placeholder[placeholder <= 0] <- 0.000001
            placeholder <- log(placeholder)
            optimism_summary <- cbind.data.frame(exp(mean(placeholder)), t(exp(bca(placeholder))))
          }
        } else {
          optimism_summary <- cbind.data.frame(mean(placeholder), t(bca(placeholder)))
        }
        colnames(optimism_summary) <- c("est", "lci", "uci")
      } else {
        optimism_summary <- cbind.data.frame(est = NA, lci = NA, uci = NA)
      }
      return(optimism_summary)
    })
    names(average_optimism_calibration_adjusted) <- relevant_parameters$parameters
    optimism_corrected_performance_calibration_adjusted <- lapply(1:nrow(relevant_parameters), function(x) {
      apparent_performance_calibration_adjusted[1,relevant_parameters$parameters[x]] - average_optimism_calibration_adjusted[[relevant_parameters$parameters[x]]]$est
    })
    names(optimism_corrected_performance_calibration_adjusted) <- relevant_parameters$parameters
    optimism_corrected_performance_with_CI_calibration_adjusted <- lapply(1:nrow(relevant_parameters), function(x) {
      placeholder <- apparent_performance_calibration_adjusted[1,relevant_parameters$parameters[x]] - optimism_calibration_adjusted[, relevant_parameters$parameters[x]]
      placeholder <- placeholder[!is.na(placeholder)]
      if (length(placeholder) > 0) {
        if (! is.na(relevant_parameters$transformation[x])) {
          if (relevant_parameters$transformation[x] == "logit") {
            placeholder[placeholder <= 0] <- 0.000001
            placeholder[placeholder >= 1] <- 1- 0.000001
            placeholder <- qlogis(placeholder)
            performance_summary <- cbind.data.frame(plogis(mean(placeholder)), t(plogis(bca(placeholder))))
          } else if (relevant_parameters$transformation[x] == "log") {
            placeholder[placeholder <= 0] <- 0.000001
            placeholder <- log(placeholder)
            performance_summary <- cbind.data.frame(exp(mean(placeholder)), t(exp(bca(placeholder))))
          }
        } else {
          performance_summary <- cbind.data.frame(mean(placeholder), t(bca(placeholder)))
        }
        colnames(performance_summary) <- c("est", "lci", "uci")
      } else {
        performance_summary <- cbind.data.frame(est = NA, lci = NA, uci = NA)
      }
      return(performance_summary)
    })
    names(optimism_corrected_performance_with_CI_calibration_adjusted) <- relevant_parameters$parameters
    # Calculate out-of-sample calibration-adjusted performance ####
    out_of_sample_performance_summary_calibration_adjusted <- lapply(1:nrow(relevant_parameters), function(x) {
      placeholder <- out_of_sample_performance_calibration_adjusted[, relevant_parameters$parameters[x]]
      placeholder <- placeholder[!is.na(placeholder)]
      if (length(placeholder) > 0) {
        if (! is.na(relevant_parameters$transformation[x])) {
          if (relevant_parameters$transformation[x] == "logit") {
            placeholder[placeholder <= 0] <- 0.000001
            placeholder[placeholder >= 1] <- 1- 0.000001
            placeholder <- qlogis(placeholder)
            performance_summary <- cbind.data.frame(plogis(mean(placeholder)), t(plogis(bca(placeholder))))
          } else if (relevant_parameters$transformation[x] == "log") {
            placeholder[placeholder <= 0] <- 0.000001
            placeholder <- log(placeholder)
            performance_summary <- cbind.data.frame(exp(mean(placeholder)), t(exp(bca(placeholder))))
          }
        } else {
          performance_summary <- cbind.data.frame(mean(placeholder), t(bca(placeholder)))
        }
        colnames(performance_summary) <- c("est", "lci", "uci")
      } else {
        performance_summary <- cbind.data.frame(est = NA, lci = NA, uci = NA)
      }
      return(performance_summary)
    })
    names(out_of_sample_performance_summary_calibration_adjusted) <- relevant_parameters$parameters
  } else {
    optimism_calibration_adjusted <- NA
    average_optimism_calibration_adjusted <- NA
    optimism_corrected_performance_calibration_adjusted <- NA
    optimism_corrected_performance_with_CI_calibration_adjusted <- NA
    out_of_sample_performance_summary_calibration_adjusted <- NA
  }
  # Output ####
  output <- list(
    apparent_performance = apparent_performance,
    bootstrap_performance = bootstrap_performance,
    test_performance = test_performance,
    out_of_sample_performance = out_of_sample_performance,
    optimism = optimism,
    average_optimism = average_optimism,
    optimism_corrected_performance = optimism_corrected_performance,
    optimism_corrected_performance_with_CI = optimism_corrected_performance_with_CI,
    out_of_sample_performance_summary = out_of_sample_performance_summary,
    apparent_performance_calibration_adjusted = apparent_performance_calibration_adjusted,
    bootstrap_performance_calibration_adjusted = bootstrap_performance_calibration_adjusted,
    test_performance_calibration_adjusted = test_performance_calibration_adjusted,
    out_of_sample_performance_calibration_adjusted = out_of_sample_performance_calibration_adjusted,
    optimism_calibration_adjusted = optimism_calibration_adjusted,
    average_optimism_calibration_adjusted = average_optimism_calibration_adjusted,
    optimism_corrected_performance_calibration_adjusted = optimism_corrected_performance_calibration_adjusted,
    optimism_corrected_performance_with_CI_calibration_adjusted = optimism_corrected_performance_with_CI_calibration_adjusted,
    out_of_sample_performance_summary_calibration_adjusted = out_of_sample_performance_summary_calibration_adjusted,
    actual_predicted_results_apparent = actual_predicted_results_apparent,
    average_lp_all_subjects = average_lp_all_subjects
  )
}
