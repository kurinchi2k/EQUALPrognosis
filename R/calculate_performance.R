calculate_performance <- function(outcome_type, time, outcome_count, actual, predicted, develop_model = TRUE, lp) {
  calculate_accuracy <- function(actual, predicted) {
    actual <- as.factor(actual)
    predicted <- as.factor(predicted)
    contingency_table <- as.data.frame.matrix(table(actual, predicted))
    # Standard calculations only when there are two rows and columns (indicating that there are events and no events in both the actual and predicted)
    if ((nrow(contingency_table) > 1) & (ncol(contingency_table) > 1)) {
      correct <- contingency_table[1,1] + contingency_table[2,2]
      FP <- contingency_table[1,2]
      FN <- contingency_table[2,1]
      accuracy <- correct/(correct + FP + FN)
    } else {
      if ((nlevels(actual) == 1) & (nlevels(predicted) == 1)) {
        # If there is only one level in both
        if (is.na(match(unique(actual), levels(predicted)))) {
          # If no match, accuracy = 0, otherwise it is perfect match
          accuracy <- 0
        } else {
          accuracy <- 1
        }
        # Otherwise either the actual have only one level or the predicted have only one level
      } else if (nlevels(actual) == 1) {
        # If there is only one level for actual, find if all the actuals are event or no event
        if (match(unique(actual), levels(predicted))[! is.na(match(unique(actual), levels(predicted)))] == 1) {
          # If it is all no events, no one has the event but some are predicted to have the event, i.e, only false positives
          correct = contingency_table[1,1]
          FP = contingency_table[1,2]
          FN = 0
        } else {
          # If it is all events, everyone has the event, but some are predicted to have no event i.e., only false negatives
          correct = contingency_table[1,2]
          FP = 0
          FN = contingency_table[1,1]
        }
        accuracy <- correct/(correct + FP + FN)
      } else if (nlevels(predicted) == 1) {
        # If number of levels of predicted is 1, find if all the predictions are event or no event
        if (match(unique(predicted), levels(actual))[! is.na(match(unique(predicted), levels(actual)))] == 1) {
          # If no one is predicted to have the event, there are only correct and false negatives
          correct = contingency_table[1,1]
          FP = 0
          FN = contingency_table[2,1]
        } else {
          # If everyone is predicted to have the event, there are only correct only false positives
          correct = contingency_table[2,1]
          FP = contingency_table[1,1]
          FN = 0
        }
        accuracy <- correct/(correct + FP + FN)
      }
    }
    return(accuracy)
  }
  calculate_c_statistic <- function(actual, predicted) {
    actual <- as.numeric(as.factor(actual))
    predicted <- as.numeric(as.factor(predicted))
    # This uses the pROC package to obtain the c-statistic. C-statistic is provided directly
    AUROC <- try(pROC::roc(actual,predicted,ci = FALSE, quiet = TRUE), silent = TRUE)
    # When there is a perfect or no match between actual and predicted, it can result in error. The AUROC when there is a perfect match is 1 (or 100%) and when there is no match is 0 (or 0%)
    if (TRUE %in% (class(AUROC) == "try-error")) {
      if (identical(actual, predicted)) {
        c_statistic <- 1
      } else {
        c_statistic <- 0
      }
    } else {
      c_statistic <- AUROC$auc
    }
    return(c_statistic)
  }
  calculate_error <- function(actual, predicted) {
    errors <- abs(as.numeric(actual) - as.numeric(predicted))
    mean_error <- mean(errors, na.rm = TRUE)
    return(mean_error)
  }
  calculate_observed_expected_ratio <- function(actual, predicted) {
    observed = length(actual[
      (! is.na(actual)) &
        (actual == levels(actual)[2])
    ])
    expected = length(predicted[
      (! is.na(predicted)) &
        (predicted == levels(predicted)[2])
    ])
    if ((observed == 0) | (expected == 0)) {
      observed = observed + 0.5
      expected = expected + 0.5
    }
    OE_ratio <- observed/expected
    return(OE_ratio)
  }
  create_intercept_slope_assessment_text <- function(outcome_type, outcome_count) {
    if (outcome_type %in% c("binary", "time-to-event", "quantitative")) {
      output <- {paste0(
        ifelse(outcome_type == 'time-to-event',
               'coxph(Surv(time, actual)',
               'glm(actual'
        ),
        ' ~ back_transformed_lp' ,
        ', data = df',
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
    } else {
      output <- "Unrecognised outcome_type. Accepted outcome_types are: 'binary', 'time-to-event', 'quantitative'"
    }
    return(output)
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
  outcome <- "Successful"
  message <- "All analyses performed as planned."
  if (develop_model == TRUE) {
    intercept_slope_text <- create_intercept_slope_assessment_text(outcome_type = outcome_type, outcome_count = outcome_count)
    # Back transformation
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
    df <- cbind.data.frame(actual = actual, predicted = predicted, back_transformed_lp = back_transformed_lp)
    if (outcome_type == "time-to-event") {
      df <- cbind.data.frame(df, time = time)
      df$actual <- as.numeric(df$actual)
    }
    intercept_slope <- suppressWarnings(try(eval(parse(text = intercept_slope_text)), silent = TRUE))
    if (TRUE %in% (class(intercept_slope)) == "try-error") {
      calibration_intercept <- NA
      calibration_slope <- NA
      m_calibration_intercept <- NA
      m_calibration_slope <- NA
      outcome <- "Unsuccessful"
      message <- "The calibration model did not converge."
    } else {
      calibration_coef <- cbind.data.frame(
        names = names(coef(intercept_slope)),
        coef = coef(intercept_slope)
      )
      if (outcome_type == "time-to-event") {
        calibration_intercept <- NA
        m_calibration_intercept <- NA
      } else {
        calibration_intercept <- calibration_coef$coef[calibration_coef$names == "(Intercept)"]
        m_calibration_intercept <- calculate_modified_calibration_intercept(calibration_intercept)
      }
      calibration_slope <- calibration_coef$coef[calibration_coef$names == "back_transformed_lp"]
      m_calibration_slope <- calculate_modified_calibration_slope(calibration_slope)
    }
  } else {
    calibration_intercept <- NA
    calibration_slope <- NA
    m_calibration_intercept <- NA
    m_calibration_slope <- NA
  }
  # Calculate other performance ####
  if (outcome_type == "quantitative") {
    error <- calculate_error(actual = actual, predicted = predicted)
    accuracy <- NA
    c_statistic <- NA
    OE_ratio <- NA
    m_OE_ratio <- NA
  } else {
    error <- NA
    accuracy <- calculate_accuracy(actual = actual, predicted = predicted)
    c_statistic <- calculate_c_statistic(actual = actual, predicted = predicted)
    OE_ratio <- calculate_observed_expected_ratio(actual = actual, predicted = predicted)
    m_OE_ratio <- calculate_modified_observed_expected_ratio(OE_ratio)
  }
  # Output ####
  output <- {list(
    outcome = outcome,
    message = message,
    calibration_intercept = calibration_intercept,
    calibration_slope = calibration_slope,
    m_calibration_intercept = m_calibration_intercept,
    m_calibration_slope = m_calibration_slope,
    error = error,
    accuracy = accuracy,
    c_statistic = c_statistic,
    OE_ratio = OE_ratio,
    m_OE_ratio = m_OE_ratio
  )}
  output <- cbind.data.frame(output)
  return(output)
}
