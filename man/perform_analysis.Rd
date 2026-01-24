\name{perform_analysis}
\alias{perform_analysis}
\title{
Perform analysis
}
\description{
This uses the \code{\link{calculate_actual_predicted}} to develop and run models and
\code{\link{calculate_performance}} to calculate the mean and confidence intervals
of the performance (please see details below).
}
\usage{
perform_analysis(generic_input_parameters,
specific_input_parameters_each_analysis, prepared_datasets, verbose)
}
\arguments{
\item{generic_input_parameters}{This is a list that contains common information
across models. If one or more items are missing or incorrect, this may result in error.
Therefore, we recommend that you use the \code{\link{create_generic_input_parameters}}
function to create this input.}
\item{specific_input_parameters_each_analysis}{This corresponds to each analysis, i.e.,
a model or scoring system. If one or more items are missing or incorrect,
this may result in error. Therefore, we recommend that you use the
\code{\link{create_specific_input_parameters}}.}
\item{prepared_datasets}{Datasets prepared using the
\code{\link{prepare_datasets}}.}
\item{verbose}{\emph{TRUE} if the progress must be displayed and \emph{FALSE} otherwise.}
}
\details{
\strong{Preparing datasets for each simulation}
Please see \code{\link{prepare_datasets}}.

\strong{Calculation of actual and predicted values}
Please see \code{\link{calculate_actual_predicted}}, particulary for details of
apparent performance, bootstrap performance, test performance, optimism as described
by Collins et al, 2024.

\strong{Calculation of performance measures}
Please see \code{\link{calculate_performance}}.

\strong{Calculation of means and confidence intervals}
For calculating the average performance measures and their confidence intervals
across multiple simulations, appropriate transformations were performed first. After
this, the bias-corrected accelerated confidence intervals were calculated based on
the "bca" function from coxed package, which is not maintained anymore (R, 2025).
The bias-corrected accelerated confidence intervals of the transformed data were
then back transformed.

The "enhanced bootstrapping internal validation approach" method described by Collins
et al., 2024 provides only the mean optimism-corrected performance. However, we have
optimism from multiple simulations. Therefore, rather than calculating the average
and then subtracting it from the apparent performance, the optimism from each
simulation was subtracted from the apparent performance. This allowed calculation of
the confidence intervals of the optimism-corrected performance using the bca function
(after appropriate transformation).

The performance measures of the calibration intercept-slope adjusted models were also
assessed by the same method. We have also presented the performance of the models in the
'out-of-sample subjects', i.e., the subjects who were not included in the bootstrap sample.
}
\value{
\item{apparent_performance}{Model is developed in the entire dataset and performance
evaluated in the same sample.}
\item{bootstrap_performance}{Model is developed in a subset of data (training set) and
evaluated in the training dataset}
\item{test_performance}{Model developed in the training set is evaluated in the entire
dataset.}
\item{out_of_sample_performance}{Performance in the sample that was not included in the
training dataset}
\item{optimism}{Test performance - bootstrap performance}
\item{average_optimism}{Average of the optimism}
\item{optimism_corrected_performance}{Apparent performance - average optimism}
\item{optimism_corrected_performance_with_CI}{Please see details above.}
\item{out_of_sample_performance_summary}{Please see details above.}
\item{apparent_performance_calibration_adjusted}{For details of calibration adjustment
see \code{\link{calculate_actual_predicted}}}
\item{bootstrap_performance_calibration_adjusted}{As above}
\item{test_performance_calibration_adjusted}{As above}
\item{out_of_sample_performance_calibration_adjusted}{As abovee}
\item{optimism_calibration_adjusted}{As above}
\item{average_optimism_calibration_adjusted}{As above}
\item{optimism_corrected_performance_calibration_adjusted}{As above}
\item{optimism_corrected_performance_with_CI_calibration_adjusted}{As above}
\item{out_of_sample_performance_summary_calibration_adjusted}{Summary of out-of-sample
performance}
\item{apparent_performance_adjusted_mandatory_predictors_only}{For details of this
model, used only for research purposes, see \code{\link{calculate_actual_predicted}},
section, 'Model with with only the mandatory predictors but based on the coefficients
of the entire model'.}
\item{bootstrap_performance_adjusted_mandatory_predictors_only}{As above}
\item{test_performance_adjusted_mandatory_predictors_only}{As above}
\item{out_of_sample_performance_adjusted_mandatory_predictors_only}{As abovee}
\item{optimism_adjusted_mandatory_predictors_only}{As above}
\item{average_optimism_adjusted_mandatory_predictors_only}{As above}
\item{optimism_corrected_performance_adjusted_mandatory_predictors_only}{As above}
\item{optimism_corrected_performance_with_CI_adjusted_mandatory_predictors_only}{As above}
\item{out_of_sample_performance_summary_adjusted_mandatory_predictors_only}{Summary of
out-of-sample performance}
\item{actual_predicted_results_apparent}{Output from
\code{\link{calculate_actual_predicted}} retained for some later calculations.}
\item{average_lp_all_subjects}{Output from \code{\link{calculate_actual_predicted}}
retained for some later calculations.}
}
\references{
Collins GS, Dhiman P, Ma J, Schlussel MM, Archer L, Van Calster B, et al.
Evaluation of clinical prediction models (part 1): from development to external
validation. Bmj. 2024;384:e074819.
}
\author{
Kurinchi Gurusamy
}
\seealso{
\code{\link{prepare_datasets}}
\code{\link{calculate_actual_predicted}}
\code{\link{calculate_performance}}
}
\examples{
  library(survival)
  colon$status <- factor(as.character(colon$status))
  # For testing, only 5 simulations are used here. Usually at least 300 to 500
  # simulations are a minimum. Increasing the simulations leads to more reliable results.
  # The default value of 2000 simulations should provide reasonably reliable results.
  generic_input_parameters <- create_generic_input_parameters(
    general_title = "Prediction of colon cancer death", simulations = 5,
    simulations_per_file = 20, seed = 1, df = colon, outcome_name = "status",
    outcome_type = "time-to-event", outcome_time = "time", outcome_count = FALSE,
    verbose = FALSE)$generic_input_parameters
  analysis_details <- cbind.data.frame(
    name = c('age', 'single_mandatory_predictor', 'complex_models',
             'complex_models_only_optional_predictors', 'predetermined_model_text'),
    analysis_title = c('Simple cut-off based on age', 'Single mandatory predictor (rx)',
                       'Multiple mandatory and optional predictors',
                       'Multiple optional predictors only', 'Predetermined model text'),
    develop_model = c(FALSE, TRUE, TRUE, TRUE, TRUE),
    predetermined_model_text = c(NA, NA, NA, NA,
    "cph(Surv(time, status) ~ rx * age, data = df_training_complete, x = TRUE, y = TRUE)"),
    mandatory_predictors = c(NA, 'rx', 'rx; differ; perfor; adhere; extent', NA, "rx; age"),
    optional_predictors = c(NA, NA, 'sex; age; nodes', 'rx; differ; perfor', NA),
    mandatory_interactions = c(NA, NA, 'rx; differ; extent', NA, NA),
    optional_interactions = c(NA, NA, 'perfor; adhere; sex; age; nodes', 'rx; differ', NA),
    model_threshold_method = c(NA, 'youden', 'youden', 'youden', 'youden'),
    scoring_system = c('age', NA, NA, NA, NA),
    predetermined_threshold = c('60', NA, NA, NA, NA),
    higher_values_event = c(TRUE, NA, NA, NA, NA)
  )
  write.csv(analysis_details, paste0(tempdir(), "/analysis_details.csv"),
            row.names = FALSE, na = "")
  analysis_details_path <- paste0(tempdir(), "/analysis_details.csv")
  # verbose is TRUE as default. If you do not want the outcome displayed, you can
  # change this to FALSE, as shown here
  results <- create_specific_input_parameters(
    generic_input_parameters = generic_input_parameters,
    analysis_details_path = analysis_details_path, verbose = FALSE)
  specific_input_parameters <- results$specific_input_parameters
  # Set a seed for reproducibility - Please see details above
  set.seed(generic_input_parameters$seed)
  prepared_datasets <- {prepare_datasets(
    df = generic_input_parameters$df,
    simulations = generic_input_parameters$simulations,
    outcome_name = generic_input_parameters$outcome_name,
    outcome_type = generic_input_parameters$outcome_type,
    outcome_time = generic_input_parameters$outcome_time,
    verbose = FALSE)}
  # There is no usually no requirement to call this function directly. This is used
  # by the perform_analysis function to create the actual and predicted values.
  specific_input_parameters_each_analysis <- specific_input_parameters[[1]]
  results <- perform_analysis(generic_input_parameters,
  specific_input_parameters_each_analysis, prepared_datasets, verbose = FALSE)
  results$apparent_performance
}

