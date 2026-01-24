\name{compile_results}
\alias{compile_results}
\title{
Run the analysis and compile the results.
}
\description{
This is wrapper function that takes generic_input_parameters created using
\code{\link{create_generic_input_parameters}}, specific_input_parameters created
using \code{\link{create_specific_input_parameters}} and the datasets prepared using
\code{\link{prepare_datasets}}, and provides the results. For details, please see
below.
}
\usage{
compile_results(generic_input_parameters, specific_input_parameters,
prepared_datasets, verbose)
}
\arguments{
\item{generic_input_parameters}{This is a list that contains common information
across models. If one or more items are missing or incorrect, this may result in error.
Therefore, we recommend that you use the \code{\link{create_generic_input_parameters}}
function to create this input.}
\item{specific_input_parameters}{This is a list that contains information related
to each model or scoring system. If one or more items are missing or incorrect,
this may result in error. Therefore, we recommend that you use the
\code{\link{create_specific_input_parameters}}.}
\item{prepared_datasets}{Datasets prepared using the \code{\link{prepare_datasets}}.}
\item{verbose}{\emph{TRUE} if the progress must be displayed and \emph{FALSE} otherwise.}
}
\details{
\strong{Overview}
This is wrapper function that takes generic_input_parameters created using
\code{\link{create_generic_input_parameters}}, specific_input_parameters created
using \code{\link{create_specific_input_parameters}}, and the datasets prepared using
\code{\link{prepare_datasets}}, and provides the results. It uses the
\code{\link{perform_analysis}} function (which itself uses the
\code{\link{calculate_actual_predicted}} to develop and run models and
\code{\link{calculate_performance}} to calculate the performance). It also
creates calibration curves for each model and summarises the information across the
different scoring systems and models.

The following steps must be done sequentially. Please see example below which provides
the details of how to perform the analysis.
1. \emph{Process the data: }The dataset should be prepared correctly for the functions
to work correctly. This can be done using \code{\link{process_data}}.
2. \emph{Create the generic input parameters: }The generic input parameters should be
provided. You can check that the input parameters are correct using
\code{\link{create_generic_input_parameters}}.
3. \emph{Create the specific input parameters: }The specific input parameters should be
provided. You can check that the input parameters are correct using
\code{\link{create_specific_input_parameters}}.
4. \emph{Prepare the datasets: }The datasets for each simulation are prepared using
\code{\link{prepare_datasets}}.
5. \emph{Provide the correct parameters to this function: }Input 'generic_input_parameters',
'specific_input_parameters', and 'prepared_datasets' to obtain the results.

\strong{Preparing datasets for each simulation}
Please see \code{\link{prepare_datasets}}.

\strong{Calculation of actual and predicted values}
Please see \code{\link{calculate_actual_predicted}}.

\strong{Calculation of performance measures}
Please see \code{\link{calculate_performance}}.

\strong{Calculation of means and confidence intervals}
Please see \code{\link{perform_analysis}}.

\strong{Calibration curves}
Flexible calibration curves were created using the package
\code{\link[CalibrationCurves]{CalibrationCurves}}, based on the paper by
Van Calster et al, 2016.

The calibration curves have been described for only one simulation. It is impractical
to interpret each calibration curve (one for each simulation). Therefore, the calibration
curves are shown only for the apparent performance. As an experimental method,
the linear predictors were averaged across the simulations and calibration
curves for the average test performance are also provided.

When the flexible calibration curves resulted in errors or were not possible, for
example, the average linear predictors across simulations for time-to-event outcomes,
\code{\link[predtools]{calibration_plot}} was used to create the calibration plots.
Appropriate backtransformations were performed to reverse the transformations applied
while calculating the linear predictors. For further information about transformations,
please see \code{\link{calculate_actual_predicted}}.
}
\value{
\item{results}{These include the performance results for each model.}
\item{summary_results}{These include the summary performance results.}
\item{html_file_location}{This provides the location of the html file, which
contains the results in html format. This can be downloaded and the content copied to
word document or PowerPoint presentations.}
}
\references{
Collins GS, Dhiman P, Ma J, Schlussel MM, Archer L, Van Calster B, et al.
Evaluation of clinical prediction models (part 1): from development to external
validation. Bmj. 2024;384:e074819.

R Package "coxed". https://CRAN.R-project.org/package=coxed)

Van Calster B, Nieboer D, Vergouwe Y, De Cock B, Pencina MJ, Steyerberg EW. A
calibration hierarchy for risk models was defined: from utopia to empirical data.
Journal of Clinical Epidemiology. 2016;74:167-76.
}
\author{
Kurinchi Gurusamy
}
\seealso{
\code{\link{process_data}}
\code{\link{prepare_datasets}}
\code{\link{perform_analysis}}
\code{\link{calculate_actual_predicted}}
\code{\link{calculate_performance}}
\code{\link[CalibrationCurves]{CalibrationCurves}}
}
\examples{
  # Load packages ####
  library(base64enc)
  library(mime)
  library(pROC)
  library(survival)
  library(ggplot2)
  library(CalibrationCurves)
  library(predtools)
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
  # For the demonstration, only the first row is run. Please remove the line
  # below i.e., analysis_details <- analysis_details[1,]
  analysis_details <- analysis_details[1,]
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
  results <- compile_results(generic_input_parameters, specific_input_parameters,
                             prepared_datasets, verbose = FALSE)
  # Results html_file_location
  results$html_file_location
}

