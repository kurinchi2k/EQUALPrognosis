\name{calculate_performance}
\alias{calculate_performance}
\title{
Calculate performance of prognostic models
}
\description{
This function calculates the different performance measures of prognostic models
and factors. Please see below for more details.
}
\usage{
calculate_performance(outcome_type, time, outcome_count, actual, predicted,
develop_model, lp)
}
\arguments{
\item{outcome_type}{One of \emph{'binary'}, \emph{'time-to-event'}, \emph{'quantitative'}.
Count outcomes are included in \emph{'quantitative'} outcome type and can be differentiated
from continuous outcomes by specifying outcome_count as \emph{TRUE}. Please see examples
below.}
\item{time}{Times at which the outcome was measured. This is applicable only
for \emph{'time-to-event'} outcome. For other outcome types, enter \emph{NA}.}
\item{outcome_count}{\emph{TRUE} if the outcome was a count outcome and \emph{FALSE}
otherwise.}
\item{actual}{A vector of actual values.}
\item{predicted}{A vector of predicted values.}
\item{develop_model}{\emph{TRUE}, if you a model was developed; \emph{FALSE},
if a scoring system with a predetermined threshold (if applicable) was used.}
\item{lp}{A vector of linear predictors (applicable only if you have developed a model.)}
}
\details{
\strong{General comment}
Most of the input parameters are already available from the generic and specific input
parameters created using \code{\link{create_generic_input_parameters}} and
\code{\link{create_specific_input_parameters}}. This function is used by the
\code{\link{compile_results}} function which provides the correct input parameters
based on the entries while using \code{\link{create_generic_input_parameters}} and
\code{\link{create_specific_input_parameters}} and the output from
\code{\link{calculate_actual_predicted}}.

\strong{Performance measures}
The performance was measured by the following parameters.
\strong{Accuracy}
Number of correct predictions/number of participants in whom the predictions
were made (Rainio et al., 2024).
\strong{Calibration}
Three measures of calibration are used.
\emph{Observed/expected ratio}
Please see Riley et al., 2024.
Values closer to 1 are better; ratios < 1 indicate overestimation of risk by the
model while ratios > 1 indicate underestimation of risk by the model
(Riley et al., 2024).

We treated underestimation and overestimation equally, i.e.,
an observed-expected ratio of 0.8 was considered equivalent to 1/0.8 = 1.25.
Therefore, we converted the observed-expected ratios to be in the same direction
(less than 1) \emph{('modified observed-expected ratio')}. This ensured that while
calculating the test performance and bootstrap performance, lower numbers
consistently indicated worse test performance (as they are more distant from 1)
and higher numbers consistently indicated better performance (noting that the
maximum value of the modified observed-expected ratio was 1). This modification
also helps in interpretation of comparison for different models, some of which may
overestimate the risk while others might underestimate the risk.

For assessing the calibration, when the expected events were zero, 0.5 were
added to both the observed events and expected events.

\emph{Calibration intercept} and \emph{calibration slope}:
\emph{Calibration slope} quantifies the spread of the risk probabilities in relation to
the observed events (Stevens et al., 2020). We used the methods described by
Riley et al, 2024 to calculate the calibration intercept and slope for all outcomes
other than time-to-event outcomes. Essentially, this involves the following
regression equation: \emph{Y = calibration intercept + coefficient * linear predictor},
where \emph{'Y'} is the log odds of observed event, log risk of observed event,
and the untransformed outcomes for binary, count, and continuous outcomes respectively.

Estimation in time-to-event is lot more uncertain and should be considered
experimental. Please note that that Cox regression does not have an intercept separately,
as the intercept is included in the baseline hazard (SAS Support, 2017).

Values closer to 1 indicate better performance when the intercept is close to 0;
values further away from 1 indicate that the predictions are incorrect in some ranges
(Van Calster et al., 2019; Riley et al., 2024; Stevens e al., 2020). The further
away from 1, the worse the relationship between the log odds of observed event,
log hazard, log risk of observed event, and the untransformed outcome with the
linear predictor.

To allow easy comparison with lower values indicating closer
to 1 and higher values indicating further away from 1, this function also calculates the
\emph{'modified calibration slope'} using the following formula:
\emph{'Modified calibration slope = absolute value (1-calibration slope)'}.

\emph{Calibration intercept} (also called "calibration-in-the-large") in the
calibration regression equation evaluates whether the observed event proportion
equals the average predicted risk (Van Calster et al, 2019).

This function also calculates the \emph{'modified calibration intercept'} as
the absolute value of calibration intercept to allow lower values indicating
closer to 0 and higher values indicating further away from 0.

\emph{C-statistic}
This is the area under the ROC curve and a measure of discrimination (Riley et al.
2025). This was calculated using \code{\link[pROC]{roc}}. Higher values indicate
better performance.
}
\value{
\item{output}{A dataframe of the calculated evaluation parameters}
}
\references{
Rainio O, Teuho J, Klén R. Evaluation metrics and statistical tests for machine
learning. Scientific Reports. 2024;14(1):6086.

Riley RD, Archer L, Snell KIE, Ensor J, Dhiman P, Martin GP, et al. Evaluation
of clinical prediction models (part 2): how to undertake an external validation
study. BMJ. 2024;384:e074820.

SAS Support. https://support.sas.com/kb/24/457.html (accessed on 16 January 2026).

Stevens RJ, Poppe KK. Validation of clinical prediction models: what does the
"calibration slope" really measure? Journal of Clinical Epidemiology. 2020;118:93-9.

Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW,
Bossuyt P, et al. Calibration: the Achilles heel of predictive analytics.
BMC Medicine. 2019;17(1):230.
}
\author{
Kurinchi Gurusamy
}
\seealso{
\code{\link[pROC]{roc}}
}
\examples{
\dontrun{
library(survival)
colon$status <- factor(as.character(colon$status))
generic_input_parameters <- create_generic_input_parameters(
  general_title = "Prediction of colon cancer death", simulations = 2000,
  simulations_per_file = 100, seed = 1, df = colon, outcome_name = "status",
  outcome_type = "time-to-event", outcome_time = "time", outcome_count = FALSE)
analysis_details <- cbind.data.frame(
  name = c('age', 'single_mandatory_predictor', 'complex_models',
           'complex_models_only_optional_predictors'),
  analysis_title = c('Simple cut-off based on age', 'Single mandatory predictor (rx)',
                     'Multiple mandatory and optional predictors',
                     'Multiple optional predictors only'),
  develop_model = c(FALSE, TRUE, TRUE, TRUE),
  mandatory_predictors = c(NA, 'rx', 'rx; differ; perfor; adhere; extent', NA),
  optional_predictors = c(NA, NA, 'sex; age; nodes', 'rx; differ; perfor'),
  mandatory_interactions = c(NA, NA, 'rx; differ; extent', NA),
  optional_interactions = c(NA, NA, 'perfor; adhere; sex; age; nodes', 'rx; differ'),
  model_threshold_method = c(NA, 'youden', 'youden', 'youden'),
  scoring_system = c('age', NA, NA, NA),
  predetermined_threshold = c('60', NA, NA, NA),
  higher_values_event = c(TRUE, NA, NA, NA)
)
write.csv(analysis_details, paste0(tempdir(), "/analysis_details.csv"),
          row.names = FALSE, na = "")
analysis_details_path <- paste0(tempdir(), "/analysis_details.csv")
specific_input_parameters <- create_specific_input_parameters(
  generic_input_parameters = generic_input_parameters,
  analysis_details_path = analysis_details_path)
prepared_datasets <- {prepare_datasets(
  df = generic_input_parameters$df,
  simulations = generic_input_parameters$simulations,
  seed = generic_input_parameters$seed,
  outcome_name = generic_input_parameters$outcome_name,
  outcome_type = generic_input_parameters$outcome_type,
  outcome_time = generic_input_parameters$outcome_time)
}
specific_input_parameters_each_analysis <- specific_input_parameters[[1]]
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
  bootstrap_results <- lapply(1:generic_input_parameters$simulations,
  function(each_simulation) {
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
}
}
