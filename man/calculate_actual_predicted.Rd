\name{calculate_actual_predicted}
\alias{calculate_actual_predicted}
\title{
Calculate actual and predicted values
}
\description{
This takes the datasets, prepared using the \code{\link{prepare_datasets}}
function and the input parameters described below, creates a model, runs the
model, and predicts the outcome.
}
\usage{
calculate_actual_predicted(prepared_datasets, outcome_name, outcome_type,
outcome_time, outcome_count, develop_model, predetermined_model_text,
mandatory_predictors, optional_predictors, mandatory_interactions,
optional_interactions, model_threshold_method, scoring_system,
predetermined_threshold, higher_values_event, each_simulation,
bootstrap_sample, verbose)
}
\arguments{
\item{prepared_datasets}{Datasets prepared using the \code{\link{prepare_datasets}}.}
\item{outcome_name}{Name of the colummn that contains the outcome data. This
must be a column name in the \emph{'df'} provided as input.}
\item{outcome_type}{One of \emph{'binary'}, \emph{'time-to-event'}, \emph{'quantitative'}.
Count outcomes are included in \emph{'quantitative'} outcome type and can be differentiated
from continuous outcomes by specifying outcome_count as \emph{TRUE}. Please see examples
below.}
\item{outcome_time}{The name of the column that provides the follow-up time. This
is applicable only for \emph{'time-to-event'} outcome. For other outcome types,
enter \emph{NA}.}
\item{outcome_count}{\emph{TRUE} if the outcome was a count outcome and \emph{FALSE}
otherwise.}
\item{develop_model}{\emph{TRUE}, if you want to develop a model; \emph{FALSE},
if you want to use a scoring system with a predetermined threshold (if applicable).}
\item{predetermined_model_text}{You can create the model text from the mandatory
and optional predictors and interactions or finer control of the model, you can provide
the model text directly.}
\item{mandatory_predictors}{Predictors that must be included in the model. These
should be provided even if you provide the 'predetermined_model_text'.}
\item{optional_predictors}{Optional predictors that may be included in the model
by \code{\link{step}}. These should be provided even if you provide the
'predetermined_model_text'.}
\item{mandatory_interactions}{Interactions that must be included in the model. These
should be provided even if you provide the 'predetermined_model_text'.}
\item{optional_interactions}{Optional interactions that may be included in the
model by \code{\link{step}}. These should be provided even if you provide the
'predetermined_model_text'.}
\item{model_threshold_method}{One of \emph{'youden'}, \emph{'topleft'},
\emph{'heuristic'}. Please see description below.}
\item{scoring_system}{Name of the pre-existing scoring system. This is ignored
if \emph{develop_model} is \emph{TRUE}.}
\item{predetermined_threshold}{Pre-determined threshold of the pre-existing
scoring system. This is mandatory when \emph{develop_model} is \emph{FALSE} and when the
\emph{outcome_type} is \emph{'binary'} or \emph{'time-to-event'}. This is ignored if
\emph{develop_model} is \emph{TRUE} or when the \emph{outcome_type} is
\emph{'quantitative'}.}
\item{higher_values_event}{\emph{TRUE} if higher values of the pre-existing system
indicates event and \emph{FALSE} otherwise. This is mandatory when \emph{develop_model} is
\emph{FALSE} and when the \emph{outcome_type} is \emph{'binary'} or \emph{'time-to-event'}.
This is ignored if \emph{develop_model} is \emph{TRUE} or when the \emph{outcome_type} is
\emph{'quantitative'}.}
\item{each_simulation}{The number of the simulation in the prepared datasets.
Please see \code{\link{prepare_datasets}}.}
\item{bootstrap_sample}{TRUE if you are calculating the bootstrap and test
performance and FALSE if you are calculating the apparent performance. Please
see below and Collins et al, 2024.}
\item{verbose}{\emph{TRUE} if the progress must be displayed and \emph{FALSE} otherwise.}
}
\details{
\strong{General comment}
Most of the input parameters are already available from the generic and specific input
parameters created using \code{\link{create_generic_input_parameters}} and
\code{\link{create_specific_input_parameters}}. This function is used by the
\code{\link{perform_analysis}} function which provides the correct input parameters
based on the entries provided in \code{\link{create_generic_input_parameters}} and
\code{\link{create_specific_input_parameters}}.

\strong{Overview}
This is a form of enhanced bootstrapping internal validation approach to
calculate the optimism-corrected performance measures described by Collins et al.,
2024. This involves calculating the apparent performance by developing the model
in the entire dataset, repeated sampling with replacement (bootstrap sample),
evaluating the performance of the model in each simulation of the bootstrap
sample (bootstrap performance), evaluating the performance of the model
(developed in the bootstrap sample of each simulation) on the ‘test sample’
i.e., all the subjects in the dataset from which the bootstrap sample was
obtained (test performance), calculating the optimism as the difference between
bootstrap performance and test performance in each simulation,
calculating the average optimism, and finally subtracting the average optimism
from the apparent performance to calculate the optimism-corrected performance
measures (Collins et al, 2024).

The model development is performed using \code{\link[stats]{glm}} for all outcomes
other than time-to-event outcomes and \code{\link[survival]{coxph}} for time-to-event
outcomes. You can either provide a model text that you have developed for finer
control of the interactions to be considered or included or you can let the computer
build the model text based on the mandatory and optional predictors and interactions.

\strong{Linear predictors}
The linear predictor describes the relationship between the outcome and the
predictors, and is a function of the covariate (predictor) values and
coefficients of the regression. It can be described by the following relation.
\emph{Linear predictor = alpha + beta_predictors * predictors +
beta_(predictors_interactions) (if interactions between predictors are included) + error}.

However, except for linear regression, the linear predictor must be transformed to
obtain the outcome. This is because of the way generalised linear regression attempts to
create a linear relationship between the outcome and predictors.

If \emph{'Y'} is the outcome, the linear predictor is \emph{'logit Y'} for binary
outcomes; therefore, inverse logit transformation must be performed to convert the
linear predictor to obtain the probability of an outcome. For count outcomes, the
linear predictor is \emph{'log Y'}; therefore, exponential transformation is
required to obtain the predicted number of events.

For time-to-event outcomes, the linear predictor gives the hazard of an event at
various time points for a subject at the given covariate levels. The function
\code{\link[survival]{basehaz}} provides a more clinically meaningful cumulative
hazard of an event by time \emph{'t'}, denoted as \emph{'H(t)'}. The function
\code{\link[survival]{basehaz}} provides the cumulative hazard of the event at each time
point denoted as \emph{'t'} for each subject in the 'training' set. The cumulative
hazard of the event by time \emph{'t'} of a new subject can be calculated using
the relation mentioned in the description of \code{\link[survival]{basehaz}}
(please see the section on calculating H(t;x)). Using this relation, one can find
the closest time point of a 'new' subject to the time points in the output of
\code{\link[survival]{basehaz}}, the corresponding cumulative hazard for the 'first'
subject (or any other subject for whom the cumulative hazard at each time point
is available), and the differences in covariate values between the 'new' subject
and the 'first' subject to calculate the cumulative hazard of event by the follow-up
time of the new subject. The survival probability is exp(-H(t)) (Simon et al. 2024);
therefore, one can calculate the probability of event by time \emph{'t'} as
\emph{1 - exp(-H(t))}.

For continous outcomes, no transformation of linear predictor is required to obtain
the outcome.

\strong{Obtaining linear predictors}
In regression models, we can get the get the transformed values of the linear
predictors (lp) (i.e., inverse logit transformation for binary outcome and exponential
for count outcomes) based on the regression model directly. For example, using the
type = "response" in \code{\link[stats]{predict}} function gives this information
directly for all outcome types other than time-to-event outcomes, which are analysed
with \code{\link[survival]{coxph}}. For time-to-event outcomes, the type = "expected"
gives the cumulative hazard by \emph{'t'} after adjusting for the covariates
(\code{\link[survival]{predict.coxph}}), from which one can estimate the probability
of the event by the time \emph{'t'} using the relations described above.

\strong{Missing linear predictors}
When predicting using the regression models directly as described above, there
must be no missing data for the predictors included in the model. One possibility
is to not make a prediction at all. However, in real life some of these predictors
will be missing but a decision must be made. One possibility is multiple imputation.
However, some assumptions about the missing data can be difficult to verify (Heymans
et al, 2022). Another possibility is to exclude the missing predictor (whose value
is missing) from the regression equation. Although the coefficient values would have
been different without the predictor, it is impossible to develop and validate for
all scenarios of missing predictors. This function calculates the linear predictor
by excluding the predictors which contain missing data (for that subject) using the
regression model developed on subjects without missing data. If the coefficient values
in the model indicates \emph{NA} (which should alert people to overfitting the data
or levels with sparse data), the variable level itself is removed from calculating the
linear predictor. To a large extent, the method used assumes that external validation
will be performed before changing clinical practice and the application of this method
compared to other methods of handling missing data must be compared as part of
external validation.

\strong{Conversion of probabilities of event (linear predictors) for binary and
time-to-event outcomes to event versus no event}
There are multiple ways of converting probabilities of event (linear predictors)
for binary and time-to-event outcomes to event versus no event. For example, one
can consider that the probabilities of event are from binomial distribution for binary
outcomes. Alternatively, one can choose an 'optimal threshold' (on the training set)
using the \code{\link[pROC]{roc}} and \code{\link[pROC]{coords}} functions. There are two
types of threshold calculated by the \code{\link[pROC]{coords}} function:
\emph{'Youden'} and \emph{'closest to top left'}. For further information, please
see \code{\link[pROC]{coords}}. Occasionally, it may not be possible to obtain
the threshold using \code{\link[pROC]{roc}} and \code{\link[pROC]{coords}}. A function
that performs a rough estimation of the threshold based on prevalence is included in the
source code of this function (please see 'calculate_heuristic_threshold' function,
included as part of this function).

\strong{Intercept-slope adjustment}
In regression models, the intercept and slope can be adjusted (Van Calster et al.,
2019). The calibration intercept and slope are calculated according to the
supplement of Van Calster et al., 2019. The paper provides details only for
logistic regression, but the procedures are based on glm, i.e., they are applicable
in glm models. The relation is regression equation is
\emph{Y = calibration_intercept + beta * linear predictor}.

Note that the linear predictor in the equation is used as variable rather than as
an offset term as with calculation of calibration intercept only. The linear predictors
must be back-transformed to the original scale before their use in the calibration
regression equation.

Robust methods for calibration slope adjustment for time-to-events are still being
developed. Until such methods become widely available, this function uses similar
principles as that described for binary outcomes for time-to-event outcomes. These
should be considered experimental until further evaluation of the performance of
calibration adjustment in external samples. It should be noted however, that for
time-to-event outcomes, Cox regression does not have an intercept separately, as
the intercept is included in the baseline hazard (SAS support 2017). Therefore,
with regards to time-to-event outcomes, there is no change to the intercept, but
there is a change to the slope when calibration adjusted models are created.

\strong{Model with with only the mandatory predictors but based on the coefficients
of the entire model}
This is solely for research purposes. Potential use of such a model with only
the mandatory predictors, but based on the coefficients of the entire model will
be to find the added value of measurement of optional predictors, particularly
when there is a single mandatory predictor, for example, a treatment. It will be
practically impossible to develop all the possible models with missing optional
predictors. This model has the potential to provide predictions in this situation.
}
\value{
\item{actual_training}{Actual values in the training sample.}
\item{predicted_training}{Predicted values in the training sample.}
\item{predicted_training_calibration_adjusted}{Predicted values after calibration
adjustment.}
\item{predicted_training_adjusted_mandatory_predictors_only}{Predicted values of
a model with only the mandatory predictors, but based on the coefficients of the
entire model.}
\item{actual_only_validation}{Actual values in the 'out-of-sample' subjects, i.e.,
the subjects excluded from the model development in each simulation.}
\item{predicted_only_validation}{Predicted values in the 'out-of-sample' subjects}
\item{predicted_only_validation_calibration_adjusted}{Predicted values in the
out-of-sample subjects after calibration adjustment}
\item{predicted_only_validation_adjusted_mandatory_predictors_only}{Predicted values
in the out-of-sample subjects using a model with only the mandatory predictors,
but based on the coefficients of the entire model.}
\item{actual_all_subjects}{Actual values in all subjects with outcomes.}
\item{predicted_all_subjects}{Predicted values in all subjects with outcomes.}
\item{predicted_all_subjects_calibration_adjusted}{Predicted values in all subjects
with outcomes after calibration adjustment.}
\item{predicted_all_subjects_adjusted_mandatory_predictors_only}{Predicted values
in all subjects using a model with only the mandatory predictors, but based on the
coefficients of the entire model.}
\item{lp_training}{Linear predictors in the 'training' sample.}
\item{lp_only_validation}{Linear predictors in the 'out-of-sample' subjects.}
\item{lp_all_subjects}{Linear predictors in all subjects with outcomes.}
\item{lp_training_calibration_adjusted}{Linear predictors in the training sample
after calibration adjustment.}
\item{lp_only_validation_calibration_adjusted}{Linear predictors in the
'out-of-sample' subjects after calibration adjustment.}
\item{lp_all_subjects_calibration_adjusted}{Linear predictors in all subjects
with outcomes after calibration adjustment.}
\item{lp_training_adjusted_mandatory_predictors_only}{Linear predictors in the
training sample  using a model with only the mandatory predictors, but based on
the coefficients of the entire model.}
\item{lp_only_validation_adjusted_mandatory_predictors_only}{Linear predictors in the
'out-of-sample' subjects using a model with only the mandatory predictors, but
based on the coefficients of the entire model.}
\item{lp_all_subjects_adjusted_mandatory_predictors_only}{Linear predictors in
all subjects with outcomes using a model with only the mandatory predictors,
but based on the coefficients of the entire model.}
\item{time_training}{Follow-up time in training sample (applicable only for
time-to-event outcomes.)}
\item{time_only_validation}{Follow-up time in 'out-of-sample' subjects
(applicable only for time-to-event outcomes.)}
\item{time_all_subjects}{Follow-up time in all subjects with outcomes (applicable
only for time-to-event outcomes.)}
\item{regression_model}{The regression model}
\item{html_file}{Some output in html format, which will be used for final output.}
\item{outcome}{Whether calculations could be made.}
}
\references{
Collins GS, Dhiman P, Ma J, Schlussel MM, Archer L, Van Calster B, et al.
Evaluation of clinical prediction models (part 1): from development to external
validation. Bmj. 2024;384:e074819.

Heymans MW, Twisk JWR. Handling missing data in clinical research.
J Clin Epidemiol. 2022 Nov;151:185-188.

SAS Support. https://support.sas.com/kb/24/457.html (accessed on 16 January 2026).

Simon G, Aliferis C. Appendix A: Models for Time-to-Event Outcomes. In:
Simon GJ, Aliferis C, editors. Artificial Intelligence and Machine Learning in
Health Care and Medical Sciences: Best Practices and Pitfalls [Internet].
Cham (CH): Springer. https://www.ncbi.nlm.nih.gov/books/NBK610554/
(accessed on 13 December 2025). 2024.

Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW,
Bossuyt P, et al. Calibration: the Achilles heel of predictive analytics.
BMC Medicine. 2019;17(1):230.
}
\author{
Kurinchi Gurusamy
}
\seealso{
\code{\link{prepare_datasets}}
\code{\link[stats]{glm}}
\code{\link[stats]{predict}}
\code{\link[survival]{coxph}}
\code{\link[survival]{basehaz}}
\code{\link[survival]{predict.coxph}}
\code{\link[pROC]{roc}}
\code{\link[pROC]{coords}}
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
# change this to FALSE
results <- create_specific_input_parameters(
  generic_input_parameters = generic_input_parameters,
  analysis_details_path = analysis_details_path, verbose = TRUE)
specific_input_parameters <- results$specific_input_parameters
# Set a seed for reproducibility - Please see details above
set.seed(generic_input_parameters$seed)
prepared_datasets <- {prepare_datasets(
  df = generic_input_parameters$df,
  simulations = generic_input_parameters$simulations,
  outcome_name = generic_input_parameters$outcome_name,
  outcome_type = generic_input_parameters$outcome_type,
  outcome_time = generic_input_parameters$outcome_time,
  verbose = TRUE)}
# There is no usually no requirement to call this function directly. This is used
# by the perform_analysis function to create the actual and predicted values.
specific_input_parameters_each_analysis <- specific_input_parameters[[1]]
actual_predicted_results_apparent <- {calculate_actual_predicted(
      prepared_datasets = prepared_datasets,
      outcome_name = generic_input_parameters$outcome_name,
      outcome_type = generic_input_parameters$outcome_type,
      outcome_time = generic_input_parameters$outcome_time,
      outcome_count = generic_input_parameters$outcome_count,
      develop_model = specific_input_parameters_each_analysis$develop_model,
      predetermined_model_text =
      specific_input_parameters_each_analysis$predetermined_model_text,
      mandatory_predictors = specific_input_parameters_each_analysis$mandatory_predictors,
      optional_predictors = specific_input_parameters_each_analysis$optional_predictors,
      mandatory_interactions = specific_input_parameters_each_analysis$mandatory_interactions,
      optional_interactions = specific_input_parameters_each_analysis$optional_interactions,
      model_threshold_method = specific_input_parameters_each_analysis$model_threshold_method,
      scoring_system = specific_input_parameters_each_analysis$scoring_system,
      predetermined_threshold = specific_input_parameters_each_analysis$predetermined_threshold,
      higher_values_event = specific_input_parameters_each_analysis$higher_values_event,
      each_simulation = 1, bootstrap_sample = FALSE, verbose = TRUE
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
    predetermined_model_text =
      specific_input_parameters_each_analysis$predetermined_model_text,
    mandatory_predictors = specific_input_parameters_each_analysis$mandatory_predictors,
    optional_predictors = specific_input_parameters_each_analysis$optional_predictors,
    mandatory_interactions = specific_input_parameters_each_analysis$mandatory_interactions,
    optional_interactions = specific_input_parameters_each_analysis$optional_interactions,
    model_threshold_method = specific_input_parameters_each_analysis$model_threshold_method,
    scoring_system = specific_input_parameters_each_analysis$scoring_system,
    predetermined_threshold = specific_input_parameters_each_analysis$predetermined_threshold,
    higher_values_event = specific_input_parameters_each_analysis$higher_values_event,
    each_simulation = each_simulation, bootstrap_sample = TRUE, verbose = TRUE
  )
})
}
