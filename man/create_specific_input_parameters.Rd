\name{create_specific_input_parameters}
\alias{create_specific_input_parameters}
\title{
Create specific input parameters
}
\description{
This function converts the analysis details that you provide in a \emph{'csv'} file to a
list that will be used for analysis. Only valid rows are included in this list.
If there are invalid rows, the reasons for rejecting the rows are provided. Any revisions
that make the rows valid are performed when possible.
}
\usage{
create_specific_input_parameters(generic_input_parameters, analysis_details_path,
verbose)
}
\arguments{
\item{generic_input_parameters}{The generic input parameters that you generated
with \code{\link{create_generic_input_parameters}} (recommended) or manually.}
\item{analysis_details_path}{The path to the \emph{'csv'} file containing the following
columns.
\emph{'name'}: The name of the analysis. Provide a short name. This will be
displayed on graphs and will also be used for naming the analysis. The name must be
unique. If the name does not follow the naming convention for R objects,
a suitable name is created.
\emph{'analysis_title'}: The title of the analysis to be displayed. If the title
is missing, the \emph{'name'} is used as the \emph{'analysis_title'}.
\emph{'develop_model'}: \emph{TRUE}, if you want to develop a model; \emph{FALSE},
if you want to use a scoring system with a predetermined threshold (if applicable).
\emph{'predetermined_model_text'}: You can create the model text from the mandatory
and optional predictors and interactions or for finer control of the model, you
can provide the model text directly.
\emph{'mandatory_predictors'}: Predictors that must be included in the model. These
should be provided even if you provide the 'predetermined_model_text'.
\emph{'optional_predictors'}: Optional predictors that may be included in the model
by \code{\link{step}}. These should be provided even if you provide the
'predetermined_model_text'.
\emph{'mandatory_interactions'}: Interactions that must be included in the model. These
should be provided even if you provide the 'predetermined_model_text'.
\emph{'optional_interactions'}: Optional interactions that may be included in the
model by \code{\link{step}}. These should be provided even if you provide the
'predetermined_model_text'.
\emph{'model_threshold_method'}: One of \emph{'youden'}, \emph{'topleft'},
\emph{'heuristic'}. Please see description in \code{\link{calculate_actual_predicted}}.
\emph{'scoring_system'}: Name of the pre-existing scoring system. This is ignored
if \emph{develop_model} is \emph{TRUE}.
\emph{'predetermined_threshold'}: Pre-determined threshold of the pre-existing
scoring system. This is mandatory when \emph{develop_model} is \emph{FALSE} and when the
\emph{outcome_type} is \emph{'binary'} or \emph{'time-to-event'}. This is ignored if
\emph{develop_model} is \emph{TRUE} or when the \emph{outcome_type} is \emph{'quantitative'}.
\emph{'higher_values_event'}: \emph{TRUE} if higher values of the pre-existing system
indicates event and \emph{FALSE} otherwise. This is mandatory when \emph{develop_model} is
\emph{FALSE} and when the \emph{outcome_type} is \emph{'binary'} or \emph{'time-to-event'}.
This is ignored if \emph{develop_model} is \emph{TRUE} or when the \emph{outcome_type} is
\emph{'quantitative'}.
}
\item{verbose}{\emph{TRUE} if the outcome message must be displayed and \emph{FALSE}
otherwise.}
}
\value{
\item{outcome}{The outcome containing the processing details. If some corrections
were made, the corrections are included in the outcome. If some rows were not valid,
the reason for the row not being valid is provided.}
\item{specific_input_parameters}{A list with the information for further analyses.}
}
\author{
Kurinchi Gurusamy
}
\examples{
library(survival)
colon$status <- factor(as.character(colon$status))
generic_input_parameters <- create_generic_input_parameters(
  general_title = "Prediction of colon cancer death", simulations = 2000,
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
}
