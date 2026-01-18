\name{prepare_datasets}
\alias{prepare_datasets}
\title{
Prepare simulated datasets from the entire dataset
}
\description{
This takes the dataset, prepared using the \code{\link{prepare_datasets}}
function and takes subsets of data for each simulation. For details, please see
below.
}
\usage{
prepare_datasets(df, simulations, seed, outcome_name, outcome_type,
outcome_time)
}
\arguments{
\item{df}{The dataset used for the analysis. This must be provided as a
dataframe. Data in files can be converted to dataframes with appropriate field types
using \code{\link{process_data}}.}
\item{simulations}{The number of simulations required. Usually at least 300 to 500
simulations are a minimum. Increasing the simulations leads to more reliable results.
The default value of 2000 simulations should provide reasonably reliable results.}
\item{seed}{Please see below for details.}
\item{outcome_name}{Name of the colummn that contains the outcome data. This
must be a column name in the \emph{'df'} provided as input.}
\item{outcome_type}{One of \emph{'binary'}, \emph{'time-to-event'}, \emph{'quantitative'}.
Count outcomes are included in \emph{'quantitative'} outcome type and can be differentiated
from continuous outcomes by specifying outcome_count as \emph{TRUE}. Please see examples
below.}
\item{outcome_time}{The name of the column that provides the follow-up time. This
is applicable only for \emph{'time-to-event'} outcome. For other outcome types,
enter \emph{NA}.}
}
\details{
\strong{Overview}
The input parameters are part of the generic_input_parameters created with
\code{\link{create_generic_input_parameters}}. In the first step, it excludes all
rows where the outcome is not available. For 'time-to-event' outcomes, the rows
without \emph{outcome_time} are also excluded. This forms the basis for the 'all_subjects'
dataset.

In the next step, subjects are sampled from the 'all_subjects' dataset. The sampling is
done using random methods. The starting point used in the random number generator is called a
'seed'. This determines all the subsequent numbers generated. The size of the sample is
the same as the original data set. This is done by sampling with replacement. The
subjects included in this sample are included for model development and the sample is called
'training' dataset.

Some subjects are not included in the 'training' dataset. These 'out-of-sample' subjects are
used only for validation. The dataset that includes only 'out-of-sample' subjects is called
'only_validation' dataset. It must be noted that some subjects in the 'all_subjects' will be
included more than once in the 'training' dataset because of the nature of the sampling.

While sampling, if all the subjects have the same outcome, this is not suitable for
our evaluation. Therefore, such simulations cannot be included in the analysis and
therefore, excluded. Such excluded seeds are called 'seeds_resulting_in_error'.
The 'seeds' that did not cause this error were called 'working_seeds'. These are provided
for reproducibility of results.

In addition to sampling, some additional processing is performed. In the 'out-of-sample'
evaluation, when there are ordinal factors absent in the 'training' dataset but present
in the 'out-of-sample' dataset, it results in errors. To avoid this, some levels of ordinal
factors are combined.

By default, seed 1 is used for initiating the random sequence. Choosing a different seed
might give slightly different results compared to seed 1.
}
\value{
\item{working_seeds}{Please see above}
\item{seeds_resulting_in_error}{Please see above}
\item{df_training_list}{A list of 'training' datasets, one for each simulation.}
\item{df_only_validation_list}{A list of datasets containing 'out-of-sample' subjects,
one for each simulation.}
\item{df_all_subjects_list}{This is the same for all simulations.}
}
\author{
Kurinchi Gurusamy
}
\seealso{
\code{\link{Random}}
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
  'Multiple mandatory and optional predictors', 'Multiple optional predictors only'),
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
  outcome_time = generic_input_parameters$outcome_time)}

}
}
