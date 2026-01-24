\name{create_generic_input_parameters}
\alias{create_generic_input_parameters}
\title{
Create generic input parameters
}
\description{
This function simply checks whether the input parameters are correct and if correct,
creates a list from the input parameters. This also makes some corrections when
possible (i.e., when there were minor correctable issues in the input parameters).
}
\usage{
create_generic_input_parameters(general_title, simulations, simulations_per_file,
seed, df, outcome_name, outcome_type, outcome_time, outcome_count, verbose)
}
\arguments{
\item{general_title}{A general title for your analysis}
\item{simulations}{The number of simulations required. Usually at least 300 to 500
simulations are a minimum. Increasing the simulations leads to more reliable results.
The default value of 2000 simulations should provide reasonably reliable results.}
\item{simulations_per_file}{This is to manage the memory requirements. The default
value of 20 simulations per file should work in most instances.}
\item{seed}{Please see \code{\link{prepare_datasets}} for details.}
\item{df}{The dataset used for the analysis. This must be provided as a
dataframe. Data in files can be converted to dataframes with appropriate field types
using \code{\link{process_data}}.}
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
\item{verbose}{\emph{TRUE} if the outcome message must be displayed and \emph{FALSE}
otherwise.}
}
\value{
\item{outcome}{The outcome containing the processing details. If some corrections
were made, the corrections are included in the outcome. If there was a fatal error,
the reason for the fatal error is provided.}
\item{generic_input_parameters}{A list with information for further analyses. If there
was a fatal error, the reason for the fatal error is displayed and
generic_input_parameters is \emph{NULL}.}
}
\author{
Kurinchi Gurusamy
}
\seealso{
\code{\link{process_data}}
\code{\link{prepare_datasets}}
}
\examples{
# Correct parameters ####
# Binary outcome
# verbose is TRUE, therefore, the outcome message will be displayed
results <- create_generic_input_parameters(
general_title = "Prediction of penguin species", simulations = 2000,
simulations_per_file = 20, seed = 1, df = penguins, outcome_name = "species",
outcome_type = "binary", outcome_time = NA, outcome_count = FALSE, verbose = TRUE)
generic_input_parameters <- results$generic_input_parameters
generic_input_parameters

# Time-to-event outcome
library(survival)
# The field 'status' is provided as numeric. This must be converted to factor. In
# this example, we can convert this to factor using a command. For conversion of more
# columns, please use process_data function.
colon$status <- factor(as.character(colon$status))
# verbose is FALSE, therefore, the outcome message will not be displayed, but the
# outcome is stored.
results <- create_generic_input_parameters(
general_title = "Prediction of colon cancer death", simulations = 2000,
simulations_per_file = 20, seed = 1, df = colon, outcome_name = "status",
outcome_type = "time-to-event", outcome_time = "time", outcome_count = FALSE,
verbose = FALSE)
# Display outcome
results$outcome
# Display generic_input_parameters
generic_input_parameters <- results$generic_input_parameters
generic_input_parameters

# Continuous outcome
# verbose is not supplied, therefore, the outcome message will be displayed as
# this is the default.
results <- create_generic_input_parameters(
general_title = "Prediction of iris petal length", simulations = 2000,
simulations_per_file = 20, seed = 1, df = iris, outcome_name = "Petal.Length",
outcome_type = "quantitative", outcome_time = NA, outcome_count = FALSE)
generic_input_parameters <- results$generic_input_parameters
generic_input_parameters

# Count outcomes
results <- create_generic_input_parameters(
general_title = "Prediction of warp breaks", simulations = 2000,
simulations_per_file = 20, seed = 1, df = warpbreaks, outcome_name = "breaks",
outcome_type = "quantitative", outcome_time = NA, outcome_count = TRUE)
generic_input_parameters <- results$generic_input_parameters
generic_input_parameters

# Non fatal errors ####
results <- create_generic_input_parameters(
general_title = "", simulations = "Use default",
simulations_per_file = "Use default", seed = "Use default",
df = warpbreaks, outcome_name = "breaks",
outcome_type = "quantitative", outcome_time = "Use default", outcome_count = TRUE,
verbose = TRUE)
generic_input_parameters <- results$generic_input_parameters
generic_input_parameters

# Fatal error ####
# Note the dataframe name supplied within quotes.
results <- create_generic_input_parameters(
general_title = "", simulations = "Use default",
simulations_per_file = "Use default", seed = "Use default",
df = "warpbreaks", outcome_name = "breaks", outcome_type = "quantitative",
outcome_time = "Use default", outcome_count = TRUE, verbose = TRUE)
generic_input_parameters <- results$generic_input_parameters
generic_input_parameters
}
