\name{get_outcome_status_at_specific_time}
\alias{get_outcome_status_at_specific_time}
\title{
Get outcome status at specific time
}
\description{
Some may want to calculate the outcome status at a specified time to predict the
outcome at a specific time. This function gets the outcome status at the specific
time.
}
\usage{
get_outcome_status_at_specific_time(df, status_field, time_field, specific_time)
}
\arguments{
\item{df}{The dataframe which contains the data.}
\item{status_field}{The field in the dataframe that indicates the status, i.e.,
event or no event.}
\item{time_field}{The field in the dataframe that indicates the follow-up time.}
\item{specific_time}{The time point at which the outcome status should be
calculated.}
}
\value{
\item{outcome}{Whether the operation was successfully performed}
\item{message}{Any information, particularly when the operation fails.}
\item{new_data}{The data with the existing status and time fields replaced by the
updated status and time. The original status and time are available with a prefix
"unmodified_".}
}
\author{
Kurinchi Gurusamy
}
\examples{
library(survival)
# Use the dataset colon as example
# Replace existing status with the status at 365 days
results <- get_outcome_status_at_specific_time(df = colon,
status_field = "status", time_field = "time", specific_time = 365)
results$outcome
results$message
# Display first 10 rows to show how the status and time have been modified
results$new_data[1:10, c("status", "time", "unmodified_status", "unmodified_time")]
}

