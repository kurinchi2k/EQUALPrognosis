\name{process_data}
\alias{process_data}
\title{
Process the data
}
\description{
This takes a dataset and the metadata for the dataset and creates R data frames
in a format required for the subsequent steps.
}
\usage{
process_data(data_file_path, metadata_file_path)
}
\arguments{
\item{data_file_path}{Path to the dataset}
\item{metadata_file_path}{Path to the metadata file}
}
\details{
The metadata should contain the following information as a minimum.
\emph{variable}: this is the name of the variable and should match the column names
of the dataset.
\emph{data_type}: \emph{'numerical'} for continuous variables, \emph{'count'} for count
variables, \emph{'binary} for binary categorical variables, \emph{'nominal'} for unordered
categorical variables with more than 2 levels, \emph{'ordinal'} for ordered categorical
variables, \emph{'date'} for variables stored as date, and \emph{'time'} for variables
stored containing the time of the day.

Optional information includes the following.
\emph{reference}: Reference category for \emph{binary} and \emph{nominal} variables. This
should be a category existing in the variable.
\emph{ordinal levels}: the levels of ordinal data from lower to higher order, separated
by ";". This must include all the levels in the data.

You can use \code{\link{guess_data_types}} as a starting point for the metadata, which is
included in the output list of the \code{\link{guess_data_types}} function.
}
\value{
\item{outcome}{Whether the operation was successfully performed}
\item{message}{Any information, particularly when the operation fails.}
\item{data_processed}{The data which has been modifed according to the metadata when
correct parameters are provided}.
\item{any_type}{All fields.}
\item{quantitative}{Fields recognised as quantitative.}
\item{numerical}{Fields recognised as continuous.}
\item{count}{Fields recognised as count.}
\item{categorical}{Fields recognised as categorical data.}
\item{nominal}{Fields recognised as nominal data}
\item{binary}{Fields recognised as binary data.}
\item{ordinal}{Fields recognised as ordinal data.}
\item{date}{Fields recognised as date.}
\item{time}{Fields recognised as time.}
}
\author{
Kurinchi Gurusamy
}
\seealso{
\code{\link{guess_data_types}}
}
\examples{
library(survival)
# Use the dataset colon as example
# Select only the survival for these examples (etype == 2)
data_file_path <- paste0(tempdir(), "/df.csv")
write.csv(colon[colon$etype == 2, ], data_file_path, row.names = FALSE, na = "")
metadata <- {data.frame(
  variable = c("id","study","rx","sex","age",
               "obstruct","perfor","adhere","nodes","status",
               "differ","extent","surg","node4","time",
               "etype"),
  data_type = c("nominal", "nominal", "nominal", "binary", "numerical",
                "binary", "binary", "binary", "count", "binary",
                "ordinal", "ordinal", "binary", "binary", "numerical",
                "nominal"),
  reference = c(NA, NA, "Obs", 0, NA,
                0, 0, 0, NA, 0,
                NA, NA, 0, 0, NA,
                NA),
  ordinal_levels = c(NA, NA, NA, NA, NA,
                     NA, NA, NA, NA, NA,
                     "1;2;3", "1;2;3;4", NA, NA, NA,
                     NA),
  comments = NA
)}
metadata_file_path <- paste0(tempdir(), "/metadata.csv")
write.csv(metadata, metadata_file_path, row.names = FALSE, na = "")
processed_data <- process_data(data_file_path, metadata_file_path)
}
