\name{guess_data_types}
\alias{guess_data_types}
\title{
Guess data types
}
\description{
This function removes any columns where there is no data and makes guesses on
the data type. However, this relies on the data not being coded already. If the
data has been coded, the metadata generated can be used as a template that can be
modified and provided as input for \code{\link{process_data}}.
}
\usage{
guess_data_types(data_file_path)
}
\arguments{
\item{data_file_path}{The path for the data file where guessing the data types of
columns is necessary.}
}
\value{
\item{outcome}{Whether the operation was successfully performed}
\item{message}{Any information, particularly when the operation fails.}
\item{data}{The data after removing the columns without any data.}
\item{metadata}{Automated metadata is created based on the data. However, this relies
on data not previously coded, for example, if the status is coded as 0 and 1 rather
than 'absent' and 'present', the variable will be recognised as a quantitative
variable rather than categorical variable.}
\item{any_type}{All fields with data.}
\item{quantitative}{Fields recognised as quantitative.}
\item{numerical}{Fields recognised as continuous.}
\item{count}{Fields recognised as count. Count data is recognised from the field
name. If a field name starts with 'Number of', it is considered as count data.}
\item{categorical}{Fields recognised as categorical data.}
\item{nominal}{Fields recognised as nominal data. All categorical data with more
than two levels are recognised as nominal data.}
\item{binary}{Fields recognised as binary data. All categorical data with only
two levels are recognised as binary data.}
\item{ordinal}{Fields recognised as ordinal data. Any categorical data with
more than two levels and with the second character of all the levels being an '_'
are recognised as ordinal data.}
\item{date}{Fields recognised as date.}
\item{time}{Fields recognised as time.}
}
\author{
Kurinchi Gurusamy
}
\seealso{
\code{\link{process_data}}
}
\examples{
\dontrun{
data_file_path <- paste0(tempdir(), "/df.csv")
write.csv(penguins, data_file_path, row.names = FALSE, na = "")
guessed_data_types <- guess_data_types(data_file_path = data_file_path)
}
}
