% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functional_utils.R
\name{get_feat_tables}
\alias{get_feat_tables}
\title{Extract Feature Tables from Functional Profile}
\usage{
get_feat_tables(func_profile)
}
\arguments{
\item{func_profile}{A list containing functional profiles and metadata.}
}
\value{
A list containing only the feature tables (e.g., UniRef, EC, KO, pathways) without metadata.
}
\description{
This function extracts and returns the feature tables (excluding metadata) from a functional profile object.
}
\examples{
# Example usage:
feat_tables <- get_feat_tables(func_profile)

}
\keyword{internal}
