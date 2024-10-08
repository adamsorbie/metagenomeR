% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomic_utils.R
\name{return_level}
\alias{return_level}
\title{Map Taxonomic Level to Abbreviation}
\usage{
return_level(level)
}
\arguments{
\item{level}{A character string representing the taxonomic level (e.g., "Species", "Genus", "Family").}
}
\value{
A character string representing the abbreviation of the specified taxonomic level.
}
\description{
This function maps a specified taxonomic level to its corresponding abbreviation.
}
\examples{
# Example usage:
level_abbreviation <- return_level("Genus")

}
