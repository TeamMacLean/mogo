% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{enricher_to_david}
\alias{enricher_to_david}
\title{conver enricher object to DAVID format}
\usage{
enricher_to_david(enrich, termsfile = "inst/extdata/MG8_mart_export.txt")
}
\arguments{
\item{enrich}{enricher object from do_enrich}

\item{termsfile}{TSV file from biomart export of GO terms}
}
\value{
DAVID format dataframe
}
\description{
conver enricher object to DAVID format
}
