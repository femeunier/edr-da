% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/specfile.R
\name{filter_specfile}
\alias{filter_specfile}
\title{Extract spectral data from HDF5 file based on criteria}
\usage{
filter_specfile(aviris_h5, ...)
}
\arguments{
\item{aviris_h5}{HDF5 file object}

\item{...}{Filters, similar to those used in `dplyr::filter`}
}
\value{
`spectra` object (matrix) containing AVIRIS data
}
\description{
Extract spectral data from HDF5 file based on criteria
}
