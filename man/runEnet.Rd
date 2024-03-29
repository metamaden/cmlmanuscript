% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/machinelearning.R
\name{runEnet}
\alias{runEnet}
\title{Run elastic net regression}
\usage{
runEnet(seset, seed = 2019)
}
\arguments{
\item{seset}{A SummarizedExperiment object.}
}
\value{
A `resultslist` object containing the fitted model and evaluation info.
}
\description{
This function runs elastic net regression (with alpha 0.5) using a SummarizedExperiment object. Credit base code: Jenny Smith.
}
