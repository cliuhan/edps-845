% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equip.R
\name{equip}
\alias{equip}
\title{Equipercentile Equating Analytic Standard Error}
\usage{
equip(x, y, ly = min(scales(y)), ky = max(scales(y)))
}
\arguments{
\item{x}{and y are frequency tables.}
}
\description{
This function is used to calcuate the analytic standard error
of equipercentile equating without smoothing under single- or random-
groups design.
}
\details{
The delta method is a commonly used to compute standard errors.
}
\examples{

# Get the analytic standard errors. 

library(equate)
rx <- as.freqtab(ACTmath[, 1:2])
ry <- as.freqtab(ACTmath[, c(1, 3)])
mod <- equip(x, y)
se <- mod$se

}

