% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R:CBHM.R
\name{decidePar}
\alias{decidePar}
\title{Function to decide tuning parameter a and b}
\usage{
decidePar(cohortsize, ntype, ntrial, p0, p1, var.small, var.big)
}
\arguments{
\item{cohortsize}{vector format, sample size for each cohort per disease type}

\item{ntype}{number of disease types/subgroups}

\item{ntrial}{the number of simulated trials}

\item{p0}{the value of null hypothesis rate}

\item{p1}{the value of alternative hypothesis rate}

\item{var.small}{the prespecified small value of shrinkage parameter sigma^2 for calibrating a and b}

\item{var.big}{the prespecified big value of shrinkage parameter sigma^2 for calibrating a and b}
}
\value{
tuning parameters a and b
}
\description{
Function to decide tuning parameter a and b
}
