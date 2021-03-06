% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bulkTotalVar1vs1.R
\docType{data}
\name{bulkTotalVar1vs1}
\alias{bulkTotalVar1vs1}
\title{bulkTotalVar1vs1}
\format{
list of lists, where the inner list contains the expression interval and the variation.
}
\source{
We downloaded 8 bulk RNA-Seq samples (FPKM) from the BLUEPRINT Epigenome Project31.
The samples were taken from the project EGAD00001001173 and have the following sample IDs:
S002EV11, S004M711, S007DD11, S007G711, S008H111, S009W411, S0018A13, and S0041C11.
The FPKM data were normalized using TMM32 before use, and scaled to an average count of 10^6 per sample.
}
\usage{
bulkTotalVar1vs1
}
\description{
Average total variation between two bulk samples at different expression intervals (reasonably comparable
to TPM).
}
\keyword{datasets}
