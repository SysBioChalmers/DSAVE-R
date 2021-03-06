% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DSAVEGenerateSNODataset.R
\name{DSAVEGenerateSNODataset}
\alias{DSAVEGenerateSNODataset}
\title{DSAVEGenerateSNODataset}
\usage{
DSAVEGenerateSNODataset(
  templDs,
  numCells = NULL,
  noiseLevel = 0,
  templDSForProfile = NULL
)
}
\arguments{
\item{templDs}{The input dataset (cell population), numeric matrix (can be sparse)}

\item{numCells}{(optional) Can be used to specify the number of cells. Defaults to the number of cells in the input dataset.}

\item{noiseLevel}{(optional) The noise level to add; defaults to 0 (no noise)}

\item{templDSForProfile}{(optional) can be used if you want to generate data from a different cell type -
defaults to templDs - the genes in the datasets need to be synchronized}
}
\value{
a matrix containing the generated dataset
}
\description{
Generates a Sampling Noise Only (SNO) dataset
}
\details{
This function uses a template dataset to compute the probability p that
a molecule from a certain gene should be picked each time a new UMI is
found. This is calculated from the counts value, i.e. counts/sum of all
counts. Counts are then sampled from a multinomial distribution
drawing the same number of UMIs that are in the original set using the
probabilities for the genes from the mean of the dataset sent in.
Important that this is really UMI counts in this case, not TPM!
Random multiplicative noise will be added if desired; (0 == no noise)
noiseLevel should be 0 or greater. A standard random normal distributed
noise multiplied by noiseLevel will be multiplied to the probabilities.
}
\author{
Juan Inda, <inda@chalmers.se>
}
