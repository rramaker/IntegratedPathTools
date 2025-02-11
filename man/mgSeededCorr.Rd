% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/mgSeededCorr.R
\name{mgSeededCorr}
\alias{mgSeededCorr}
\title{This function takes transcript data and metabolite data and seeded correlation networks based on your metabolite or gene of interest}
\usage{
mgSeededCorr(geneData, metaboliteData, geneSeed = NULL,
  metaboliteSeed = NULL, correlationMethod = "pearson", minCoef = 0.5,
  maxNodes = 20)
}
\arguments{
\item{geneData}{Dataframe of transcript data with samplie IDs as rownames and transcript IDs as column headers.  Sample IDs should match metabolite data input.}

\item{metaboliteData}{Dataframe of metabolite data with samplie IDs as rownames and metabolite IDs as column headers.  Sample IDs should match transcript data input.}

\item{geneSeed}{Gene of interest to be used as seed of correlation network.  Be sure to use metaboliteSeed if you want to seed with a metabolite.  Defaults to NULL}

\item{metaboliteSeed}{Metabolite of interest to be used as seed of correlation network.  Be sure to use geneSeed if you want to seed with a gene.  Defaults to NULL}

\item{correlationMethod}{Indcate which correlation method to use.  Can be any method accepted by cor function.  See ?cor for details. Defaults to "pearson".}

\item{minCoef}{Number between 0 and 1 indicating the minimum correlation coeficient cutoff to be included in network. Defaults to 0.5.}

\item{maxNodes}{Number indicating the maximum number of nodes desired in your output.}
}
\value{
Plots seeded correlation networks to default environment.  Returns a list of dataframe with edgelist for the plot.
}
\description{
This function takes transcript data and metabolite data and seeded correlation networks based on your metabolite or gene of interest
}
\examples{
data(kDataRaw)
data(rDataRaw)
result<-mgSeededCorr(kDataRaw, rDataRaw, metaboliteSeed = "C00334")
}
\keyword{Correlation}
\keyword{Network}

