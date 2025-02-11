% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/mgCorr.R
\name{mgCorr}
\alias{mgCorr}
\title{This function takes transcript data and metabolite data and generates unsupervised gene-metabolite networks}
\usage{
mgCorr(geneData, metaboliteData, moduleDriver = c("metabolites", "genes"),
  correlationMethod = "pearson", minCoef = 0.5, minModuleSize = 10,
  maxPlotNodes = 500, numNetworks = 2, plotHeatmap = FALSE)
}
\arguments{
\item{geneData}{Dataframe of transcript data with samplie IDs as rownames and transcript IDs as column headers.  Sample IDs should match metabolite data input.}

\item{metaboliteData}{Dataframe of metabolite data with samplie IDs as rownames and metabolite IDs as column headers.  Sample IDs should match transcript data input.}

\item{moduleDriver}{Input "metabolite" or "gene" to dictate whether networks should be driven by gene or metabolite data.  This script is optimized for using metabolite data as the moduleDrive.  If using "gene" is much more computationally intensive.  Be sure to use a heavily filtered gene list of less than 10000 genes for reasonable computation time.}

\item{correlationMethod}{Indcate which correlation method to use.  Can be any method accepted by cor function.  See ?cor for details. Defaults to "pearson".}

\item{minCoef}{Number between 0 and 1 indicating the minimum correlation coeficient cutoff to be included in network. Defaults to 0.5.}

\item{minModuleSize}{Number indicating the minimum size of driver modules used for network construction. See minClusterSize in ?cutreeDynamic. Defaults to 10}

\item{numNetworks}{Number indicating the number of networks to construct.  Networks will be constructed in order module size.}

\item{plotHeatmap}{Logical indicating whether to plot metabolite-gene correlation heatmap.  This will greatly increase computation time.  Defaults to FALSE}

\item{maxPlotEdges}{Number indicating the maximum number of edges to draw on correlation plot.  Priority is given to highest correlation.  Defaults to 500.}
}
\value{
Plots correlation networks to default environment.  Returns a list of dataframes.  The first returned dataframe gives the connectivity of each driving element of each module.  The following dataframes contain the full edgelist information for modules in order of module size.
}
\description{
This function takes transcript data and metabolite data and generates unsupervised gene-metabolite networks
}
\examples{
data(kDataRaw)
data(rDataRaw)
result<-mgCorr(kDataRaw, rDataRaw, numNetworks=1)
}
\keyword{Correlation}
\keyword{Network}

