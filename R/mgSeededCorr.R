#' This function takes transcript data and metabolite data and seeded correlation networks based on your metabolite or gene of interest
#' @param geneData Dataframe of transcript data with samplie IDs as rownames and transcript IDs as column headers.  Sample IDs should match metabolite data input.
#' @param metaboliteData Dataframe of metabolite data with samplie IDs as rownames and metabolite IDs as column headers.  Sample IDs should match transcript data input.
#' @param geneSeed Gene of interest to be used as seed of correlation network.  Be sure to use metaboliteSeed if you want to seed with a metabolite.  Defaults to NULL
#' @param metaboliteSeed Metabolite of interest to be used as seed of correlation network.  Be sure to use geneSeed if you want to seed with a gene.  Defaults to NULL
#' @param correlationMethod Indcate which correlation method to use.  Can be any method accepted by cor function.  See ?cor for details. Defaults to "pearson".
#' @param minCoef Number between 0 and 1 indicating the minimum correlation coeficient cutoff to be included in network. Defaults to 0.5.
#' @param maxNodes Number indicating the maximum number of nodes desired in your output.
#' @import dynamicTreeCut
#' @import qgraph
#' @keywords Correlation Network
#' @export
#' @return Plots seeded correlation networks to default environment.  Returns a list of dataframe with edgelist for the plot.
#' @examples
#' data(kDataRaw)
#' data(rDataRaw)
#' result<-mgSeededCorr(kDataRaw, rDataRaw, metaboliteSeed = "C00334")



mgSeededCorr<-function(geneData, metaboliteData, geneSeed=NULL, metaboliteSeed=NULL, correlationMethod="pearson", minCoef=0.5, maxNodes=20){
  #format data input
  metabolites<-metaboliteData[order(row.names(metaboliteData)),]
  genes<-geneData[order(row.names(geneData)),]
  
  #Find seed values
  if(!is.null(metaboliteSeed)){
    seed<-metabolites[,metaboliteSeed,drop=FALSE]
  }
  
  if(!is.null(geneSeed)){
    seed<-genes[,geneSeed,drop=FALSE]
  }
  
  
  colnames(metabolites)<-paste("M_", colnames(metabolites))
  colnames(genes)<-paste("G_", colnames(genes))
  
    
  #Find gene/metabolite correlations with seed
  metaboliteCorr<-cor(seed, metabolites, method=correlationMethod)
  metaboliteCorr<-metaboliteCorr[,order(-abs(metaboliteCorr[1,])), drop=FALSE]
  metaboliteCorr<-metaboliteCorr[,c(1:min(ncol(metaboliteCorr),maxNodes)), drop=FALSE]
  geneCorr<-cor(seed, genes, method=correlationMethod)
  geneCorr<-geneCorr[,order(-abs(geneCorr[1,])), drop=FALSE]
  geneCorr<-geneCorr[,c(1:min(ncol(geneCorr),maxNodes)), drop=FALSE]
  metFilt<-metabolites[,colnames(metaboliteCorr[,which(abs(metaboliteCorr[1,,drop=FALSE])>minCoef),drop=FALSE]),drop=FALSE]
  geneFilt<-genes[,colnames(geneCorr[,which(abs(geneCorr[1,,drop=FALSE])>minCoef),drop=FALSE]), drop=FALSE]
  
  #merge pertinant gene metabolite data and generate filtered correlation matrix
  filtMerge<-merge(metFilt, geneFilt, by.x="row.names", by.y="row.names")
  filtMerge<-data.frame(filtMerge[,-1], row.names=filtMerge[,1])
  filtCorr<-cor(filtMerge, method=correlationMethod)
  
  #Set gene/metabolite specific colors
  sampIndex<-row.names(filtCorr)
  colIndex<-replace(replace(sampIndex, grep("M_", sampIndex), "blue"), grep("G_", sampIndex), "yellow")
  colnames(filtCorr)<-gsub("M_.", "", colnames(filtCorr))
  colnames(filtCorr)<-gsub("G_.", "", colnames(filtCorr))
  
  #Create edgeList
  Graph<-qgraph(filtCorr, threshold=0.5, fade=TRUE, colFactor=3, esize=5, color = colIndex)
  plot(Graph)
  edgeList<-cbind(row.names(filtCorr)[z$Edgelist$from], row.names(filtCorr)[z$Edgelist$to], z$Edgelist$weight)
  colnames(edgeList)<-c("To", "From", "Correlation_Coef")
  return(edgeList)
  
}