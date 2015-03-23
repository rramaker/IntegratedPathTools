#' This function takes transcript data and metabolite data and generates unsupervised gene-metabolite networks
#' @param geneData Dataframe of transcript data with samplie IDs as rownames and transcript IDs as column headers.  Sample IDs should match metabolite data input.
#' @param metaboliteData Dataframe of metabolite data with samplie IDs as rownames and metabolite IDs as column headers.  Sample IDs should match transcript data input.
#' @param moduleDriver Input "metabolite" or "gene" to dictate whether networks should be driven by gene or metabolite data.  This script is optimized for using metabolite data as the moduleDrive.  If using "gene" is much more computationally intensive.  Be sure to use a heavily filtered gene list of less than 10000 genes for reasonable computation time.
#' @param correlationMethod Indcate which correlation method to use.  Can be any method accepted by cor function.  See ?cor for details. Defaults to "pearson".
#' @param minCoef Number between 0 and 1 indicating the minimum correlation coeficient cutoff to be included in network. Defaults to 0.5.
#' @param minModuleSize Number indicating the minimum size of driver modules used for network construction. See minClusterSize in ?cutreeDynamic. Defaults to 10
#' @param maxPlotEdges Number indicating the maximum number of edges to draw on correlation plot.  Priority is given to highest correlation.  Defaults to 500.
#' @param numNetworks Number indicating the number of networks to construct.  Networks will be constructed in order module size.
#' @param plotHeatmap Logical indicating whether to plot metabolite-gene correlation heatmap.  This will greatly increase computation time.  Defaults to FALSE
#' @import gplots
#' @import dynamicTreeCut
#' @import qgraph
#' @import RColorBrewer
#' @keywords Correlation Network
#' @export
#' @return Plots correlation networks to default environment.  Returns a list of dataframes.  The first returned dataframe gives the connectivity of each driving element of each module.  The following dataframes contain the full edgelist information for modules in order of module size.
#' @examples
#' data(kDataRaw)
#' data(rDataRaw)
#' result<-mgCorr(kDataRaw, rDataRaw, numNetworks=1)

mgCorr<-function(geneData, metaboliteData, moduleDriver = c("metabolites", "genes"), correlationMethod = "pearson", minCoef = 0.5, minModuleSize=10, maxPlotNodes=500, numNetworks = 2, plotHeatmap=FALSE ){
  #format data inputted
  metabolites<-metaboliteData[order(row.names(metaboliteData)),]
  genes<-geneData[order(row.names(geneData)),]
  colnames(metabolites)<-paste("M_", colnames(metabolites))
  colnames(genes)<-paste("G_", colnames(genes))
  
  #Generate correlation matrix according to moduleDriver
  if(moduleDriver[1] == "metabolites"){
    fullCorMatrix<-cor(metabolites, genes, method=correlationMethod)
    corDist<-dist(fullCorMatrix)
    clust<-hclust(corDist)
  }
  
  if(moduleDriver[1]== "genes"){
    fullCorMatrix<-cor(genes, metabolites, method=correlationMethod)
    corDist<-dist(fullCorMatrix)
    clust<-hclust(corDist)
  }
  
  #Identify cluster modules with cutreeDynamic
  moduleResult<-cutreeDynamic(clust, minClusterSize = minModuleSize, distM=as.matrix(corDist), verbose=0)
  resultFrame<-data.frame(moduleResult, row.names=clust$labels)
  
  #Construct heatmap of correlation network cluster
  if(plotHeatmap==TRUE){
    print("constructing heatmap, to reduce run time set plotHeatmap to FALSE")
    toMap<-merge(resultFrame, fullCorMatrix, by.x="row.names", by.y="row.names")
    toMap<-data.frame(toMap[,-1], row.names=toMap[,1])
    toMapMatrix<-as.matrix(t(toMap[,-1]))
    colors<-rainbow(max(toMap[,1]+1))
    colSideColors<-c(colors[(toMap[,1]+1)])
    heatColors=c("red", "black", "green")
    heatmap.2( x = toMapMatrix, hclustfun = function(d) hclust(d), trace = "none", breaks = quantile(toMapMatrix, probs = seq(0,1, 0.05), na.rm=T), col = colorRampPalette(colors = heatColors), ColSideColors = colSideColors, margins = c(4,10), offsetRow = 0, offsetCol = 0, key = T, keysize = 0.5, labRow=F, labCol=F, dendrogram='col')
  }
  
  #empty list to return final result
  networkBind<-vector("list", numNetworks)
  #empty dataframe to input connectivity scores
  networkScoreFrame<-data.frame()
  
  for(n in 1:numNetworks){
    print(paste("building network ",  n), sep="")
    
    #Find genes correlated with metabolites in a module
    currentModule<-row.names(resultFrame[which(resultFrame[,1]==n),,drop=F])
    geneFind<-fullCorMatrix[currentModule,]
    geneHits<-apply(geneFind, 2, function(x) max(abs(x)))  
    moduleGenes<-colnames(geneFind[,which(geneHits>minCoef)])
    if(length(moduleGenes)==0){
      stop("minCoef too stringent to draw the number of networks desired")
    }
    
    #Subset correlation matrix to pertinant module metabolites and genes
    zFilt<-fullCorMatrix[currentModule, moduleGenes]
    
    #Generate edge list
    edgeInfo<-data.frame()
    for(i in 1:ncol(zFilt)){
      x<-row.names(zFilt[which(abs(zFilt[,i])>minCoef),,drop=F])
      y<-data.frame()
      for(t in 1:length(x)){
        y[t,1]<-colnames(zFilt[,i, drop=F])
        y[t,2]<-x[t]
        y[t,3]<-zFilt[x[t],i]
      }
      edgeInfo<-rbind(edgeInfo, y)
    }
    
    #Order and name edgelist
    finalEdge<-edgeInfo[order(-abs(edgeInfo[,3])),]
    colnames(finalEdge)<-c(paste("Module_", n,"_Gene_ID", sep=""), paste("Module_", n,"_Metabolite_ID", sep=""), paste("Module_", n,"_Corr_Coef", sep=""))
    
    #Calculate connectivity scores for module
    currentScoreFrame<-data.frame()
    for(q in 1:length(currentModule)){
      currentScoreFrame[q,1]<- currentModule[q]
      currentScoreFrame[q,2]<- nrow(finalEdge[which(finalEdge[,2]==currentModule[q]),])
      currentScoreFrame[q,3]<- n    
    }
    currentScoreFrame<-currentScoreFrame[order(-currentScoreFrame[,2]),]
    networkScoreFrame<-rbind(networkScoreFrame, currentScoreFrame)
    
    #Trim edge list according to plotting parameters
    currentPlot<-finalEdge[1:min(c(nrow(finalEdge),maxPlotNodes)),]
    trimMet<-unique(currentPlot[,2])
    for(s in 1:length(trimMet)){
      trimGene<-currentPlot[which(currentPlot[,2]==trimMet[s]),1]
      if(sum(currentPlot[,1]%in%trimGene)==length(trimGene)){
        currentPlot<-currentPlot[-which(currentPlot[,2]==trimMet[s]),]
      }
    }
    
    #Color code metabolite vs. gene nodes
    sampIndex<-unique(c(currentPlot[,1], currentPlot[,2]))
    colIndex<-replace(replace(sampIndex, grep("M_", sampIndex), "blue"), grep("G_", sampIndex), "yellow")
    
    #Remove metabolite/gene markers
    finalEdge[,1]<-gsub("G_", "", finalEdge[,1])
    finalEdge[,2]<-gsub("M_", "", finalEdge[,2])
    currentPlot[,1]<-gsub("G_", "", currentPlot[,1])
    currentPlot[,2]<-gsub("M_", "", currentPlot[,2])
    
    #plot network
    qgraph(currentPlot, color=colIndex, edge.width=abs(currentPlot[,3]), arrows=FALSE)
    #Add edge list to returned result
    networkBind[[n+1]]<-finalEdge
  }
  colnames(networkScoreFrame)<-c("Metabolite", "Connectivity", "Module_Number")
  #Add connectivity values to returned result
  networkBind[[1]]<-networkScoreFrame
  
  return(networkBind)
}
