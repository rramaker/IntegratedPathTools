#' View metabolite foldchanges overlaid on KEGG pathways
#'
#' This function is a wrapper for the pathview package which allows you to plot data on KEGG pathways
#' @param geneResult Dataframe with transcript IDs as rownames, transcript signifigance values (between 0 and 1) as the first column, and optionally transcript fold change values as the second column. Defaults to NULL
#' @param proteinResult Dataframe with protein accession numbers as rownames, protein signifigance values (between 0 and 1) as the first column, and optionally protein fold change values as the second column. Defaults to NULL
#' @param metaboliteResult Dataframe with metabolite IDs as rownames, metabolite signifigance values (between 0 and 1) as the first column, and optionally metabolite fold change values as the second column. Defaults to NULL
#' @param cutoff The FDR to use for plotting. Defaults to 1
#' @param pathway.id The name of the KEGG pathway to plot (e.g. "hsa04727")
#' @param geneColors  A vector of three colors to use to represent gene/protein fold changes. Defaults to c("green", "grey","red").  
#' @param metaboliteColors  A vector of three colors to use to represent metabolite fold changes.  Defaults to c("blue","grey","yellow")
#' @keywords view kegg path
#' @import pathview
#' @export
#' @return Exports a png file to your working directory
#' @examples
#' data(kData)
#' data(rData)
#' library(pathview)
#' data(demo.paths)
#' viewPath(geneResult=kData, metaboliteResult=rData, pathway.id=demo.paths$sel.paths[1])
#' @note
#' viewPath is powered by the following open source databases.  Commercial use and/or redistribution may restricted.  Please see respective terms of use pages and citations for more details.
#'
#' #' KEGG
#' 
#' Terms of Use: http://www.kegg.jp/kegg/legal.html
#' 
#' Citations:
#' 
#' Kanehisa, M., Goto, S., Sato, Y., Kawashima, M., Furumichi, M., and Tanabe, M.  Data, information, knowledge and principle: back to metabolism in KEGG. Nucleic Acids Res. 42, D199-D205 (2014).
#' 
#' Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30 (2000).



viewPath = function(geneResult=NULL, proteinResult=NULL, metaboliteResult=NULL, pathway.id, geneColors = c("green", "grey","red"), metaboliteColors=c("blue","grey","yellow"), cutoff = 1){
  data(KEGG_Chebi, envir=environment())
  KEGG_Chebi<-KEGG_Chebi
  
  
  #error checking
  if(!is.null(proteinResult)&!is.null(geneResult)){
    stop("Can only provide geneResult OR proteinResult, not both")
  }
  #generate cpd.data for pathview function if provided
  if(!is.null(metaboliteResult)){
    KEGGfilt<-KEGG_Chebi[which(row.names(KEGG_Chebi)%in%row.names(metaboliteResult)),,drop=F]
    metaboliteResult<-merge(KEGGfilt, metaboliteResult, by.x="row.names", by.y="row.names")
    metaboliteResult<-data.frame(metaboliteResult[,-c(1:2)], row.names=metaboliteResult[,2])
    mToView <- metaboliteResult[order(metaboliteResult[,1]),]
    mToView <- as.matrix(subset(mToView, mToView[,1] <= cutoff))[,2]    
  }
  #generate gene.data for pathview function if provided
  if(!is.null(geneResult)){
    row.names(geneResult)<-gsub("X", "", row.names(geneResult))
    gToView <- geneResult[order(geneResult[,1]),]
    gene_proteinResult <- as.matrix(subset(gToView, gToView[,1] <= cutoff))[,2]
    #execute pathview function
    pathway <- pathview::pathview(gene.data = gene_proteinResult, cpd.data=mToView, pathway.id=pathway.id, cpd.idtype="kegg", gene.idtype="entrez", low=list(gene=geneColors[1], cpd = metaboliteColors[1]), mid = list(gene=geneColors[2], cpd=metaboliteColors[2]), high = list(gene=geneColors[3], cpd=metaboliteColors[3]))
  }
    
  if(!is.null(proteinResult)){
    pToView <- proteinResult[order(proteinResult[,1]),]
    gene_proteinResult <- as.matrix(subset(pToView, pToView[,1] <= cutoff))[,2]
    #execute pathview function
    pathway <- pathview::pathview(gene.data = gene_proteinResult, cpd.data=mToView, pathway.id=pathway.id, cpd.idtype="kegg", gene.idtype="ACCNUM", low=list(gene=geneColors[1], cpd = metaboliteColors[1]), mid = list(gene=geneColors[2], cpd=metaboliteColors[2]), high = list(gene=geneColors[3], cpd=metaboliteColors[3]))
  }
}