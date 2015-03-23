#' This function takes dataframes with gene or protein info and finds drug targets
#' @param geneResult Dataframe with transcript IDs as rownames, transcript signifigance values (between 0 and 1) as the first column, and optionally transcript fold change values as the second column. Defaults to NULL
#' @param proteinResult Dataframe with Uniprot Accession IDs as rownames, Protein signifigance values (between 0 and 1) as the first column, and optionally protein fold change values as the second column. Defaults to NULL
#' @param cutoff Significance cutoff to use for selecting which genes/proteins to use for drug target analysis.
#' @export
#' @return Returns a dataframe with gene or protein targets listed for each drug capable of targeting significant genes or proteins.
#' @examples
#'
#' data(kData)
#' findDrugTargets(geneResult = kData, cutoff=0.1)
#' 
#' @note
#' findDrugTargets is powered by the following open source databases.  Commercial use and/or redistribution may restricted.  Please see respective terms of use pages and citations for more details.
#' 
#' DrugBank:
#' 
#' Terms of Use: http://www.drugbank.ca/about
#' 
#' Citations:
#' 
#' DrugBank 4.0: shedding new light on drug metabolism. Law V, Knox C, Djoumbou Y, Jewison T, Guo AC, Liu Y, Maciejewski A, Arndt D, Wilson M, Neveu V, Tang A, Gabriel G, Ly C, Adamjee S, Dame ZT, Han B, Zhou Y, Wishart DS.Nucleic Acids Res. 2014 Jan 1;42(1):D1091-7. PubMed ID: 24203711
#' 
#' DrugBank 3.0: a comprehensive resource for 'omics' research on drugs. Knox C, Law V, Jewison T, Liu P, Ly S, Frolkis A, Pon A, Banco K, Mak C, Neveu V, Djoumbou Y, Eisner R, Guo AC, Wishart DS.Nucleic Acids Res. 2011 Jan;39(Database issue):D1035-41.  PubMed ID: 21059682
#' 
#' DrugBank: a knowledgebase for drugs, drug actions and drug targets. Wishart DS, Knox C, Guo AC, Cheng D, Shrivastava S, Tzur D, Gautam B, Hassanali M.Nucleic Acids Res. 2008 Jan;36(Database issue):D901-6. PubMed ID: 18048412
#' 
#' DrugBank: a comprehensive resource for in silico drug discovery and exploration. Wishart DS, Knox C, Guo AC, Shrivastava S, Hassanali M, Stothard P, Chang Z, Woolsey J.Nucleic Acids Res. 2006 Jan 1;34(Database issue):D668-72.  PubMed ID: 16381955



findDrugTargets<-function(geneResult=NULL, proteinResult=NULL, cutoff){
  data(DrugBank_GeneTarget, envir=environment())
  data(DrugBank_ProteinTarget, envir=environment())
  DrugBank_GeneTarget<-DrugBank_GeneTarget
  DrugBank_ProteinTarget<-DrugBank_ProteinTarget
  
  
  if(!is.null(geneResult)){
    drugData<-DrugBank_GeneTarget
    testResult<-geneResult
  }
  if(!is.null(proteinResult)){
    drugData<-DrugBank_ProteinTarget
    testResult<-proteinResult
  }
  
  toORA <- testResult[order(testResult[,1]),]
  #Filter for only markers below the specified FDR cutoff
  sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
  sigFrame<- toORA[sigList,]
  drugFilt<-drugData[row.names(drugData)%in%sigList,]
  drugFilt<-drugFilt[,which(colSums(drugFilt, na.rm=TRUE)>0)]
  result<-data.frame()
  for(i in 1:ncol(drugFilt)){
    currentTargets<-row.names(drugFilt[which(drugFilt[,i]==1),, drop=FALSE])
    targetCell<-c()
    for( t in 1:length(currentTargets)){
      targetCell<-paste(targetCell, currentTargets[t], sep= ", ")
    }
    targetCell = gsub("^, ", "", targetCell)
    result[i,1]<-targetCell
    result[i,2]<-(colSums(drugData[,match(colnames(drugFilt[,i,drop=FALSE]), colnames(drugData)), drop=FALSE])*(1/(1 + length(currentTargets)^(-log(length(sigList)/nrow(testResult))))))
  }
  row.names(result)<-colnames(drugFilt)
  result<-result[order(result[,2]),]
  colnames(result)<-c("Drug_Targets", "Promiscuity_Score")
  return(result)
}