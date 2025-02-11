#' This takes dataframes with rna and metabolite info and performs pathway analysis
#' @param geneResult Dataframe with transcript IDs as rownames, transcript signifigance values (between 0 and 1) as the first column, and optionally transcript fold change values as the second column. Defaults to NULL
#' @param metaboliteResult Dataframe with metabolite IDs as rownames, metabolite signifigance values (between 0 and 1) as the first column, and optionally metabolite fold change values as the second column. Defaults to NULL
#' @param combine Logical indicated whether to perform combined transcirpt and metabolite analysis.  Defaults to FALSE
#' @param method Method of pathway analysis to perform. Options include "fisher.exact",  "EASE", "mean.significance", and "hypergeometric.test". Defaults to "fisher.exact.
#' @param geneCutoff Significance cutoff to use for transcript analysis
#' @param metaboliteCutoff Significance cutoff to use for metabolite analysis
#' @param pathwayType Pathway database to use for pathway analysis.  Options include "KEGG", "SMPDB", "Reactome", or "All".  Defaults to "All"
#' @param calculateFoldChange Logical indicating whether to perform cumulative fold change calculations.  If true, must have a second column with fold change numbers in geneResult and/or metaboliteResult input.  Defaults to TRUE.
#' @param countThreshold Integer indicating the minimum number of pathway members to require for a valid result. Defaults to 2.
#' @param p.adjustMethod Character indicating the type of multiple hypothesis correction to perform.  Defaults to "BH". See ?p.adjust() for details
#' @param plotResult Logical indicating whether to barplot negative log FDR values of pathway results.
#' @param plotNumber Integer indicating number of pathways, in order of significance, to include in barplot if plotResult == TRUE
#' @keywords IPA
#' @export
#' @return Returns a dataframe with pathway FDR corrected pathway significance and cumulative fold change for each pathway.
#' @examples
#' ## Parallel pathway analysis on individual datasets using Fisher's exact test.
#' data(kData)
#' data(rData)
#' ipa(geneResult=kData, metaboliteResult=rData, combine=FALSE, geneCutoff=0.1,
#'      metaboliteCutoff=0.219)
#' 
#' ##Integrated pathway analysis using Fisher's exact test.
#' data(kData)
#' data(rData)
#' ipa(geneResult=kData, metaboliteResult=rData, combine=TRUE, geneCutoff=0.1,
#'      metaboliteCutoff=0.219)
#' @note
#' ipa is powered by the following open source databases.  Commercial use and/or redistribution may restricted.  Please see respective terms of use pages and citations for more details.
#' 
#' KEGG
#' 
#' Terms of Use: http://www.kegg.jp/kegg/legal.html
#' 
#' Citations:
#' 
#' Kanehisa, M., Goto, S., Sato, Y., Kawashima, M., Furumichi, M., and Tanabe, M.  Data, information, knowledge and principle: back to metabolism in KEGG. Nucleic Acids Res. 42, D199-D205 (2014).
#' 
#' Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30 (2000).
#'
#' SMPDB
#' 
#' Terms of Use: http://smpdb.ca/about
#' 
#' Citations:
#' 
#' Wishart DS, Frolkis A, Knox C, et al. SMPDB: The Small Molecule Pathway Database. Nucleic Acids Res. 2010 Jan;38(Database issue):D480-7.
#' 
#' Jewison T, Su Y, Disfany FM, et al. SMPDB 2.0: Big Improvements to the Small Molecule Pathway Database Nucleic Acids Res. 2014 Jan;42(Database issue):D478-84.
#' 
#' Reactome
#' 
#' Terms of Use: http://www.reactome.org/pages/about/license-agreement/
#' 
#' Citations:
#' 
#' Croft et al. 2014 PMID: 24243840
#' 
#' Milacic et al. 2012 PMID:24213504

ipa<-function(geneResult=NULL, metaboliteResult=NULL, combine=FALSE, 
              method=c("fisher.exact", "EASE", "mean.significance", "hypergeometric"), 
              geneCutoff=NULL, metaboliteCutoff=NULL, 
              pathwayType=c("All", "KEGG", "SMPDB", "Reactome"), 
              calculateFoldChange=TRUE, countThreshold=2, p.adjustMethod="BH",
              plotResult = TRUE, plotNumber = 10){
  #Import database info
  data(genes, envir=environment())
  data(metabolites, envir=environment())
  data(pathIndex, envir=environment())
  pathIndex<-pathIndex
  #Subset database info by pathwayType
  if(pathwayType[1] == "All"){ 
    genes<-genes
    metabolites<-metabolites
  }
  else if(pathwayType[1] == "KEGG"){
    genes<-genes[,colnames(genes)[which(grepl("hsa", colnames(genes)))]]
    metabolites<-metabolites[,colnames(metabolites)[which(grepl("hsa", colnames(metabolites)))]]
  }
  else if(pathwayType[1] == "SMPDB"){
    genes<-genes[,colnames(genes)[which(grepl("SMP", colnames(genes)))]]
    metabolites<-metabolites[,colnames(metabolites)[which(grepl("SMP", colnames(metabolites)))]]
  }
  else if(pathwayType[1] == "Reactome"){
    genes<-genes[,colnames(genes)[which(grepl("REACT", colnames(genes)))]]
    metabolites<-metabolites[,colnames(metabolites)[which(grepl("REACT", colnames(metabolites)))]]
  }
  
  #establish result index
  geneR<-NULL
  metabR<-NULL
  
  #Funcions for Analysis:
  #Function for calculating pathway cumulative fold change efficiently
  calcFoldChange<-function(pathwayFilt, toORA){
    FCvalues<-toORA[row.names(pathwayFilt),]
    toFCsum<-FCvalues[,2]*pathwayFilt
    FCsum<-data.frame(apply(toFCsum, 2, function(x) sum(x, na.rm=T)))
    return(FCsum)
  }
  
  #Function for calculating Fisher's Exact Test
  calcFisherExact<-function(pathwayFilt, sigFilter, background, p.adjustMethod){
    #Determine number of significant markers overlapping with each pathway
    k = data.frame(colSums(pathwayFilt[sigFilter,], na.rm=TRUE))
    #Determine which background markers are found in specified pathways
    backgroundFilter <- row.names(pathwayFilt)%in%background
    #Determine number of background markers overlapping with each pathway
    m = data.frame(colSums(pathwayFilt[backgroundFilter,], na.rm=TRUE))
    #Perform fisher Exact Test
    N = sum(backgroundFilter, na.rm=TRUE)
    n = sum(sigFilter, na.rm=TRUE)
    #Calculate Fold Enrichment
    FE<-data.frame((k/n)/(m/N))
    colnames(FE)<-"Fold_Enrichment"
    toFish<-cbind(k,m,n,N,FE)
    #Function for calculating fisher exact test
    Fapply<-function(x){
      result<-c(fisher.test(as.matrix(data.frame(x=c(x[[1]], (x[[2]]-x[[1]])), y = c((x[[3]]-x[[1]]), (x[[4]]+x[[1]]-x[[3]]-x[[2]])))))$p.value, x[[1]], x[[2]], x[[5]]) 
      return(result)
    }
    result<-data.frame(t(apply(toFish, 1, function(x) Fapply(x))))
    
    
    #Perform tail-based (BH) Fdr analysis
    toFDR = data.frame(replace(result[,1,drop=FALSE], result[,1,drop=FALSE]>1, 1))
    FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
    row.names(FDR)<-row.names(result)
    #Create final dataframe
    finalResult = cbind(FDR, result)
    return(finalResult)
  }
  
  #Function for calculating EASE Score
  calcEASE<-function(pathwayFilt, sigFilter, background, p.adjustMethod){
    k = data.frame(colSums(pathwayFilt[sigFilter,], na.rm=TRUE))
    #Determine which background markers are found in specified gene pathways
    backgroundFilter <- row.names(pathwayFilt)%in%background
    #Determine number of background markers overlapping with each gene pathway
    m = data.frame(colSums(pathwayFilt[backgroundFilter,], na.rm=TRUE))
    N = sum(backgroundFilter, na.rm=TRUE)
    n = sum(sigFilter, na.rm=TRUE)
    #Calculate Fold Enrichment
    FE<-data.frame((k/n)/(m/N))
    colnames(FE)<-"Fold_Enrichment"
    toEASE<-cbind(k,m,n,N,FE)
    #Function for calculating EASE Score
    Eapply<-function(x){
      result<-c(fisher.test(as.matrix(data.frame(x=c(max(x[[1]]-1,0), (x[[2]]-x[[1]])), y = c((x[[3]]-x[[1]]), (x[[4]]+x[[1]]-x[[3]]-x[[2]])))))$p.value, x[[1]], x[[2]], x[[5]]) 
      return(result)
    }
    result<-data.frame(t(apply(toEASE, 1, function(x) Eapply(x))))
    
    #Perform tail-based (BH) Fdr analysis
    toFDR = data.frame(replace(result[,1], result[,1]>1, 1))
    FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
    row.names(FDR)<-row.names(result)
    #Create final dataframe
    finalResult = cbind(FDR, result)
    return(finalResult)
  }
  
  #Function for calculating Hypergeometric Test
  calcHypergeometric<-function(pathwayFilt, sigFilter, background, p.adjustMethod){
    #Determine number of significant markers overlapping with each pathway
    q = data.frame(colSums(pathwayFilt[sigFilter,], na.rm=TRUE))
    #Determine which background markers are found in specified pathways
    backgroundFilter <- row.names(pathwayFilt)%in%background
    #Determine number of background markers overlapping with each pathway
    m = data.frame(colSums(pathwayFilt[backgroundFilter,], na.rm=TRUE))
    N = sum(backgroundFilter, na.rm=TRUE)
    #Perform Hypergeometric Test
    k = sum(sigFilter, na.rm=TRUE)
    FE<-data.frame((q/k)/(m/N))
    toHG<-cbind(q,m,k,N,FE)
    Happly<-function(x){
      result<-c(phyper(q=x[[1]], m=x[[2]], k=x[[3]], n=x[[4]]-x[[2]], lower.tail=FALSE), x[[1]], x[[2]],x[[5]]) 
      return(result)
    }
    result<-data.frame(t(apply(toHG, 1, function(x) Happly(x))))
    
    
    #Perform tail-based (BH) Fdr analysis
    toFDR = data.frame(replace(result[,1], result[,1]>1 | result[,1]==0, 1))
    FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
    #Create final dataframe
    finalResult = cbind(FDR, result)
    return(finalResult)
  }
  
  #Function for calculating mean pathway significance
  calcMeanSig<-function(pathwayFilt, background, testResult){
    m = data.frame(colSums(pathwayFilt[background,], na.rm=TRUE))
    meanIndex<-testResult[row.names(pathwayFilt),]
    toMean<-meanIndex[,1]*pathwayFilt
    meanFrame<-data.frame(apply(toMean, 2, function(x) sum(x,na.rm=TRUE)/sum(x>0, na.rm=TRUE)))
    result<-cbind(meanFrame, m)
    return(result)
  }
  
  #Function for calculating gene analysis
  gene.test <- function(testResult, background = row.names(testResult), genes, cutoff, method, calculateFoldChange = TRUE){
    pathwayFilt<-genes[,-which(colSums(genes[background, ], na.rm=TRUE)<countThreshold)]
    if(method[1] == "fisher.exact"){
      if(calculateFoldChange == FALSE){
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(pathwayFilt)%in%sigList
        
        #Calculate Fisher's exact test and other pathway info
        finalResult<-calcFisherExact(pathwayFilt, sigFilter, background, p.adjustMethod)
        
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " genes out of ", length(background), " total genes to ", ncol(pathwayFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        
        colnames(finalResult)<- c("Pathway_Gene_FDR","Pathway_Gene_p_value", "Sig_Gene_Number", "Total_Gene_Number", "Gene_Fold_Enrichment")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
      }
      else{
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(pathwayFilt)%in%sigList
        
        #Calculate cumulative Fold Change for pathway directionality
        FCsum<-calcFoldChange(pathwayFilt, toORA)
        
        #Calculate Fisher's exact test and other pathway info
        finalResult<-calcFisherExact(pathwayFilt, sigFilter, background, p.adjustMethod)
        
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " genes out of ", length(background), " total genes to ", ncol(pathwayFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        
        #Add significant fold change information
        finalResult<-cbind(finalResult, FCsum)
        colnames(finalResult)<- c("Pathway_Gene_FDR","Pathway_Gene_p_value", "Sig_Gene_Number", "Total_Gene_Number", "Gene_Fold_Enrichment", "Gene_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
    }
    
    else if(method[1] == "EASE"){
      if(calculateFoldChange == FALSE){
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(pathwayFilt)%in%sigList
        #calculate EASE score and other pathway info
        finalResult<-calcEASE(pathwayFilt, sigFilter, background, p.adjustMethod)
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " genes out of ", length(background), " total genes to ", ncol(pathwayFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        colnames(finalResult)<- c("Pathway_Gene_FDR","Pathway_Gene_p_value", "Sig_Gene_Number", "Total_Gene_Number", "Gene_Fold_Enrichment")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
      }
      else{
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(pathwayFilt)%in%sigList
        
        #Calculate cumulative Fold Change for pathway directionality
        FCsum<-calcFoldChange(pathwayFilt, toORA)
        
        #calculate EASE score and other pathway info
        finalResult<-calcEASE(pathwayFilt, sigFilter, background, p.adjustMethod)
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " genes out of ", length(background), " total genes to ", ncol(pathwayFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        
        #Add significant fold change information
        finalResult<-cbind(finalResult, FCsum)
        colnames(finalResult)<- c("Pathway_Gene_FDR","Pathway_Gene_p_value", "Sig_Gene_Number", "Total_Gene_Number", "Gene_Fold_Enrichment", "Gene_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
    }
    
    else if(method[1] == "hypergeometric"){
      if(calculateFoldChange == FALSE){
        
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(pathwayFilt)%in%sigList
        #Calculate hypergeometric test and other pathway info
        finalResult<-calcHypergeometric(pathwayFilt, sigFilter, background, p.adjustMethod)
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " genes out of ", length(background), " total genes to ", ncol(pathwayFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        colnames(finalResult)<- c("Pathway_Gene_FDR","Pathway_Gene_p_value", "Sig_Gene_Number", "Total_Gene_Number", "Fold_Enrichment")
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
      else{
        
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(pathwayFilt)%in%sigList
        
        #Calculate cumulative Fold Change for pathway directionality
        FCsum<-calcFoldChange(pathwayFilt, toORA)
        
        #Calculate hypergeometric test and other pathway info
        finalResult<-calcHypergeometric(pathwayFilt, sigFilter, background, p.adjustMethod)
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " genes out of ", length(background), " total genes to ", ncol(pathwayFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        
        #Add significant fold change information
        finalResult<-cbind(finalResult, FCsum)
        colnames(finalResult)<- c("Pathway_Gene_FDR","Pathway_Gene_p_value", "Sig_Gene_Number", "Total_Gene_Number", "Fold_Enrichment", "Gene_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
      
    }
    
    else if(method[1] == "mean.significance"){
      if(calculateFoldChange == FALSE){
        #Calculate mean significance values and other pathway info
        result<-calcMeanSig(pathwayFilt, background, testResult)
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " genes out of ", length(background), " total genes to ", ncol(pathwayFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        colnames(result)=c("Mean_Pathway_Significance", "Sig_Gene_Number")
        result<-result[order(result[,1]),]
        return(result)
        
      }else{
        
        #Calculate pathway fold changes
        #Generate dataframe for FC analysis
        toORA <- testResult[order(testResult[,1]),]
        FCsum<-calcFoldChange(pathwayFilt, toORA)
        
        #Calculate mean significance values and other pathway info
        result<-calcMeanSig(pathwayFilt, background, testResult)
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " genes out of ", length(background), " total genes to ", ncol(pathwayFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        result<-cbind(result, FCsum)
        colnames(result)=c("Mean_Pathway_Significance", "Sig_Gene_Number", "Gene_Cumulative_Fold_Change")
        result<-result[order(result[,1]),]
        return(result)
      }
    }
  }
  
  #Function for calculating metabolite analysis
  metab.test <- function(testResult, background = row.names(testResult), metabolites, cutoff, method, calculateFoldChange = TRUE){
    pathwayFilt<-metabolites[,-which(colSums(metabolites[background, ], na.rm=TRUE)<countThreshold)]
    if(method[1] == "fisher.exact"){
      if(calculateFoldChange == FALSE) {
        
        toORA <- testResult[order(testResult[,1]),]
        
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        
        sigFilter <- row.names(pathwayFilt)%in%sigList
        
        #Calculate Fisher's exact test and other pathway info
        finalResult<-calcFisherExact(pathwayFilt, sigFilter, background, p.adjustMethod)
        
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " metabolites out of ", length(background), " total metabolites to ", ncol(pathwayFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        
        colnames(finalResult)<- c("Pathway_Metabolite_FDR", "Pathway_Metabolite_p_value", "Sig_Metabolite_Number", "Total_Metabolite_Number", "Metabolite_Fold_Enrichment")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
      }
      
      else {
        
        toORA <- testResult[order(testResult[,1]),]
        
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        
        sigFilter <- row.names(pathwayFilt)%in%sigList
        FCsum<-calcFoldChange(pathwayFilt, toORA)
        #Calculate Fisher's exact test and other pathway info
        finalResult<-calcFisherExact(pathwayFilt, sigFilter, background, p.adjustMethod)
        
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " metabolites out of ", length(background), " total metabolites to ", ncol(pathwayFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        
        #Add significant fold change information
        finalResult<-cbind(finalResult, FCsum)
        colnames(finalResult)<- c("Pathway_Metabolite_FDR","Pathway_Metabolite_p_value", "Sig_Metabolite_Number", "Total_Metabolite_Number", "Metabolite_Fold_Enrichment", "Metabolite_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
      }
    }
    else if(method[1] == "EASE"){
      if(calculateFoldChange == FALSE) {
        
        toORA <- testResult[order(testResult[,1]),]
        
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        
        sigFilter <- row.names(pathwayFilt)%in%sigList
        
        #calculate EASE score and other pathway info
        finalResult<-calcEASE(pathwayFilt, sigFilter, background, p.adjustMethod)
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " metabolites out of ", length(background), " total metabolites to ", ncol(pathwayFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        colnames(finalResult)<- c("Pathway_Metabolite_FDR", "Pathway_Metabolite_p_value", "Sig_Metabolite_Number", "Total_Metabolite_Number", "Metabolite_Fold_Enrichment")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
      }
      
      else {
        
        toORA <- testResult[order(testResult[,1]),]
        
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        
        sigFilter <- row.names(pathwayFilt)%in%sigList
        FCsum<-calcFoldChange(pathwayFilt, toORA)
        
        #calculate EASE score and other pathway info
        finalResult<-calcEASE(pathwayFilt, sigFilter, background, p.adjustMethod)
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " metabolites out of ", length(background), " total metabolites to ", ncol(pathwayFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        
        #Add significant fold change information
        finalResult<-cbind(finalResult, FCsum)
        colnames(finalResult)<- c("Pathway_Metabolite_FDR","Pathway_Metabolite_p_value", "Sig_Metabolite_Number", "Total_Metabolite_Number", "Metabolite_Fold_Enrichment", "Metabolite_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
      }
    }
    else if(method[1] == "hypergeometric"){
      if(calculateFoldChange == FALSE){
        
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(pathwayFilt)%in%sigList
        
        #Calculate hypergeometric test and other pathway info
        finalResult<-calcHypergeometric(pathwayFilt, sigFilter, background, p.adjustMethod)
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " metabolites out of ", length(background), " total metabolites to ", ncol(pathwayFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        colnames(finalResult)<- c("Pathway_Metabolite_FDR","Pathway_Metabolite_p_value", "Sig_Metabolite_Number", "Total_Metabolite_Number", "Fold_Enrichment")
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
      else{
        
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(pathwayFilt)%in%sigList
        
        #Calculate cumulative Fold Change for pathway directionality
        FCsum<-calcFoldChange(pathwayFilt, toORA)
        
        #Calculate hypergeometric test and other pathway info
        finalResult<-calcHypergeometric(pathwayFilt, sigFilter, background, p.adjustMethod)
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " metabolites out of ", length(background), " total metabolites to ", ncol(pathwayFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        
        #Add significant fold change information
        finalResult<-cbind(finalResult, FCsum)
        colnames(finalResult)<- c("Pathway_Metabolite_FDR","Pathway_Metabolite_p_value", "Sig_Metabolite_Number", "Total_Metabolite_Number", "Fold_Enrichment", "Metabolite_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
      
    }
    else if(method[1] == "mean.significance"){
      if(calculateFoldChange == FALSE){
        #Calculate mean significance values and other pathway info
        result<-calcMeanSig(pathwayFilt, background, testResult)
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " metabolites out of ", length(background), " total metabolites to ", ncol(pathwayFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        colnames(result)=c("Mean_Pathway_Significance", "Sig_Metabolite_Number")
        result<-result[order(result[,1]),]
        return(result)
        
      }else{
        
        #Calculate pathway fold changes
        #Generate dataframe for FC analysis
        toORA <- testResult[order(testResult[,1]),]
        FCsum<-calcFoldChange(pathwayFilt, toORA)
        
        #Calculate mean significance values and other pathway info
        result<-calcMeanSig(pathwayFilt, background, testResult)
        message(paste("Matched ", sum(row.names(pathwayFilt)%in%background, na.rm=TRUE), " metabolites out of ", length(background), " total metabolites to ", ncol(pathwayFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        result<-cbind(result, FCsum)
        colnames(result)=c("Mean_Pathway_Significance", "Sig_Metabolite_Number", "Metabolite_Cumulative_Fold_Change")
        result<-result[order(result[,1]),]
        return(result)
      }
    }
  }
  
  #Function for calculating combined gene-metabolite analysis
  combined.test <- function(gene_result=NULL, metab_result=NULL, calculateFoldChange = TRUE, method){
    #define fisher function
    Fisher.test <- function(p) {
      Xsq <- -2*sum(log(p))
      p.val <- pchisq(Xsq, df = 2*length(p),lower.tail=FALSE)
      return(p.value = p.val)
    }
    
    if(method[1]=="mean.significance"){
      if(calculateFoldChange == FALSE){  
        #filter gene_result for only pathways which include metabolites
        genes_w_metab<-row.names(gene_result)[which(row.names(gene_result)%in%row.names(metab_result))]
        gene_result<-gene_result[genes_w_metab, ]
        metab_result<-metab_result[genes_w_metab, ]
        
        result = matrix(nrow=nrow(metab_result), ncol=3)
        for(i in 1:nrow(metab_result)){
          result[i,1] = sum(gene_result[i,1]*(gene_result[i,2]/(gene_result[i,2]+metab_result[i,2])) + metab_result[i,1]*(metab_result[i,2]/(metab_result[i,2]+gene_result[i,2])))
          result[i,2] = gene_result[i,2]
          result[i,3] = metab_result[i,2]
        }
        row.names(result)<-row.names(gene_result)
        colnames(result)<- c("Mean_Combined_Pathway_Significance", "Total_Gene_Number", "Total_Metabolite_Number")
        
        result<-result[order(result[,1]),]
      }
      else{  
        #filter gene_result for only pathways which include metabolites
        genes_w_metab<-row.names(gene_result)[which(row.names(gene_result)%in%row.names(metab_result))]
        gene_result<-gene_result[genes_w_metab, ]
        metab_result<-metab_result[genes_w_metab, ]
        
        result = matrix(nrow=nrow(metab_result), ncol=5)
        for(i in 1:nrow(metab_result)){
          result[i,1] = sum(gene_result[i,1]*(gene_result[i,2]/(gene_result[i,2]+metab_result[i,2])) + metab_result[i,1]*(metab_result[i,2]/(metab_result[i,2]+gene_result[i,2])))
          result[i,2] = gene_result[i,2]
          result[i,3] = gene_result[i,3]
          result[i,4] = metab_result[i,2]
          result[i,5] = metab_result[i,3]
        }
        row.names(result)<-row.names(gene_result)
        colnames(result)<- c("Mean_Combined_Pathway_Significance", "Total_Gene_Number", "Gene_Cumulative_Fold_Change","Total_Metabolite_Number", "Metabolite_Cumulative_Fold_Change")
        
        result<-result[order(result[,1]),]
      }
    }
    else{
      if(calculateFoldChange == FALSE){  
        #filter gene_result for only pathways which include metabolites
        genes_w_metab<-row.names(gene_result)[which(row.names(gene_result)%in%row.names(metab_result))]
        gene_result<-gene_result[genes_w_metab, ]
        metab_result<-metab_result[genes_w_metab, ]
        
        #perform fisher method for joint p-values on all pathways
        result = matrix(nrow=length(metab_result[,1]), ncol=9)
        for(i in 1:length(metab_result[,1])){
          result[i,1] = Fisher.test(p = c(gene_result[i,2], metab_result[i,2]))
          result[i,2] = gene_result[i,2]
          result[i,3] = gene_result[i,3]
          result[i,4] = gene_result[i,4]
          result[i,5] = gene_result[i,5]
          result[i,6] = metab_result[i,2]
          result[i,7] = metab_result[i,3]
          result[i,8] = metab_result[i,4]
          result[i,9] = metab_result[i,5]
          
        }
        
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1 | result[,1]==0, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #label row and column names appropriately
        row.names(finalResult)<-row.names(metab_result)
        colnames(finalResult) <- c("Pathway_Combined_FDR","Pathway_Combined_pValue",colnames(gene_result[-1]), colnames(metab_result[-1]))
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
      else{
        #filter gene_result for only pathways which include metabolites
        genes_w_metab<-row.names(gene_result)[which(row.names(gene_result)%in%row.names(metab_result))]
        gene_result<-gene_result[genes_w_metab, ]
        metab_result<-metab_result[genes_w_metab, ]
        
        #perform fisher method for joint p-values on all pathways
        result = matrix(nrow=length(metab_result[,1]), ncol=11)
        for(i in 1:length(metab_result[,1])){
          result[i,1] = Fisher.test(p = c(gene_result[i,2], metab_result[i,2]))
          result[i,2] = gene_result[i,2]
          result[i,3] = gene_result[i,3]
          result[i,4] = gene_result[i,4]
          result[i,5] = gene_result[i,5]
          result[i,6] = gene_result[i,6]
          result[i,7] = metab_result[i,2]
          result[i,8] = metab_result[i,3]
          result[i,9] = metab_result[i,4]
          result[i,10] = metab_result[i,5]
          result[i,11] = metab_result[i,6]
          
        }
        
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1 | result[,1]==0, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #label row and column names appropriately
        row.names(finalResult)<-row.names(metab_result)
        colnames(finalResult) <- c("Pathway_Combined_FDR","Pathway_Combined_pValue",colnames(gene_result[-1]), colnames(metab_result[-1]))
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
    }
  }
  
  #Function for plotting analysis result
  plotIPA = function(x, cutoff=10){
    
    toPlot <- x[order(x[,3]),][1:cutoff,]
    Annotation = factor(row.names(toPlot), levels = unique(row.names(toPlot)))
    if(length(grep("Fold_Change", colnames(toPlot)))>0){
      mid=0
      p1 = ggplot2::qplot(Annotation, -log(toPlot[,2], base=2), geom= "bar", stat="identity", binwidth= 0.00001, fill = toPlot[,grep("Fold_Change", colnames(toPlot))[1]])
      p2 = p1 + ggplot2::scale_fill_gradient2(midpoint = mid, "Fold Change", limits = c(min(x[,grep("Fold_Change", colnames(toPlot))[1]]),max(x[,grep("Fold_Change", colnames(toPlot))[1]])), low = "blue",mid = "black", high = "red", space = "Lab" )
      p3 = p2 + ggplot2::ylab("Negative Log FDR")+ggplot2::xlab("")+ggplot2::coord_flip() + ggplot2::theme( panel.background = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black", size = 1), axis.text = ggplot2::element_text(face = "bold", colour = "black"),  axis.ticks = ggplot2::element_line(colour = "black"), axis.title = ggplot2::element_text(face="bold", colour = "black"))
      return(p3)
    }
    else{
      p1 = ggplot2::qplot(Annotation, -log(toPlot[,2], base=2), geom= "bar", stat="identity", binwidth= 0.00001)
      p2 = p1 + ggplot2::ylab("Negative Log FDR")+ggplot2::xlab("")+ggplot2::coord_flip() + ggplot2::theme( panel.background = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black", size = 1), axis.text = ggplot2::element_text(face = "bold", colour = "black"),  axis.ticks = ggplot2::element_line(colour = "black"), axis.title = ggplot2::element_text(face="bold", colour = "black"))
      return(p2)
    }
  }
  
  if(!is.null(geneResult)){
    
    gene<-gene.test(testResult=geneResult, method=method, genes=genes, cutoff=geneCutoff, calculateFoldChange=calculateFoldChange)
    geneR<-gene[which(gene[,"Sig_Gene_Number"]>=countThreshold),]
    geneR<-merge(pathIndex, geneR, by.x="row.names", by.y="row.names")
    geneR<-data.frame(geneR[,-1], row.names=geneR[,1])
    geneR<-geneR[order(geneR[,3]),]
  }
  
  if(!is.null(metaboliteResult)){
    
    metab<-metab.test(testResult=metaboliteResult, method=method, metabolites=metabolites, cutoff=metaboliteCutoff, calculateFoldChange=calculateFoldChange)
    metabR<-metab[which(metab[,"Sig_Metabolite_Number"]>=countThreshold),]
    metabR<-merge(pathIndex, metabR, by.x="row.names", by.y="row.names")
    metabR<-data.frame(metabR[,-1], row.names=metabR[,1])
    metabR<-metabR[order(metabR[,3]),]
  }
  
  result<-list(geneR, metabR)
  names(result)<-c("Gene_Analysis", "Metabolite_Analysis")
  
  if(combine==TRUE){
    combineR<-combined.test(gene_result=geneR[,-1], metab_result=metabR[,-1], calculateFoldChange=calculateFoldChange, method=method)
    combineR<-merge(pathIndex, combineR, by.x="row.names", by.y="row.names")
    combineR<-data.frame(combineR[,-1], row.names=combineR[,1])
    combineR<-combineR[order(combineR[,3]),]
    if(plotResult==TRUE){
      PathPlot<-plotIPA(combineR, cutoff=plotNumber)
      print(PathPlot)
    }
    return(combineR)
  }
  else if(combine==FALSE){
    if (plotResult==TRUE){
      if(!is.null(result[[1]])){
        PathPlot<-plotIPA(result[[1]], cutoff=plotNumber)
        print(PathPlot)
      }
      if(!is.null(result[[2]])){
        PathPlot<-plotIPA(result[[2]], cutoff=plotNumber)
        print(PathPlot)
      }
    }
    return(result)
  }
}
