#' This takes dataframes with rna and metabolite info and performs pathway analysis
#' @param geneResult Dataframe with transcript IDs as rownames, transcript signifigance values (between 0 and 1) as the first column, and optionally transcript fold change values as the second column. Defaults to NULL
#' @param proteinResult Dataframe with Uniprot Accession IDs as rownames, Protein signifigance values (between 0 and 1) as the first column, and optionally protein fold change values as the second column. Defaults to NULL
#' @param metaboliteResult Dataframe with metabolite IDs as rownames, metabolite signifigance values (between 0 and 1) as the first column, and optionally metabolite fold change values as the second column. Defaults to NULL
#' @param combine Logical indicated whether to perform combined transcirpt and metabolite analysis.  Defaults to FALSE
#' @param method Method of pathway analysis to perform. Options include "fisher.exact",  "EASE", "mean.significance", and "hypergeometric.test". Defaults to "fisher.exact.
#' @param geneCutoff Significance cutoff to use for transcript analysis
#' @param proteinCutoff Significance cutoff to use for protein analysis
#' @param metaboliteCutoff Significance cutoff to use for metabolite analysis
#' @param pathwayType Pathway database to use for pathway analysis.  Options include "KEGG", "SMPDB", "Reactome", or "All".  Defaults to "All"
#' @param calculateFoldChange Logical indicating whether to perform cumulative fold change calculations.  If true, must have a second column with fold change numbers in geneResult and/or metaboliteResult input.  Defaults to TRUE.
#' @param countThreshold Integer indicating the minimum number of pathway members to require for a valid result. Defaults to 2.
#' @param p.adjustMethod Character indicating the type of multiple hypothesis correction to perform.  Defaults to "BH". See ?p.adjust() for details
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


ipa<-function(geneResult=NULL, proteinResult=NULL, metaboliteResult=NULL, combine=FALSE, method=c("fisher.exact", "EASE", "mean.significance", "hypergeometric"), geneCutoff=NULL, proteinCutoff=NULL, metaboliteCutoff=NULL, pathwayType=c("All", "KEGG", "SMPDB", "Reactome"), calculateFoldChange=TRUE, countThreshold=2, p.adjustMethod="BH"){
  
  data(genes, envir=environment())
  data(metabolites, envir=environment())
  data(proteins, envir=environment())
  data(pathIndex, envir=environment())
  pathIndex<-pathIndex
  
  if(pathwayType[1] == "All"){ 
    genes<-genes
    metabolites<-metabolites
    proteins<-proteins
  }
  if(pathwayType[1] == "KEGG"){
    genes<-genes[,colnames(genes)[which(grepl("hsa", colnames(genes)))]]
    metabolites<-metabolites[,colnames(metabolites)[which(grepl("hsa", colnames(metabolites)))]]
    proteins<-proteins[,colnames(metabolites)[which(grepl("hsa", colnames(metabolites)))]]    
  }
  if(pathwayType[1] == "SMPDB"){
    genes<-genes[,colnames(genes)[which(grepl("SMP", colnames(genes)))]]
    metabolites<-metabolites[,colnames(metabolites)[which(grepl("SMP", colnames(metabolites)))]]
    proteins<-proteins[,colnames(metabolites)[which(grepl("SMP", colnames(metabolites)))]]
  }
  if(pathwayType[1] == "Reactome"){
    genes<-genes[,colnames(genes)[which(grepl("REACT", colnames(genes)))]]
    metabolites<-metabolites[,colnames(metabolites)[which(grepl("REACT", colnames(metabolites)))]]
    proteins<-proteins[,colnames(metabolites)[which(grepl("REACT", colnames(metabolites)))]]
  }
  
  #establish result index
  geneR<-NULL
  proteinR<-NULL
  metabR<-NULL
  
  
  gene.test <- function(testResult, background = row.names(testResult), genes, cutoff, method, calculateFoldChange = TRUE){
    genesFilt<-genes[,-which(colSums(genes[background, ], na.rm=TRUE)<countThreshold)]
    if(method[1] == "fisher.exact"){
      if(calculateFoldChange == FALSE){
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(genesFilt)%in%sigList
        #Determine number of significant markers overlapping with each gene pathway
        k = data.frame(colSums(genesFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified gene pathways
        backgroundFilter <- row.names(genesFilt)%in%background
        #Determine number of background markers overlapping with each gene pathway
        m = data.frame(colSums(genesFilt[backgroundFilter,], na.rm=TRUE))
        
        #Perform fisher Exact Test
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " genes out of ", length(background), " total genes to ", ncol(genesFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        n = sum(sigFilter, na.rm=TRUE)
        #Calculate Fold Enrichment
        FE<-data.frame((k/n)/(m/N))
        result = data.frame()
        for (i in 1:length(k[,1])){
          result[i,1] = fisher.test(as.matrix(data.frame(x = c(k[i,], (m[i,]-k[i,])), y = c((n-k[i,]), (N+k[i,]-n-m[i,])))))$p.value
          result[i,2] = k[i,1]
          result[i,3] = m[i,1]
          result[i,4] = FE[i,1]
        }
        
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(k)
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
        sigFilter <- row.names(genesFilt)%in%sigList
        
        #Calculate cumulative Fold Change for pathway directionality
        toFCsum = data.frame(genesFilt[sigFilter,])
        toFCsum<-toFCsum[,-which(colSums(toFCsum, na.rm=TRUE)<countThreshold)]
        FCsum = data.frame()
        for(i in 1:ncol(toFCsum)){
          filter = row.names(toFCsum[which(toFCsum[,i] == 1),])
          FCsum[i,1] = sum(toORA[,2][which(row.names(toORA)%in%filter)])
        }
        row.names(FCsum)<-colnames(toFCsum)
        
        #Determine number of significant markers overlapping with each gene pathway
        k = data.frame(colSums(genesFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified gene pathways
        backgroundFilter <- row.names(genesFilt)%in%background
        #Determine number of background markers overlapping with each gene pathway
        m = data.frame(colSums(genesFilt[backgroundFilter,], na.rm=TRUE))
        #Perform fisher Exact Test
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " genes out of ", length(background), " total genes to ", ncol(genesFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        n = sum(sigFilter, na.rm=TRUE)
        #Calculate Fold Enrichment
        FE<-data.frame((k/n)/(m/N))
        result = data.frame()
        for (i in 1:length(k[,1])){
          result[i,1] = fisher.test(as.matrix(data.frame(x = c(k[i,], (m[i,]-k[i,])), y = c((n-k[i,]), (N+k[i,]-n-m[i,])))))$p.value
          result[i,2] = k[i,1]
          result[i,3] = m[i,1]
          result[i,4] = FE[i,1]
        }
        
        
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(k)
        
        #Add significant fold change information
        finalResult<-merge(finalResult, FCsum, by.x="row.names", by.y="row.names", all=TRUE)
        finalResult<-data.frame(finalResult[,-1], row.names=finalResult[,1])
        
        colnames(finalResult)<- c("Pathway_Gene_FDR","Pathway_Gene_p_value", "Sig_Gene_Number", "Total_Gene_Number", "Gene_Fold_Enrichment", "Gene_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
    }
    
    if(method[1] == "EASE"){
      if(calculateFoldChange == FALSE){
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(genesFilt)%in%sigList
        #Determine number of significant markers overlapping with each gene pathway
        k = data.frame(colSums(genesFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified gene pathways
        backgroundFilter <- row.names(genesFilt)%in%background
        #Determine number of background markers overlapping with each gene pathway
        m = data.frame(colSums(genesFilt[backgroundFilter,], na.rm=TRUE))
        
        #Perform fisher Exact Test
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " genes out of ", length(background), " total genes to ", ncol(genesFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        n = sum(sigFilter, na.rm=TRUE)
        #Calculate Fold Enrichment
        FE<-data.frame((k/n)/(m/N))
        result = data.frame()
        for (i in 1:length(k[,1])){
          result[i,1] = fisher.test(as.matrix(data.frame(x = c(max(k[i,]-1,0), (m[i,]-k[i,])), y = c((n-k[i,]), (N+k[i,]-n-m[i,])))))$p.value
          result[i,2] = k[i,1]
          result[i,3] = m[i,1]
          result[i,4] = FE[i,1]
        }
        
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(k)
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
        sigFilter <- row.names(genesFilt)%in%sigList
        
        #Calculate cumulative Fold Change for pathway directionality
        toFCsum = data.frame(genesFilt[sigFilter,])
        toFCsum<-toFCsum[,-which(colSums(toFCsum, na.rm=TRUE)<countThreshold)]
        FCsum = data.frame()
        for(i in 1:ncol(toFCsum)){
          filter = row.names(toFCsum[which(toFCsum[,i] == 1),])
          FCsum[i,1] = sum(toORA[,2][which(row.names(toORA)%in%filter)])
        }
        row.names(FCsum)<-colnames(toFCsum)
        
        #Determine number of significant markers overlapping with each gene pathway
        k = data.frame(colSums(genesFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified gene pathways
        backgroundFilter <- row.names(genesFilt)%in%background
        #Determine number of background markers overlapping with each gene pathway
        m = data.frame(colSums(genesFilt[backgroundFilter,], na.rm=TRUE))
        #Perform fisher Exact Test
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " genes out of ", length(background), " total genes to ", ncol(genesFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        n = sum(sigFilter, na.rm=TRUE)
        #Calculate Fold Enrichment
        FE<-data.frame((k/n)/(m/N))
        result = data.frame()
        for (i in 1:length(k[,1])){
          result[i,1] = fisher.test(as.matrix(data.frame(x = c(max(k[i,]-1,0), (m[i,]-k[i,])), y = c((n-k[i,]), (N+k[i,]-n-m[i,])))))$p.value
          result[i,2] = k[i,1]
          result[i,3] = m[i,1]
          result[i,4] = FE[i,1]
        }
        
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(k)
        
        #Add significant fold change information
        finalResult<-merge(finalResult, FCsum, by.x="row.names", by.y="row.names", all=TRUE)
        finalResult<-data.frame(finalResult[,-1], row.names=finalResult[,1])
        
        colnames(finalResult)<- c("Pathway_Gene_FDR","Pathway_Gene_p_value", "Sig_Gene_Number", "Total_Gene_Number", "Gene_Fold_Enrichment", "Gene_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
    }
    
    if(method[1] == "hypergeometric"){
      if(calculateFoldChange == FALSE){
        
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(genesFilt)%in%sigList
        
        #Determine number of significant markers overlapping with each gene pathway
        q = data.frame(colSums(genesFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified gene pathways
        backgroundFilter <- row.names(genesFilt)%in%background
        #Determine number of background markers overlapping with each gene pathway
        m = data.frame(colSums(genesFilt[backgroundFilter,], na.rm=TRUE))
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " genes out of ", length(background), " total genes to ", ncol(genesFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        #Perform Hypergeometric Test
        k = sum(sigFilter, na.rm=TRUE)
        result = data.frame()
        for (i in 1:length(q[,1])){
          result[i,1] = phyper(q = q[i,],  m=m[i,], k = k, n = sum(backgroundFilter, na.rm=TRUE)-m[i,], lower.tail = FALSE)
          result[i,2] = q[i,1]
          result[i,3] = m[i,1]
        }
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1 | result[,1]==0, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(q)
        colnames(finalResult)<- c("Pathway_Gene_FDR","Pathway_Gene_p_value", "Sig_Gene_Number", "Total_Gene_Number")
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
      else{
        
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(genesFilt)%in%sigList
        
        #Calculate cumulative Fold Change for pathway directionality
        toFCsum = data.frame(genesFilt[sigFilter,])
        toFCsum<-toFCsum[,-which(colSums(toFCsum, na.rm=TRUE)<countThreshold)]
        FCsum = data.frame()
        for(i in 1:ncol(toFCsum)){
          filter = row.names(toFCsum[which(toFCsum[,i] == 1),])
          FCsum[i,1] = sum(toORA[,2][which(row.names(toORA)%in%filter)])
        }
        row.names(FCsum)<-colnames(toFCsum)
        
        #Determine number of significant markers overlapping with each gene pathway
        q = data.frame(colSums(genesFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified gene pathways
        backgroundFilter <- row.names(genesFilt)%in%background
        #Determine number of background markers overlapping with each gene pathway
        m = data.frame(colSums(genesFilt[backgroundFilter,], na.rm=TRUE))
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " genes out of ", length(background), " total genes to ", ncol(genesFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        #Perform Hypergeometric Test
        k = sum(sigFilter, na.rm=TRUE)
        result = data.frame()
        for (i in 1:length(q[,1])){
          result[i,1] = phyper(q = q[i,],  m=m[i,], k = k, n = sum(backgroundFilter, na.rm=TRUE)-m[i,], lower.tail = FALSE)
          result[i,2] = q[i,1]
          result[i,3] = m[i,1]
        }
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1 | result[,1]==0, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(q)
        
        
        #Add significant fold change information
        finalResult<-merge(finalResult, FCsum, by.x="row.names", by.y="row.names", all=TRUE)
        finalResult<-data.frame(finalResult[,-1], row.names=finalResult[,1])
        
        colnames(finalResult)<- c("Pathway_Gene_FDR","Pathway_Gene_p_value", "Sig_Gene_Number", "Total_Gene_Number", "Gene_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
      
    }
    
    if(method[1] == "mean.significance"){
      if(calculateFoldChange == FALSE){
        
        #Calculate number of genes mapped to pathways
        N = sum(row.names(genesFilt)%in%background, na.rm=TRUE)
        m = data.frame(colSums(genesFilt[background,], na.rm=TRUE))
        message(paste("Matched ", N, " genes out of ", length(background), " total genes to ", ncol(genesFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        message("Calculating mean pathway significance values.  This method is slower than the others.  To speed up, specify pathway type, increase count threshold, or change method.")
        
        #Calculate mean signifcance values for each pathway
        meanFrame = data.frame()
        for(i in 1:ncol(genesFilt)){
          filter = row.names(genesFilt[which(genesFilt[,i] == 1),])
          meanFrame[i,1] = mean(testResult[,1][which(row.names(testResult)%in%filter)], na.rm=TRUE)
        }
        row.names(meanFrame)<-colnames(genesFilt)
        
        #Compile results into data frame
        result = data.frame()
        for (i in 1:length(meanFrame[,1])){
          result[i,1] = meanFrame[i,1]
          result[i,2] = m[i,1]
        }
        
        row.names(result)=row.names(m)
        colnames(result)=c("Mean_Pathway_Significance", "Sig_Gene_Number")
        result<-result[order(result[,1]),]
        
        return(result)
        
      }else{
        
        #Calculate number of genes mapped to pathways
        N = sum(row.names(genesFilt)%in%background, na.rm=TRUE)
        m = data.frame(colSums(genesFilt[background,], na.rm=TRUE))
        message(paste("Matched ", N, " genes out of ", length(background), " total genes to ", ncol(genesFilt), " pathways out of ", ncol(genes), " total pathways with a count threshold of ", countThreshold, sep=""))
        message("Calculating mean pathway significance values.  This method is slower than the others.  To speed up, specify pathway type, increase count threshold, or change method.")
        
        #Calculate mean signifcance values for each pathway
        meanFrame = data.frame()
        for(i in 1:ncol(genesFilt)){
          filter = row.names(genesFilt[which(genesFilt[,i] == 1),])
          meanFrame[i,1] = mean(testResult[,1][which(row.names(testResult)%in%filter)], na.rm=TRUE)
        }
        row.names(meanFrame)<-colnames(genesFilt)
        
        #Calculate pathway fold changes
        FCsum = data.frame()
        for(i in 1:ncol(genesFilt)){
          filter = row.names(genesFilt[which(genesFilt[,i] == 1),])
          FCsum[i,1] = sum(testResult[,2][which(row.names(testResult)%in%filter)], na.rm=TRUE)
        }
        row.names(FCsum)<-colnames(genesFilt)
        
        
        #Compile results into data frame
        result = data.frame()
        for (i in 1:length(meanFrame[,1])){
          result[i,1] = meanFrame[i,1]
          result[i,2] = m[i,1]
          result[i,3] = FCsum[i,1]
        }
        
        row.names(result)=row.names(m)
        colnames(result)=c("Mean_Pathway_Significance", "Sig_Gene_Number", "Gene_Cumulative_Fold_Change")
        result<-result[order(result[,1]),]
        
        return(result)
      }
    }
  }
  
  protein.test <- function(testResult, background = row.names(testResult), proteins, cutoff, method, calculateFoldChange = TRUE){
    proteinFilt<-proteins[,-which(colSums(proteins[background, ], na.rm=TRUE)<countThreshold)]
    if(method[1] == "fisher.exact"){
      if(calculateFoldChange == FALSE){
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified protein pathways
        sigFilter <- row.names(proteinFilt)%in%sigList
        #Determine number of significant markers overlapping with each protein pathway
        k = data.frame(colSums(proteinFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified protein pathways
        backgroundFilter <- row.names(proteinFilt)%in%background
        #Determine number of background markers overlapping with each protein pathway
        m = data.frame(colSums(proteinFilt[backgroundFilter,], na.rm=TRUE))
        #Perform fisher Exact Test
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " proteins out of ", length(background), " total proteins to ", ncol(proteinFilt), " pathways out of ", ncol(proteins), " total pathways with a count threshold of ", countThreshold, sep=""))
        n = sum(sigFilter, na.rm=TRUE)
        #Calculate Fold Enrichment
        FE<-data.frame((k/n)/(m/N))
        result = data.frame()
        for (i in 1:length(k[,1])){
          result[i,1] = fisher.test(as.matrix(data.frame(x = c(k[i,], (m[i,]-k[i,])), y = c((n-k[i,]), (N+k[i,]-n-m[i,])))))$p.value
          result[i,2] = k[i,1]
          result[i,3] = m[i,1]
          result[i,4] = FE[i,1]
        }
        
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(k)
        colnames(finalResult)<- c("Pathway_protein_FDR","Pathway_protein_p_value", "Sig_protein_Number", "Total_protein_Number", "Protein_Fold_Enrichment")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
      }
      else{
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified protein pathways
        sigFilter <- row.names(proteinFilt)%in%sigList
        
        #Calculate cumulative Fold Change for pathway directionality
        toFCsum = data.frame(proteinFilt[sigFilter,])
        toFCsum<-toFCsum[,-which(colSums(toFCsum, na.rm=TRUE)<countThreshold)]
        FCsum = data.frame()
        for(i in 1:ncol(toFCsum)){
          filter = row.names(toFCsum[which(toFCsum[,i] == 1),])
          FCsum[i,1] = sum(toORA[,2][which(row.names(toORA)%in%filter)])
        }
        row.names(FCsum)<-colnames(toFCsum)
        
        
        #Determine number of significant markers overlapping with each protein pathway
        k = data.frame(colSums(proteinFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified protein pathways
        backgroundFilter <- row.names(proteinFilt)%in%background
        #Determine number of background markers overlapping with each protein pathway
        m = data.frame(colSums(proteinFilt[backgroundFilter,], na.rm=TRUE))
        #Perform fisher Exact Test
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " proteins out of ", length(background), " total proteins to ", ncol(proteinFilt), " pathways out of ", ncol(proteins), " total pathways with a count threshold of ", countThreshold, sep=""))
        n = sum(sigFilter, na.rm=TRUE)
        #Calculate Fold Enrichment
        FE<-data.frame((k/n)/(m/N))
        result = data.frame()
        for (i in 1:length(k[,1])){
          result[i,1] = fisher.test(as.matrix(data.frame(x = c(k[i,], (m[i,]-k[i,])), y = c((n-k[i,]), (N+k[i,]-n-m[i,])))))$p.value
          result[i,2] = k[i,1]
          result[i,3] = m[i,1]
          result[i,4] = FE[i,1]
        }
        
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(k)
        
        #Add significant fold change information
        finalResult<-merge(finalResult, FCsum, by.x="row.names", by.y="row.names", all=TRUE)
        finalResult<-data.frame(finalResult[,-1], row.names=finalResult[,1])
        
        colnames(finalResult)<- c("Pathway_protein_FDR","Pathway_protein_p_value", "Sig_protein_Number", "Total_protein_Number", "Protein_Fold_Enrichment", "Protein_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
    }
    if(method[1] == "EASE"){
      if(calculateFoldChange == FALSE){
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified protein pathways
        sigFilter <- row.names(proteinFilt)%in%sigList
        #Determine number of significant markers overlapping with each protein pathway
        k = data.frame(colSums(proteinFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified protein pathways
        backgroundFilter <- row.names(proteinFilt)%in%background
        #Determine number of background markers overlapping with each protein pathway
        m = data.frame(colSums(proteinFilt[backgroundFilter,], na.rm=TRUE))
        #Perform fisher Exact Test
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " proteins out of ", length(background), " total proteins to ", ncol(proteinFilt), " pathways out of ", ncol(proteins), " total pathways with a count threshold of ", countThreshold, sep=""))
        n = sum(sigFilter, na.rm=TRUE)
        #Calculate Fold Enrichment
        FE<-data.frame((k/n)/(m/N))
        result = data.frame()
        for (i in 1:length(k[,1])){
          result[i,1] = fisher.test(as.matrix(data.frame(x = c(max(k[i,]-1,0), (m[i,]-k[i,])), y = c((n-k[i,]), (N+k[i,]-n-m[i,])))))$p.value
          result[i,2] = k[i,1]
          result[i,3] = m[i,1]
          result[i,4] = FE[i,1]
        }
        
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(k)
        colnames(finalResult)<- c("Pathway_protein_FDR","Pathway_protein_p_value", "Sig_protein_Number", "Total_protein_Number", "Protein_Fold_Enrichment")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
      }
      else{
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified protein pathways
        sigFilter <- row.names(proteinFilt)%in%sigList
        
        #Calculate cumulative Fold Change for pathway directionality
        toFCsum = data.frame(proteinFilt[sigFilter,])
        toFCsum<-toFCsum[,-which(colSums(toFCsum, na.rm=TRUE)<countThreshold)]
        FCsum = data.frame()
        for(i in 1:ncol(toFCsum)){
          filter = row.names(toFCsum[which(toFCsum[,i] == 1),])
          FCsum[i,1] = sum(toORA[,2][which(row.names(toORA)%in%filter)])
        }
        row.names(FCsum)<-colnames(toFCsum)
        
        
        #Determine number of significant markers overlapping with each protein pathway
        k = data.frame(colSums(proteinFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified protein pathways
        backgroundFilter <- row.names(proteinFilt)%in%background
        #Determine number of background markers overlapping with each protein pathway
        m = data.frame(colSums(proteinFilt[backgroundFilter,], na.rm=TRUE))
        #Perform fisher Exact Test
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " proteins out of ", length(background), " total proteins to ", ncol(proteinFilt), " pathways out of ", ncol(proteins), " total pathways with a count threshold of ", countThreshold, sep=""))
        n = sum(sigFilter, na.rm=TRUE)
        #Calculate Fold Enrichment
        FE<-data.frame((k/n)/(m/N))
        result = data.frame()
        for (i in 1:length(k[,1])){
          result[i,1] = fisher.test(as.matrix(data.frame(x = c(max(k[i,]-1,0), (m[i,]-k[i,])), y = c((n-k[i,]), (N+k[i,]-n-m[i,])))))$p.value
          result[i,2] = k[i,1]
          result[i,3] = m[i,1]
          result[i,4] = FE[i,1]
        }
        
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(k)
        
        #Add significant fold change information
        finalResult<-merge(finalResult, FCsum, by.x="row.names", by.y="row.names", all=TRUE)
        finalResult<-data.frame(finalResult[,-1], row.names=finalResult[,1])
        
        colnames(finalResult)<- c("Pathway_protein_FDR","Pathway_protein_p_value", "Sig_protein_Number", "Total_protein_Number", "Protein_Fold_Enrichment", "Protein_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
    }
    if(method[1] == "hypergeometric"){
      if(calculateFoldChange == FALSE){
        
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(proteinFilt)%in%sigList
        
        #Determine number of significant markers overlapping with each gene pathway
        q = data.frame(colSums(proteinFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified gene pathways
        backgroundFilter <- row.names(proteinFilt)%in%background
        #Determine number of background markers overlapping with each gene pathway
        m = data.frame(colSums(proteinFilt[backgroundFilter,], na.rm=TRUE))
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " proteins out of ", length(background), " total proteins to ", ncol(proteinFilt), " pathways out of ", ncol(proteins), " total pathways with a count threshold of ", countThreshold, sep=""))
        #Perform Hypergeometric Test
        k = sum(sigFilter, na.rm=TRUE)
        result = data.frame()
        for (i in 1:length(q[,1])){
          result[i,1] = phyper(q = q[i,],  m=m[i,], k = k, n = sum(backgroundFilter, na.rm=TRUE)-m[i,], lower.tail = FALSE)
          result[i,2] = q[i,1]
          result[i,3] = m[i,1]
        }
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1 | result[,1]==0, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(q)
        colnames(finalResult)<- c("Pathway_Protein_FDR","Pathway_Protein_p_value", "Sig_protein_Number", "Total_Protein_Number")
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
      else{
        
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(proteinFilt)%in%sigList
        
        #Calculate cumulative Fold Change for pathway directionality
        toFCsum = data.frame(proteinFilt[sigFilter,])
        toFCsum<-toFCsum[,-which(colSums(toFCsum, na.rm=TRUE)<countThreshold)]
        FCsum = data.frame()
        for(i in 1:ncol(toFCsum)){
          filter = row.names(toFCsum[which(toFCsum[,i] == 1),])
          FCsum[i,1] = sum(toORA[,2][which(row.names(toORA)%in%filter)])
        }
        row.names(FCsum)<-colnames(toFCsum)
        
        #Determine number of significant markers overlapping with each gene pathway
        q = data.frame(colSums(proteinFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified gene pathways
        backgroundFilter <- row.names(proteinFilt)%in%background
        #Determine number of background markers overlapping with each gene pathway
        m = data.frame(colSums(proteinFilt[backgroundFilter,], na.rm=TRUE))
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " proteins out of ", length(background), " total proteins to ", ncol(proteinFilt), " pathways out of ", ncol(proteins), " total pathways with a count threshold of ", countThreshold, sep=""))
        #Perform Hypergeometric Test
        k = sum(sigFilter, na.rm=TRUE)
        result = data.frame()
        for (i in 1:length(q[,1])){
          result[i,1] = phyper(q = q[i,],  m=m[i,], k = k, n = sum(backgroundFilter, na.rm=TRUE)-m[i,], lower.tail = FALSE)
          result[i,2] = q[i,1]
          result[i,3] = m[i,1]
        }
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1 | result[,1]==0, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(q)
        
        
        #Add significant fold change information
        finalResult<-merge(finalResult, FCsum, by.x="row.names", by.y="row.names", all=TRUE)
        finalResult<-data.frame(finalResult[,-1], row.names=finalResult[,1])
        
        colnames(finalResult)<- c("Pathway_Protein_FDR","Pathway_Protein_p_value", "Sig_protein_Number", "Total_Protein_Number", "Protein_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
      
    }
    if(method[1] == "mean.significance"){
      if(calculateFoldChange == FALSE){
        
        #Calculate number of proteins mapped to pathways
        N = sum(row.names(proteinFilt)%in%background, na.rm=TRUE)
        m = data.frame(colSums(proteinFilt[background,], na.rm=TRUE))
        message(paste("Matched ", N, " proteins out of ", length(background), " total proteins to ", ncol(proteinFilt), " pathways out of ", ncol(proteins), " total pathways with a count threshold of ", countThreshold, sep=""))
        message("Calculating mean pathway significance values.  This method is slower than the others.  To speed up, specify pathway type, increase count threshold, or change method.")
        
        #Calculate mean signifcance values for each pathway
        meanFrame = data.frame()
        for(i in 1:ncol(proteinFilt)){
          filter = row.names(proteinFilt[which(proteinFilt[,i] == 1),])
          meanFrame[i,1] = mean(testResult[,1][which(row.names(testResult)%in%filter)], na.rm=TRUE)
        }
        row.names(meanFrame)<-colnames(proteinFilt)
        
        #Compile results into data frame
        result = data.frame()
        for (i in 1:length(meanFrame[,1])){
          result[i,1] = meanFrame[i,1]
          result[i,2] = m[i,1]
        }
        
        row.names(result)=row.names(m)
        colnames(result)=c("Mean_Pathway_Significance", "Sig_protein_Number")
        result<-result[order(result[,1]),]
        
        return(result)
        
      }else{
        
        #Calculate number of genes mapped to pathways
        N = sum(row.names(proteinFilt)%in%background, na.rm=TRUE)
        m = data.frame(colSums(proteinFilt[background,], na.rm=TRUE))
        message(paste("Matched ", N, " proteins out of ", length(background), " total proteins to ", ncol(proteinFilt), " pathways out of ", ncol(proteins), " total pathways with a count threshold of ", countThreshold, sep=""))
        message("Calculating mean pathway significance values.  This method is slower than the others.  To speed up, specify pathway type, increase count threshold, or change method.")
        
        #Calculate mean signifcance values for each pathway
        meanFrame = data.frame()
        for(i in 1:ncol(proteinFilt)){
          filter = row.names(proteinFilt[which(proteinFilt[,i] == 1),])
          meanFrame[i,1] = mean(testResult[,1][which(row.names(testResult)%in%filter)], na.rm=TRUE)
        }
        row.names(meanFrame)<-colnames(proteinFilt)
        
        #Calculate pathway fold changes
        FCsum = data.frame()
        for(i in 1:ncol(proteinFilt)){
          filter = row.names(proteinFilt[which(proteinFilt[,i] == 1),])
          FCsum[i,1] = sum(testResult[,2][which(row.names(testResult)%in%filter)], na.rm=TRUE)
        }
        row.names(FCsum)<-colnames(proteinFilt)
        
        
        #Compile results into data frame
        result = data.frame()
        for (i in 1:length(meanFrame[,1])){
          result[i,1] = meanFrame[i,1]
          result[i,2] = m[i,1]
          result[i,3] = FCsum[i,1]
        }
        
        row.names(result)=row.names(m)
        colnames(result)=c("Mean_Pathway_Significance", "Sig_protein_Number", "Protein_Cumulative_Fold_Change")
        result<-result[order(result[,1]),]
        
        return(result)
      }
    }
  }
  
  metab.test <- function(testResult, background = row.names(testResult), metabolites, cutoff, method, calculateFoldChange = TRUE){
    metabFilt<-metabolites[,-which(colSums(metabolites[background, ], na.rm=TRUE)<countThreshold)]
    if(method[1] == "fisher.exact"){
      if(calculateFoldChange == FALSE) {
        
        toORA <- testResult[order(testResult[,1]),]
        
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        
        sigFilter <- row.names(metabFilt)%in%sigList
        
        k = data.frame(colSums(metabFilt[sigFilter,], na.rm=TRUE))
        
        backgroundFilter <- row.names(metabFilt)%in%background
        
        m = data.frame(colSums(metabFilt[backgroundFilter,], na.rm=TRUE))
        
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " metabolites out of ", length(background), " total metabolites to ", ncol(metabFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        n = sum(sigFilter, na.rm = TRUE)
        #Calculate Fold Enrichment
        FE<-data.frame((k/n)/(m/N))
        result = data.frame()
        for (i in 1:length(k[,1])){
          result[i,1] = fisher.test(as.matrix(data.frame(x = c(k[i,], (m[i,]-k[i,])), y = c((n-k[i,]), (N+k[i,]-n-m[i,])))))$p.value
          result[i,2] = k[i,1]
          result[i,3] = m[i,1]
          result[i,4] = FE[i,1]
        }
        
        
        toFDR = data.frame(replace(result[,1], result[,1]>1, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        
        finalResult = cbind(FDR, result)
        
        
        row.names(finalResult)<-row.names(k)
        colnames(finalResult)<- c("Pathway_Metabolite_FDR", "Pathway_Metabolite_p_value", "Sig_Metabolite_Number", "Total_Metabolite_Number", "Metabolite_Fold_Enrichment")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
      }
      
      else {
        
        toORA <- testResult[order(testResult[,1]),]
        
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        
        sigFilter <- row.names(metabFilt)%in%sigList
        
        
        #Calculate cumulative Fold Change for pathway directionality
        toFCsum = data.frame(metabFilt[sigFilter,])
        toFCsum<-toFCsum[,-which(colSums(toFCsum, na.rm=TRUE)<countThreshold)]
        FCsum = data.frame()
        for(i in 1:ncol(toFCsum)){
          filter = row.names(toFCsum[which(toFCsum[,i] == 1),])
          FCsum[i,1] = sum(toORA[,2][which(row.names(toORA)%in%filter)])
        }
        row.names(FCsum)<-colnames(toFCsum)
        
        k = data.frame(colSums(metabFilt[sigFilter,], na.rm=TRUE))
        
        backgroundFilter <- row.names(metabFilt)%in%background
        
        m = data.frame(colSums(metabFilt[backgroundFilter,], na.rm=TRUE))
        
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " metabolites out of ", length(background), " total metabolites to ", ncol(metabFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        n = sum(sigFilter, na.rm=TRUE)
        #Calculate Fold Enrichment
        FE<-data.frame((k/n)/(m/N))
        result = data.frame()
        for (i in 1:length(k[,1])){      
          result[i,1] = fisher.test(as.matrix(data.frame(x = c(k[i,], (m[i,]-k[i,])), y = c((n-k[i,]), (N+k[i,]-n-m[i,])))))$p.value
          result[i,2] = k[i,1]
          result[i,3] = m[i,1]
          result[i,4] = FE[i,1]
        }
        
        
        toFDR = data.frame(replace(result[,1], result[,1]>1, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        
        finalResult = cbind(FDR, result)        
        row.names(finalResult)<-row.names(k)
        
        #Add significant fold change information
        finalResult<-merge(finalResult, FCsum, by.x="row.names", by.y="row.names", all=TRUE)
        finalResult<-data.frame(finalResult[,-1], row.names=finalResult[,1])
        
        colnames(finalResult)<- c("Pathway_Metabolite_FDR","Pathway_Metabolite_p_value", "Sig_Metabolite_Number", "Total_Metabolite_Number", "Metabolite_Fold_Enrichment", "Metabolite_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
      }
    }
    if(method[1] == "EASE"){
      if(calculateFoldChange == FALSE) {
        
        toORA <- testResult[order(testResult[,1]),]
        
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        
        sigFilter <- row.names(metabFilt)%in%sigList
        
        k = data.frame(colSums(metabFilt[sigFilter,], na.rm=TRUE))
        
        backgroundFilter <- row.names(metabFilt)%in%background
        
        m = data.frame(colSums(metabFilt[backgroundFilter,], na.rm=TRUE))
        
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " metabolites out of ", length(background), " total metabolites to ", ncol(metabFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        n = sum(sigFilter, na.rm = TRUE)
        #Calculate Fold Enrichment
        FE<-data.frame((k/n)/(m/N))
        result = data.frame()
        for (i in 1:length(k[,1])){
          result[i,1] = fisher.test(as.matrix(data.frame(x = c(max(k[i,]-1,0), (m[i,]-k[i,])), y = c((n-k[i,]), (N+k[i,]-n-m[i,])))))$p.value
          result[i,2] = k[i,1]
          result[i,3] = m[i,1]
          result[i,4] = FE[i,1]
        }
        
        
        toFDR = data.frame(replace(result[,1], result[,1]>1, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        
        finalResult = cbind(FDR, result)
        
        
        row.names(finalResult)<-row.names(k)
        colnames(finalResult)<- c("Pathway_Metabolite_FDR", "Pathway_Metabolite_p_value", "Sig_Metabolite_Number", "Total_Metabolite_Number", "Metabolite_Fold_Enrichment")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
      }
      
      else {
        
        toORA <- testResult[order(testResult[,1]),]
        
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        
        sigFilter <- row.names(metabFilt)%in%sigList
        
        
        toFCsum = data.frame(metabFilt[sigFilter,])
        FCsum = data.frame()
        for(i in 1:ncol(toFCsum)){
          filter = row.names(toFCsum[which(toFCsum[,i] == 1),])
          FCsum[i,1] = sum(toORA[,2][which(row.names(toORA)%in%filter)])
        }
        
        
        k = data.frame(colSums(metabFilt[sigFilter,], na.rm=TRUE))
        
        backgroundFilter <- row.names(metabFilt)%in%background
        
        m = data.frame(colSums(metabFilt[backgroundFilter,], na.rm=TRUE))
        
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " metabolites out of ", length(background), " total metabolites to ", ncol(metabFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        n = sum(sigFilter, na.rm=TRUE)
        #Calculate Fold Enrichment
        FE<-data.frame((k/n)/(m/N))
        result = data.frame()
        for (i in 1:length(k[,1])){      
          result[i,1] = fisher.test(as.matrix(data.frame(x = c(max(k[i,]-1,0), (m[i,]-k[i,])), y = c((n-k[i,]), (N+k[i,]-n-m[i,])))))$p.value
          result[i,3] = k[i,1]
          result[i,4] = m[i,1]
          result[i,5] = FE[i,1]
        }
        
        
        toFDR = data.frame(replace(result[,1], result[,1]>1, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        
        finalResult = cbind(FDR, result)        
        row.names(finalResult)<-row.names(k)
        
        #Add significant fold change information
        finalResult<-merge(finalResult, FCsum, by.x="row.names", by.y="row.names", all=TRUE)
        finalResult<-data.frame(finalResult[,-1], row.names=finalResult[,1])
        
        colnames(finalResult)<- c("Pathway_Metabolite_FDR","Pathway_Metabolite_p_value", "Sig_Metabolite_Number", "Total_Metabolite_Number", "Metabolite_Fold_Enrichment", "Metabolite_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
      }
    }
    if(method[1] == "hypergeometric"){
      if(calculateFoldChange == FALSE){
        
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(metabFilt)%in%sigList
        
        #Determine number of significant markers overlapping with each gene pathway
        q = data.frame(colSums(metabFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified gene pathways
        backgroundFilter <- row.names(metabFilt)%in%background
        #Determine number of background markers overlapping with each gene pathway
        m = data.frame(colSums(metabFilt[backgroundFilter,], na.rm=TRUE))
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " metabolites out of ", length(background), " total metabolites to ", ncol(metabFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        #Perform Hypergeometric Test
        k = sum(sigFilter, na.rm=TRUE)
        result = data.frame()
        for (i in 1:length(q[,1])){
          result[i,1] = phyper(q = q[i,],  m=m[i,], k = k, n = sum(backgroundFilter, na.rm=TRUE)-m[i,], lower.tail = FALSE)
          result[i,2] = q[i,1]
          result[i,3] = m[i,1]
        }
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1 | result[,1]==0, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(q)
        colnames(finalResult)<- c("Pathway_Metabolite_FDR","Pathway_Metabolite_p_value", "Sig_Metabolite_Number", "Total_Metabolite_Number")
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
      else{
        
        #Generate dataframe for ORA analysis
        toORA <- testResult[order(testResult[,1]),]
        #Filter for only markers below the specified FDR cutoff
        sigList <- row.names(subset(toORA, toORA[,1] < cutoff))
        #Determine which significant markers are found in specified gene pathways
        sigFilter <- row.names(metabFilt)%in%sigList
        
        #Calculate cumulative Fold Change for pathway directionality
        toFCsum = data.frame(metabFilt[sigFilter,])
        toFCsum<-toFCsum[,-which(colSums(toFCsum, na.rm=TRUE)<countThreshold)]
        FCsum = data.frame()
        for(i in 1:ncol(toFCsum)){
          filter = row.names(toFCsum[which(toFCsum[,i] == 1),])
          FCsum[i,1] = sum(toORA[,2][which(row.names(toORA)%in%filter)])
        }
        row.names(FCsum)<-colnames(toFCsum)
        
        #Determine number of significant markers overlapping with each gene pathway
        q = data.frame(colSums(metabFilt[sigFilter,], na.rm=TRUE))
        #Determine which background markers are found in specified gene pathways
        backgroundFilter <- row.names(metabFilt)%in%background
        #Determine number of background markers overlapping with each gene pathway
        m = data.frame(colSums(metabFilt[backgroundFilter,], na.rm=TRUE))
        N = sum(backgroundFilter, na.rm=TRUE)
        message(paste("Matched ", N, " metabolites out of ", length(background), " total metabolites to ", ncol(metabFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        #Perform Hypergeometric Test
        k = sum(sigFilter, na.rm=TRUE)
        result = data.frame()
        for (i in 1:length(q[,1])){
          result[i,1] = phyper(q = q[i,],  m=m[i,], k = k, n = sum(backgroundFilter, na.rm=TRUE)-m[i,], lower.tail = FALSE)
          result[i,2] = q[i,1]
          result[i,3] = m[i,1]
        }
        #Perform tail-based (BH) Fdr analysis
        toFDR = data.frame(replace(result[,1], result[,1]>1 | result[,1]==0, 1))
        FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
        #Create final dataframe
        finalResult = cbind(FDR, result)
        
        #Add appropriate labes to row and column names
        row.names(finalResult)<-row.names(q)
        
        
        #Add significant fold change information
        finalResult<-merge(finalResult, FCsum, by.x="row.names", by.y="row.names", all=TRUE)
        finalResult<-data.frame(finalResult[,-1], row.names=finalResult[,1])
        
        colnames(finalResult)<- c("Pathway_Metabolite_FDR","Pathway_Metabolite_p_value", "Sig_Metabolite_Number", "Total_Metabolite_Number", "Metabolite_Cumulative_Fold_Change")
        
        finalResult<-finalResult[order(finalResult[,2]),]
        return(finalResult)
        
      }
      
    }
    if(method[1] == "mean.significance"){
      if(calculateFoldChange == FALSE){
        
        #Calculate number of metabolites mapped to pathways
        N = sum(row.names(metabFilt)%in%background, na.rm=TRUE)
        m = data.frame(colSums(metabFilt[background,], na.rm=TRUE))
        message(paste("Matched ", N, " metabolites out of ", length(background), " total metabolites to ", ncol(metabFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        message("Calculating mean pathway significance values.  This method is slower than the others.  To speed up, specify pathway type, increase count threshold, or change method.")
        
        #Calculate mean signifcance values for each pathway
        meanFrame = data.frame()
        for(i in 1:ncol(metabFilt)){
          filter = row.names(metabFilt[which(metabFilt[,i] == 1),])
          meanFrame[i,1] = mean(testResult[,1][which(row.names(testResult)%in%filter)], na.rm=TRUE)
        }
        row.names(meanFrame)<-colnames(metabFilt)
        
        #Compile results into data frame
        result = data.frame()
        for (i in 1:length(meanFrame[,1])){
          result[i,1] = meanFrame[i,1]
          result[i,2] = m[i,1]
        }
        
        row.names(result)=row.names(m)
        colnames(result)=c("Mean_Pathway_Significance", "Sig_Metabolite_Number")
        result<-result[order(result[,1]),]
        
        return(result)
        
      }else{
        
        #Calculate number of metabolites mapped to pathways
        N = sum(row.names(metabFilt)%in%background, na.rm=TRUE)
        m = data.frame(colSums(metabFilt[background,], na.rm=TRUE))
        message(paste("Matched ", N, " metabolites out of ", length(background), " total metabolites to ", ncol(metabFilt), " pathways out of ", ncol(metabolites), " total pathways with a count threshold of ", countThreshold, sep=""))
        message("Calculating mean pathway significance values.  This method is slower than the others.  To speed up, specify pathway type, increase count threshold, or change method.")
        
        #Calculate mean signifcance values for each pathway
        meanFrame = data.frame()
        for(i in 1:ncol(metabFilt)){
          filter = row.names(metabFilt[which(metabFilt[,i] == 1),])
          meanFrame[i,1] = mean(testResult[,1][which(row.names(testResult)%in%filter)], na.rm=TRUE)
        }
        row.names(meanFrame)<-colnames(metabFilt)
        
        #Calculate pathway fold changes
        FCsum = data.frame()
        for(i in 1:ncol(metabFilt)){
          filter = row.names(metabFilt[which(metabFilt[,i] == 1),])
          FCsum[i,1] = sum(testResult[,2][which(row.names(testResult)%in%filter)], na.rm=TRUE)
        }
        row.names(FCsum)<-colnames(metabFilt)
        
        
        #Compile results into data frame
        result = data.frame()
        for (i in 1:length(meanFrame[,1])){
          result[i,1] = meanFrame[i,1]
          result[i,2] = m[i,1]
          result[i,3] = FCsum[i,1]
        }
        
        row.names(result)=row.names(m)
        colnames(result)=c("Mean_Pathway_Significance", "Sig_Metabolite_Number", "Metabolite_Cumulative_Fold_Change")
        result<-result[order(result[,1]),]
        
        return(result)
      }
    }
  }
  
  combined.test <- function(gene_result=NULL, protein_result=NULL, metab_result=NULL, calculateFoldChange = TRUE, method){
    #define fisher function
    Fisher.test <- function(p) {
      Xsq <- -2*sum(log(p))
      p.val <- pchisq(Xsq, df = 2*length(p),lower.tail=FALSE)
      return(p.value = p.val)
    }
    
    if(method[1]=="mean.significance"){
      if(is.null(protein_result)){
        if(calculateFoldChange == FALSE){  
          #filter gene_result for only pathways which include metabolites
          genes_w_metab<-row.names(gene_result)[which(row.names(gene_result)%in%row.names(metab_result))]
          gene_result<-gene_result[genes_w_metab, ]
          metab_result<-metab_result[genes_w_metab, ]
          
          result = data.frame()
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
          
          result = data.frame()
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
      if(is.null(gene_result)){
        if(calculateFoldChange == FALSE){  
          #filter protein_result for only pathways which include metabolites
          proteins_w_metab<-row.names(protein_result)[which(row.names(protein_result)%in%row.names(metab_result))]
          protein_result<-protein_result[proteins_w_metab, ]
          metab_result<-metab_result[proteins_w_metab, ]
          
          result = data.frame()
          for(i in 1:nrow(metab_result)){
            result[i,1] = sum(protein_result[i,1]*(protein_result[i,2]/(protein_result[i,2]+metab_result[i,2])) + metab_result[i,1]*(metab_result[i,2]/(metab_result[i,2]+protein_result[i,2])))
            result[i,2] = protein_result[i,2]
            result[i,3] = metab_result[i,2]
          }
          row.names(result)<-row.names(protein_result)
          colnames(result)<- c("Mean_Combined_Pathway_Significance", "Total_Protein_Number", "Total_Metabolite_Number")
          
          result<-result[order(result[,1]),]
        }
        else{  
          #filter protein_result for only pathways which include metabolites
          proteins_w_metab<-row.names(protein_result)[which(row.names(protein_result)%in%row.names(metab_result))]
          protein_result<-protein_result[proteins_w_metab, ]
          metab_result<-metab_result[proteins_w_metab, ]
          
          result = data.frame()
          for(i in 1:nrow(metab_result)){
            result[i,1] = sum(protein_result[i,1]*(protein_result[i,2]/(protein_result[i,2]+metab_result[i,2])) + metab_result[i,1]*(metab_result[i,2]/(metab_result[i,2]+protein_result[i,2])))
            result[i,2] = protein_result[i,2]
            result[i,3] = protein_result[i,3]
            result[i,4] = metab_result[i,2]
            result[i,5] = metab_result[i,3]
          }
          row.names(result)<-row.names(protein_result)
          colnames(result)<- c("Mean_Combined_Pathway_Significance", "Total_Protein_Number", "Protein_Cumulative_Fold_Change","Total_Metabolite_Number", "Metabolite_Cumulative_Fold_Change")
          
          result<-result[order(result[,1]),]
        }
      }
      if(is.null(metab_result)){
        if(calculateFoldChange == FALSE){  
          #filter protein_result for only pathways which include genes
          proteins_w_genes<-row.names(protein_result)[which(row.names(protein_result)%in%row.names(gene_result))]
          protein_result<-protein_result[proteins_w_genes, ]
          gene_result<-gene_result[proteins_w_genes, ]
          
          result = data.frame()
          for(i in 1:nrow(gene_result)){
            result[i,1] = sum(protein_result[i,1]*(protein_result[i,2]/(protein_result[i,2]+gene_result[i,2])) + gene_result[i,1]*(gene_result[i,2]/(gene_result[i,2]+protein_result[i,2])))
            result[i,2] = protein_result[i,2]
            result[i,3] = gene_result[i,2]
          }
          row.names(result)<-row.names(protein_result)
          colnames(result)<- c("Mean_Combined_Pathway_Significance", "Total_Protein_Number", "Total_Gene_Number")
          
          result<-result[order(result[,1]),]
        }
        else{  
          #filter protein_result for only pathways which include metabolites
          proteins_w_genes<-row.names(protein_result)[which(row.names(protein_result)%in%row.names(gene_result))]
          protein_result<-protein_result[proteins_w_genes, ]
          gene_result<-metab_result[proteins_w_genes, ]
          
          result = data.frame()
          for(i in 1:nrow(gene_result)){
            result[i,1] = sum(protein_result[i,1]*(protein_result[i,2]/(protein_result[i,2]+gene_result[i,2])) + gene_result[i,1]*(gene_result[i,2]/(gene_result[i,2]+protein_result[i,2])))
            result[i,2] = protein_result[i,2]
            result[i,3] = protein_result[i,3]
            result[i,4] = gene_result[i,2]
            result[i,5] = gene_result[i,3]
          }
          row.names(result)<-row.names(protein_result)
          colnames(result)<- c("Mean_Combined_Pathway_Significance", "Total_Protein_Number", "Protein_Cumulative_Fold_Change","Total_Gene_Number", "Gene_Cumulative_Fold_Change")
          
          result<-result[order(result[,1]),]
        }
      }
      if(!is.null(gene_result)&!is.null(protein_result)&!is.null(metab_result)){
        if(calculateFoldChange == FALSE){  
          #filter for common pathways
          proteins_w_genes = row.names(protein_result)[which(row.names(protein_result)%in%row.names(gene_result))]
          proteins_w_genes_metab = proteins_w_genes[which(proteins_w_genes%in%row.names(metab_result))]
          protein_result<-protein_result[proteins_w_genes_metab, ]
          gene_result<-gene_result[proteins_w_genes_metab, ]
          metab_result<-metab_result[proteins_w_genes_metab, ]
          
          result = data.frame()
          for(i in 1:nrow(gene_result)){
            result[i,1] = sum(protein_result[i,1]*(protein_result[i,2]/(protein_result[i,2]+gene_result[i,2]+metab_result[i,2])) + gene_result[i,1]*(gene_result[i,2]/(gene_result[i,2]+protein_result[i,2]+metab_result[i,2])) + metab_result[i,2]*(metab_result[i,2]/(gene_result[i,2]+protein_result[i,2]+metab_result[i,2])))
            result[i,2] = protein_result[i,2]
            result[i,3] = gene_result[i,2]
            result[i,4] = metab_result[i,2]
          }
          row.names(result)<-row.names(protein_result)
          colnames(result)<- c("Mean_Combined_Pathway_Significance", "Total_Protein_Number", "Total_Gene_Number", "Total_Metabolite_Number")
          
          result<-result[order(result[,1]),]
        }
        else{  
          #filter for common pathways
          proteins_w_genes = row.names(protein_result)[which(row.names(protein_result)%in%row.names(gene_result))]
          proteins_w_genes_metab = proteins_w_genes[which(proteins_w_genes%in%row.names(metab_result))]
          protein_result<-protein_result[proteins_w_genes_metab, ]
          gene_result<-gene_result[proteins_w_genes_metab, ]
          metab_result<-metab_result[proteins_w_genes_metab, ]
          
          result = data.frame()
          for(i in 1:nrow(gene_result)){
            result[i,1] = sum(protein_result[i,1]*(protein_result[i,2]/(protein_result[i,2]+gene_result[i,2]+metab_result[i,2])) + gene_result[i,1]*(gene_result[i,2]/(gene_result[i,2]+protein_result[i,2]+metab_result[i,2])) + metab_result[i,2]*(metab_result[i,2]/(gene_result[i,2]+protein_result[i,2]+metab_result[i,2])))
            result[i,2] = protein_result[i,2]
            result[i,3] = protein_result[i,3]
            result[i,4] = gene_result[i,2]
            result[i,5] = gene_result[i,3]
            result[i,6] = metab_result[i,2]
            result[i,7] = metab_result[i,3]
          }
          row.names(result)<-row.names(protein_result)
          colnames(result)<- c("Mean_Combined_Pathway_Significance", "Total_Protein_Number", "Protein_Cumulative_Fold_Change","Total_Gene_Number", "Gene_Cumulative_Fold_Change", "Total_Metabolite_Number", "Metabolite_Cumulative_Fold_Change")
          
          result<-result[order(result[,1]),]
        }
      }
    }
    else{
      
      if(is.null(protein_result)){
        if(calculateFoldChange == FALSE){  
          #filter gene_result for only pathways which include metabolites
          genes_w_metab<-row.names(gene_result)[which(row.names(gene_result)%in%row.names(metab_result))]
          gene_result<-gene_result[genes_w_metab, ]
          metab_result<-metab_result[genes_w_metab, ]
          
          #perform fisher method for joint p-values on all pathways
          result = data.frame()
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
          result = data.frame()
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
      if(is.null(gene_result)){
        if(calculateFoldChange == FALSE){  
          #filter protein_result for only pathways which include metabolites
          proteins_w_metab = row.names(protein_result)[which(row.names(protein_result)%in%row.names(metab_result))]
          protein_result<-protein_result[proteins_w_metab, ]
          metab_result<-metab_result[proteins_w_metab, ]
          
          #perform fisher method for joint p-values on all pathways
          result = data.frame()
          for(i in 1:length(metab_result[,1])){
            result[i,1] = Fisher.test(p = c(protein_result[i,2], metab_result[i,2]))
            result[i,2] = protein_result[i,2]
            result[i,3] = protein_result[i,3]
            result[i,4] = protein_result[i,4]
            result[i,5] = protein_result[i,5]
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
          colnames(finalResult) <- c("Pathway_Combined_FDR","Pathway_Combined_pValue",colnames(protein_result[-1]), colnames(metab_result[-1]))
        
          finalResult<-finalResult[order(finalResult[,2]),]
          return(finalResult)
          
        }
        else{
          #filter protein_result for only pathways which include metabolites
          proteins_w_metab = row.names(protein_result)[which(row.names(protein_result)%in%row.names(metab_result))]
          protein_result<-protein_result[proteins_w_metab, ]
          metab_result<-metab_result[proteins_w_metab, ]
          
          #perform fisher method for joint p-values on all pathways
          result = data.frame()
          for(i in 1:length(metab_result[,1])){
            result[i,1] = Fisher.test(p = c(protein_result[i,2], metab_result[i,2]))
            result[i,2] = protein_result[i,2]
            result[i,3] = protein_result[i,3]
            result[i,4] = protein_result[i,4]
            result[i,5] = protein_result[i,5]
            result[i,6] = protein_result[i,6]
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
          colnames(finalResult) <- c("Pathway_Combined_FDR","Pathway_Combined_pValue",colnames(protein_result[-1]), colnames(metab_result[-1]))
          
          finalResult<-finalResult[order(finalResult[,2]),]
          return(finalResult)
          
        }
      }
      if(is.null(metab_result)){
        if(calculateFoldChange == FALSE){  
          proteins_w_genes = row.names(protein_result)[which(row.names(protein_result)%in%row.names(gene_result))]
          protein_result<-protein_result[proteins_w_genes, ]
          gene_result<-gene_result[proteins_w_genes, ]
          #perform fisher method for joint p-values on all pathways
          result = data.frame()
          for(i in 1:length(gene_result[,1])){
            result[i,1] = Fisher.test(p = c(gene_result[i,2], protein_result[i,2]))
            result[i,2] = gene_result[i,2]
            result[i,3] = gene_result[i,3]
            result[i,4] = gene_result[i,4]
            result[i,5] = gene_result[i,5]
            result[i,6] = protein_result[i,2]
            result[i,7] = protein_result[i,3]
            result[i,8] = protein_result[i,4]
            result[i,9] = protein_result[i,5]
            
          }
          
          #Perform tail-based (BH) Fdr analysis
          toFDR = data.frame(replace(result[,1], result[,1]>1 | result[,1]==0, 1))
          FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
          #Create final dataframe
          finalResult = cbind(FDR, result)
          
          #label row and column names appropriately
          row.names(finalResult)<-row.names(metab_result)
          colnames(finalResult) <- c("Pathway_Combined_FDR","Pathway_Combined_pValue",colnames(gene_result[-1]), colnames(protein_result[-1]))
          
          finalResult<-finalResult[order(finalResult[,2]),]
          return(finalResult)
          
        }
        else{
          proteins_w_genes = row.names(protein_result)[which(row.names(protein_result)%in%row.names(gene_result))]
          protein_result<-protein_result[proteins_w_genes, ]
          gene_result<-gene_result[proteins_w_genes, ]
          #perform fisher method for joint p-values on all pathways
          result = data.frame()
          for(i in 1:length(gene_result[,1])){
            result[i,1] = Fisher.test(p = c(gene_result[i,2], protein_result[i,2]))
            result[i,2] = gene_result[i,2]
            result[i,3] = gene_result[i,3]
            result[i,4] = gene_result[i,4]
            result[i,5] = gene_result[i,5]
            result[i,6] = gene_result[i,6]
            result[i,7] = protein_result[i,2]
            result[i,8] = protein_result[i,3]
            result[i,9] = protein_result[i,4]
            result[i,10] = protein_result[i,5]
            result[i,11] = protein_result[i,6]
            
          }
          
          #Perform tail-based (BH) Fdr analysis
          toFDR = data.frame(replace(result[,1], result[,1]>1 | result[,1]==0, 1))
          FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
          #Create final dataframe
          finalResult = cbind(FDR, result)
          
          #label row and column names appropriately
          row.names(finalResult)<-row.names(metab_result)
          colnames(finalResult) <- c("Pathway_Combined_FDR","Pathway_Combined_pValue",colnames(gene_result[-1]), colnames(protein_result[-1]))
          
          finalResult<-finalResult[order(finalResult[,2]),]
          return(finalResult)
          
        }
      }
      if(!is.null(gene_result)&!is.null(protein_result)&!is.null(metab_result)){
        if(calculateFoldChange == FALSE){  
          #filter for common pathways
          proteins_w_genes = row.names(protein_result)[which(row.names(protein_result)%in%row.names(gene_result))]
          proteins_w_genes_metab = proteins_w_genes[which(proteins_w_genes%in%row.names(metab_result))]
          protein_result<-protein_result[proteins_w_genes_metab, ]
          gene_result<-gene_result[proteins_w_genes_metab, ]
          metab_result<-metab_result[proteins_w_genes_metab, ]
          
          
          #perform fisher method for joint p-values on all pathways
          result = data.frame()
          for(i in 1:length(metab_result[,1])){
            result[i,1] = Fisher.test(p = c(gene_result[i,2], protein_result[i,2], metab_result[i,2]))
            result[i,2] = gene_result[i,2]
            result[i,3] = gene_result[i,3]
            result[i,4] = gene_result[i,4]
            result[i,5] = gene_result[i,5]
            result[i,6] = protein_result[i,2]
            result[i,7] = protein_result[i,3]
            result[i,8] = protein_result[i,4]
            result[i,9] = protein_result[i,5]
            result[i,10] = metab_result[i,2]
            result[i,11] = metab_result[i,3]
            result[i,12] = metab_result[i,4]
            result[i,13] = metab_result[i,5]
            
          }
          
          #Perform tail-based (BH) Fdr analysis
          toFDR = data.frame(replace(result[,1], result[,1]>1 | result[,1]==0, 1))
          FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
          #Create final dataframe
          finalResult = cbind(FDR, result)
          
          #label row and column names appropriately
          row.names(finalResult)<-row.names(metab_result)
          colnames(finalResult) <- c("Pathway_Combined_FDR","Pathway_Combined_pValue", colnames(gene_result[-1]), colnames(protein_result[-1]), colnames(metab_result[-1]))
          
          finalResult<-finalResult[order(finalResult[,2]),]
          return(finalResult)
          
        }
        else{
          #filter for common pathways
          proteins_w_genes = row.names(protein_result)[which(row.names(protein_result)%in%row.names(gene_result))]
          proteins_w_genes_metab = proteins_w_genes[which(proteins_w_genes%in%row.names(metab_result))]
          protein_result<-protein_result[proteins_w_genes_metab, ]
          gene_result<-gene_result[proteins_w_genes_metab, ]
          metab_result<-metab_result[proteins_w_genes_metab, ]
          
          
          #perform fisher method for joint p-values on all pathways
          result = data.frame()
          for(i in 1:length(metab_result[,1])){
            result[i,1] = Fisher.test(p = c(gene_result[i,2], protein_result[i,2], metab_result[i,2]))
            result[i,2] = gene_result[i,2]
            result[i,3] = gene_result[i,3]
            result[i,4] = gene_result[i,4]
            result[i,5] = gene_result[i,5]
            result[i,6] = gene_result[i,6]
            result[i,7] = protein_result[i,2]
            result[i,8] = protein_result[i,3]
            result[i,9] = protein_result[i,4]
            result[i,10] = protein_result[i,5]
            result[i,11] = protein_result[i,6]
            result[i,12] = metab_result[i,2]
            result[i,13] = metab_result[i,3]
            result[i,14] = metab_result[i,4]
            result[i,15] = metab_result[i,5]
            result[i,16] = metab_result[i,6]
            
          }
          
          #Perform tail-based (BH) Fdr analysis
          toFDR = data.frame(replace(result[,1], result[,1]>1 | result[,1]==0, 1))
          FDR = data.frame(p.adjust(toFDR[,1], method=p.adjustMethod))
          #Create final dataframe
          finalResult = cbind(FDR, result)
          
          #label row and column names appropriately
          row.names(finalResult)<-row.names(metab_result)
          colnames(finalResult) <- c("Pathway_Combined_FDR","Pathway_Combined_pValue", colnames(gene_result[-1]), colnames(protein_result[-1]), colnames(metab_result[-1]))
          
          finalResult<-finalResult[order(finalResult[,2]),]
          return(finalResult)
          
        }
      }
    }
  }
  
  
  if(!is.null(geneResult)){
    
    gene<-gene.test(testResult=geneResult, method=method, genes=genes, cutoff=geneCutoff, calculateFoldChange=calculateFoldChange)
    geneR<-gene[-which(gene[,"Sig_Gene_Number"]<countThreshold),]
    geneR<-merge(pathIndex, geneR, by.x="row.names", by.y="row.names")
    geneR<-data.frame(geneR[,-1], row.names=geneR[,1])
    geneR<-geneR[order(geneR[,3]),]
  }
  
  if(!is.null(proteinResult)){
    
    protein<-protein.test(testResult=proteinResult, method=method, proteins=proteins, cutoff=proteinCutoff, calculateFoldChange=calculateFoldChange)
    proteinR<-protein[-which(protein[,"Sig_protein_Number"]<countThreshold),]
    proteinR<-merge(pathIndex, proteinR, by.x="row.names", by.y="row.names")
    proteinR<-data.frame(proteinR[,-1], row.names=proteinR[,1])
    proteinR<-proteinR[order(proteinR[,3]),]
  }
  
  if(!is.null(metaboliteResult)){
    
    metab<-metab.test(testResult=metaboliteResult, method=method, metabolites=metabolites, cutoff=metaboliteCutoff, calculateFoldChange=calculateFoldChange)
    metabR<-metab[-which(metab[,"Sig_Metabolite_Number"]<countThreshold),]
    metabR<-merge(pathIndex, metabR, by.x="row.names", by.y="row.names")
    metabR<-data.frame(metabR[,-1], row.names=metabR[,1])
    metabR<-metabR[order(metabR[,3]),]
  }
  
  result<-list(geneR, proteinR, metabR)
  names(result)<-c("Gene_Analysis", "Protein_Analysis", "Metabolite_Analysis")
  
  if(combine==TRUE){
    combineR<-combined.test(gene_result=geneR[,-1], protein_result=proteinR[,-1], metab_result=metabR[,-1], calculateFoldChange=calculateFoldChange, method=method)
    combineR<-merge(pathIndex, combineR, by.x="row.names", by.y="row.names")
    combineR<-data.frame(combineR[,-1], row.names=combineR[,1])
    combineR<-combineR[order(combineR[,3]),]
    return(combineR)
  }
  if(combine==FALSE){
    
    return(result)
  }
}