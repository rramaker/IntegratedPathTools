#' Perform interactive pathway analysis on rna and metabolite datasets
#'
#' This takes dataframes with rna and metabolite info and performs interactive pathway analysis
#' @param geneResult Dataframe with transcript IDs as rownames, transcript signifigance values (between 0 and 1) as the first column, and optionally transcript fold change values as the second column. Defaults to NULL
#' @param metaboliteResult Dataframe with metabolite IDs as rownames, metabolite signifigance values (between 0 and 1) as the first column, and optionally metabolite fold change values as the second column. Defaults to NULL
#' @param method method Method of pathway analysis to perform. Options include "fisher.exact",  "EASE", "mean.significance", and "hypergeometric.test". Defaults to "fisher.exact.
#' @param pathwayType Pathway database to use for pathway analysis.  Options include "KEGG", "SMPDB", "Reactome", or "All".  Defaults to "All"
#' @param countThreshold Integer indicating the minimum number of pathway members to require for a valid result. Defaults to 2.
#' @param p.adjustMethod Character indicating the type of multiple hypothesis correction to perform.  See ?p.adjust() for details
#' @keywords IPA
#' @import manipulate
#' @import ggplot2
#' @export
#' @return Returns an interactive dataplot of top 5 significant pathways with ticker to adjust FDR thresholds
#' @examples 
#' data(kData)
#' data(rData)
#' ## Not run: manipulateIPA(kData, rData, method="fisher.exact") ## End(Not run)
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



manipulateIPA<-function(geneResult=NULL, metaboliteResult=NULL, method, countThreshold=2, p.adjustMethod="BH", pathwayType=c("All", "KEGG", "SMPDB", "Reactome")){
  
  data(genes, envir=environment())
  data(metabolites, envir=environment())
  
  
  if(pathwayType[1] == "All"){ 
    genes<-genes
    metabolites<-metabolites
  }
  if(pathwayType[1] == "KEGG"){
    genes<-genes[,colnames(genes)[which(grepl("hsa", colnames(genes)))]]
    metabolites<-metabolites[,colnames(metabolites)[which(grepl("hsa", colnames(metabolites)))]]
  }
  if(pathwayType[1] == "SMPDB"){
    genes<-genes[,colnames(genes)[which(grepl("SMP", colnames(genes)))]]
    metabolites<-metabolites[,colnames(metabolites)[which(grepl("SMP", colnames(metabolites)))]]
  }
  if(pathwayType[1] == "Reactome"){
    genes<-genes[,colnames(genes)[which(grepl("REACT", colnames(genes)))]]
    metabolites<-metabolites[,colnames(metabolites)[which(grepl("REACT", colnames(metabolites)))]]
  }
  
  gene.test <- function(testResult, genes, cutoff, method, calculateFoldChange=FALSE, countThreshold=2, p.adjustMethod="BH"){
      background = row.names(testResult)
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
  metab.test <- function(testResult, metabolites, cutoff, method, calculateFoldChange = FALSE, countThreshold=2, p.adjustMethod="BH"){
    background = row.names(testResult)
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
  combined.test <- function(gene_result=NULL, metab_result=NULL, calculateFoldChange = FALSE){
    #define fisher function
    Fisher.test <- function(p) {
      Xsq <- -2*sum(log(p))
      p.val <- pchisq(Xsq, df = 2*length(p),lower.tail=FALSE)
      return(p.value = p.val)
    }
    
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
      
    
  
  
  if(!is.null(geneResult)&is.null(metaboliteResult)){ 
    myBar <- function(geneCutoff){
    
      x<-gene.test(testResult=geneResult, genes=genes, cutoff=geneCutoff, method=method)
      toPlot<- x[1:5,]
      Annotation = factor(row.names(toPlot), levels = unique(row.names(toPlot)))
      p1 = ggplot2::qplot(Annotation, -log(toPlot[,1], base=2), geom= "bar", stat="identity", binwidth= 0.00001)
      p2 = p1 + ggplot2::ylab("Negative Log FDR")+ggplot2::xlab("")+ggplot2::coord_flip() + ggplot2::theme( panel.background = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black", size = 1), axis.text = ggplot2::element_text(face = "bold", colour = "black"),  axis.ticks = ggplot2::element_line(colour = "black"), axis.title = ggplot2::element_text(face="bold", colour = "black"))
      return(p2)
    }
    return(manipulate::manipulate(myBar(geneCutoff), geneCutoff=manipulate::slider((min(geneResult[,1], na.rm=T)+0.05),0.5)))
  }
  if(is.null(geneResult)&!is.null(metaboliteResult)){ 
    myBar <- function(metaboliteCutoff){
      
      x<-metab.test(testResult=metaboliteResult, metabolites=metabolites, cutoff=metaboliteCutoff, method=method)
      toPlot<- x[1:5,]
      Annotation = factor(row.names(toPlot), levels = unique(row.names(toPlot)))
      p1 = ggplot2::qplot(Annotation, -log(toPlot[,1], base=2), geom= "bar", stat="identity", binwidth= 0.00001)
      p2 = p1 + ggplot2::ylab("Negative Log FDR")+ggplot2::xlab("")+ggplot2::coord_flip() + ggplot2::theme( panel.background = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black", size = 1), axis.text = ggplot2::element_text(face = "bold", colour = "black"),  axis.ticks = ggplot2::element_line(colour = "black"), axis.title = ggplot2::element_text(face="bold", colour = "black"))
      return(p2)
    }
    return(manipulate::manipulate(myBar(metaboliteCutoff), metaboliteCutoff=manipulate::slider((min(metaboliteResult[,1], na.rm=T)+0.05),0.5))) 
  }
  if(!is.null(geneResult)&!is.null(metaboliteResult)){
    myBar <- function(geneCutoff, metaboliteCutoff){
      geneR<-gene.test(testResult=geneResult, genes=genes, cutoff=geneCutoff, method=method)
      metabR<-metab.test(testResult=metaboliteResult, metabolites=metabolites, cutoff=metaboliteCutoff, method=method)
      x<-combined.test(gene_result=geneR, metab_result=metabR)
      toPlot<- x[1:5,]
      Annotation = factor(row.names(toPlot), levels = unique(row.names(toPlot)))
      p1 = ggplot2::qplot(Annotation, -log(toPlot[,1], base=2), geom= "bar", stat="identity", binwidth= 0.00001)
      p2 = p1 + ggplot2::ylab("Negative Log FDR")+ggplot2::xlab("")+ggplot2::coord_flip() + ggplot2::theme( panel.background = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black", size = 1), axis.text = ggplot2::element_text(face = "bold", colour = "black"),  axis.ticks = ggplot2::element_line(colour = "black"), axis.title = ggplot2::element_text(face="bold", colour = "black"))
      return(p2)
    }
    return(manipulate::manipulate(myBar(geneCutoff, metaboliteCutoff), geneCutoff=manipulate::slider((min(geneResult[,1], na.rm=T)+0.05),0.5), metaboliteCutoff=manipulate::slider((min(metaboliteResult[,1], na.rm=T)+0.05),0.5)))
  }

}




