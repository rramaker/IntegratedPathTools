#' Perform interactive post hoc bootstrap analysis on pathway over-resentation
#'
#' This takes dataframes with rna and metabolite info and performs interactive simulated bootstrap analysis of pathway over-representation.
#' @param geneResult Dataframe with transcript IDs as rownames, transcript signifigance values (between 0 and 1) as the first column, and optionally transcript fold change values as the second column. Defaults to NULL
#' @param metaboliteResult Dataframe with metabolite IDs as rownames, metabolite signifigance values (between 0 and 1) as the first column, and optionally metabolite fold change values as the second column. Defaults to NULL
#' @param pathway.id KEGG pathway name to perform bootstrap on.  ("hsa04727")
#' @param seed Random number seed for reproducibility. Defaults to 1990
#' @param numSims Integer indicating the number of bootstrap simulations to perform.  Defaults to 10000
#' @keywords IPA
#' @import manipulate
#' @import ggplot2
#' @export
#' @return Returns an interactive histogram with pathway representation frequencies from bootstrap analysis with manipulate tickers for gene and metabolite significance cutoffs.
#' @examples 
#' data(kData)
#' data(rData)
#' ## Not run: manipulateBootstrapPath("hsa04727", kData, rData)
#' ## End(Not run)
#' @note
#' manipulateBootstrapPath is powered by the following open source databases.  Commercial use and/or redistribution may restricted.  Please see respective terms of use pages and citations for more details.
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



manipulateBootstrapPath=function(pathway.id, geneResult=NULL, metaboliteResult=NULL, seed=1990, numSims=1000){
  data(genes, envir=environment())
  data(metabolites, envir=environment())
  genes<-genes
  metabolites<-metabolites
  
  
  bootstrapPath = function(pathway.id, geneResult=NULL, metaboliteResult = NULL, geneCutoff=NULL, metaboliteCutoff=NULL, seed = 1990, numSims=10000){
   
    if(!is.null(geneResult)&is.null(metaboliteResult)){
      
      #filter testResult by FDR cutoff specified
      sigList <- row.names(subset(geneResult, geneResult[,1] <= geneCutoff))
      #filter pathway data for significant markers
      sigFilter <- row.names(genes)%in%sigList
      #calculate number of significant markers relevant to pathway.id
      sigNum <- sum(genes[sigFilter,pathway.id], na.rm=TRUE)
      #find total number of markers in background relevant to pathway
      totFilter <- row.names(genes)%in%row.names(geneResult)
      totNum <- sum(genes[totFilter,pathway.id], na.rm=TRUE)
      #generate a simulated data frame with dimensions identical to testResult for relevenat pathway. 1's in column two are equivalent to total number of markers relevant and zeros are all other markers.
      simulateFrame <- matrix(nrow = nrow(geneResult), ncol=2, data = c(1:nrow(geneResult), c(rep(1, totNum), rep(0, nrow(geneResult)-totNum))))
      set.seed(seed)
      #scramble simulated data frame
      simulateFrame <- simulateFrame[sample.int(nrow(simulateFrame)),]
      #sample from simulated data frame at a size equivalent to total markers above cutoff
      result<- data.frame()
      set.seed(seed)
      for(i in 1:numSims){
        result[i,1]<-sum(sample(simulateFrame[,2], length(sigList), replace=TRUE))
      }
      #calculate the probability of achieving pathway representation by chance
      probability = length(which(result[,1]>=sigNum))/nrow(result)
      #plot histogram of bootstrap analysis
      local.env<-environment()
      histogram <- ggplot2::ggplot(result, ggplot2::aes(x = result[,1]), environment=local.env) + ggplot2::geom_histogram(fill = "grey", color = "black", binwidth = 1) + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = sigNum, color = "red", lwd = 2) + ggplot2::labs(title = paste("Histogram of Bootstrap Simulation of ", pathway.id,  "\n probability = ", probability, sep=""), x = "Number of Pathway Hits", y = "Frequency") + ggplot2::xlim(0,max(max(result[,1],sigNum)))
      return(histogram)
    }
    
    if(!is.null(metaboliteResult)&is.null(geneResult)){
      
      #filter testResult by FDR cutoff specified
      sigList <- row.names(subset(metaboliteResult, metaboliteResult[,1] <= metaboliteCutoff))
      #filter pathway data for significant markers
      sigFilter <- row.names(metabolites)%in%sigList
      #calculate number of significant markers relevant to pathway.id
      sigNum <- sum(metabolites[sigFilter,pathway.id], na.rm=TRUE)
      #find total number of markers in background relevant to pathway
      totFilter <- row.names(metabolites)%in%row.names(metaboliteResult)
      totNum <- sum(genes[totFilter,pathway.id], na.rm=TRUE)
      #generate a simulated data frame with dimensions identical to testResult for relevenat pathway. 1's in column two are equivalent to total number of markers relevant and zeros are all other markers.
      simulateFrame <- matrix(nrow = nrow(metaboliteResult), ncol=2, data = c(1:nrow(metaboliteResult), c(rep(1, totNum), rep(0, nrow(metaboliteResult)-totNum))))
      set.seed(seed)
      #scramble simulated data frame
      simulateFrame <- simulateFrame[sample.int(nrow(simulateFrame)),]
      #sample from simulated data frame at a size equivalent to total markers above cutoff
      result<- data.frame()
      set.seed(seed)
      for(i in 1:numSims){
        result[i,1]<-sum(sample(simulateFrame[,2], length(sigList), replace=TRUE))
      }
      #calculate the probability of achieving pathway representation by chance
      probability = length(which(result[,1]>=sigNum))/nrow(result)
      
      #plot histogram of bootstrap analysis
      local.env<-environment()
      histogram <- ggplot2::ggplot(result, ggplot2::aes(x = result[,1]), environment=local.env) + ggplot2::geom_histogram(fill = "grey", color = "black", binwidth = 1) + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = sigNum, color = "red", lwd = 2) + ggplot2::labs(title = paste("Histogram of Bootstrap Simulation of ", pathway.id,  "\n probability = ", probability, sep=""), x = "Number of Pathway Hits", y = "Frequency") + ggplot2::xlim(0,max(max(result[,1],sigNum)))
      return(histogram)
    }
    
    if(!is.null(geneResult)&!is.null(metaboliteResult)){
      
      #filter testResult by FDR cutoff specified
      sigList1 <- row.names(subset(geneResult, geneResult[,1] <= geneCutoff))
      #filter pathway data for significant markers
      sigFilter1 <- row.names(genes)%in%sigList1
      #filter pathway data for significant markers
      sigNum1 <- sum(genes[sigFilter1,pathway.id], na.rm=TRUE)
      #find total number of markers in background relevant to pathway
      totFilter1 <- row.names(genes)%in%row.names(geneResult)
      totNum1 <- sum(genes[totFilter1,pathway.id], na.rm=TRUE)
      #generate a simulated data frame with dimensions identical to testResult for relevenat pathway. 1's in column two are equivalent to total number of markers relevant and zeros are all other markers.
      simulateFrame1 <- matrix(nrow = nrow(geneResult), ncol=2, data = c(1:nrow(geneResult), c(rep(1, totNum1), rep(0, nrow(geneResult)-totNum1))))
      set.seed(seed)
      #scramble simulated data frame
      simulateFrame1 <- simulateFrame1[sample.int(nrow(simulateFrame1)),]
      
      #filter testResult by FDR cutoff specified
      sigList2 <- row.names(subset(metaboliteResult, metaboliteResult[,1] <= metaboliteCutoff))
      #filter pathway data for significant markers
      sigFilter2 <- row.names(metabolites)%in%sigList2
      #filter pathway data for significant markers
      sigNum2 <- sum(metabolites[sigFilter2,pathway.id], na.rm=TRUE)
      #find total number of markers in background relevant to pathway
      totFilter2 <- row.names(metabolites)%in%row.names(metaboliteResult)
      totNum2 <- sum(metabolites[totFilter2,pathway.id], na.rm=TRUE)
      #generate a simulated data frame with dimensions identical to testResult for relevenat pathway. 1's in column two are equivalent to total number of markers relevant and zeros are all other markers.
      simulateFrame2 <- matrix(nrow = nrow(metaboliteResult), ncol=2, data = c(1:nrow(metaboliteResult), c(rep(1, totNum2), rep(0, nrow(metaboliteResult)-totNum2))))
      set.seed(seed)
      #scramble simulated data frame
      simulateFrame2 <- simulateFrame2[sample.int(nrow(simulateFrame2)),]
      
      #sample from simulated data frame at a size equivalent to total markers above cutoff
      result<- data.frame()
      set.seed(seed)
      for(i in 1:numSims){
        result[i,1]<-sum(c(sample(simulateFrame1[,2], length(sigList1), replace=TRUE), sample(simulateFrame2[,2], length(sigList2), replace=TRUE)))
      }
      #calculate the probability of achieving pathway representation by chance
      probability = (length(which(result[,1]>=(sigNum1+sigNum2)))/nrow(result))
      
      #plot histogram of bootstrap analysis
      local.env<-environment()
      histogram <- ggplot2::ggplot(result, ggplot2::aes(x = result[,1]), environment=local.env) + ggplot2::geom_histogram(fill = "grey", color = "black", binwidth = 1) + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = sigNum1+sigNum2, color = "red", lwd = 2) + ggplot2::labs(title = paste("Histogram of Bootstrap Simulation of ", pathway.id,  "\n probability = ", probability, sep=""), x = "Number of Pathway Hits", y = "Frequency") + ggplot2::xlim(0,max(max(result[,1],(sigNum1+sigNum2))))
      return(histogram)
      
    }
  }
  myHist = function(geneCutoff, metaboliteCutoff){
    return(bootstrapPath(pathway.id=pathway.id, geneResult=geneResult, metaboliteResult=metaboliteResult, metaboliteCutoff=metaboliteCutoff, geneCutoff=geneCutoff, seed=seed, numSims=numSims))
  }
  return(manipulate::manipulate(myHist(geneCutoff, metaboliteCutoff),geneCutoff=manipulate::slider(0.1,0.5, step=0.05), metaboliteCutoff=manipulate::slider(0.1,0.5, step=0.05)))
}