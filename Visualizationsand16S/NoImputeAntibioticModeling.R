# No Imputation Modeling
rm(list=ls())
## library("phyloseq")

elsewhere <- c("MNT A", "MNT B", "UW A", "UW B")
## UP versus elsewhere k
## FCE versus DS (DSA and DSB) k (not sure if that works)
## FCE to DS A k
## FCE to DS B k
## DS A to DS B k
## That is where the concern about multiple hypothesis testing comes from
## INF to FCE k
## Cartoon of elimination of taxa
## Those are that present of the INF
## Run some at the level of ShortBRED genes
## Map ShortBRED results to a tree
runCompare <- list(c("UP", "RURAL"),
                   c("HOSP", "RES"), c("INF", "PCE"), c("PCE", "ATE"),
                   c("ATE", "FCE"), c("FCE", "DS A"), c("FCE", "DS B"),
                   c("FCE", "DOWN"), c("DS A", "DS B"), c("INF", "FCE"),
                   c("UP", "DOWN")
                   ##c("UP", "DOWN"), c("UP A", "UP B"), c("FCE", "DOWN"),
                   ##c("HOSP", "RES"), c("HOSP", "INF"), c("RES", "INF"),
                   ##c("FCE", "DS A"), c("FCE", "DS B"), c("DS A", "DS B"),
                   ##c("INF", "PCE"), c("PCE", "ATE"), c("ATE", "FCE"), c("INF", "FCE")##,
                   
                   ##c("FCE", "UV", "Sugar Creek"), c("UV", "DS A", "Sugar Creek"),
                   ##c("INF", "PCI", "Mallard Creek"), c("PCI", "PCE", "Mallard Creek"),
                   ##"INF", "PCE", "ATE", "FCE"
)

## runCompare <- list(c("INF", "FCE"), c("PCI", "PCE", "Mallard Creek"), "FCE")
## This is still kind of redundance right now
## runCompare <- c("ATE", "FCE")
## runCompare <- list(c("UP", "RURAL"))
## The above still doesn't use mtext correctly
## runCompare <- list(c("UP", "DOWN"))
imputeValues <- 0

## This is the model we feed in
## expRel <- c("Location + Sample + Timepoint")
## Probably still helpful for MLM as the MLM term might not change much
## abundanceCols <- c((endMetadataIndex+1):ncol(myT))
## df <- myT
## "(?<=[?.!|])", perl=TRUE
## Requires all of the data to be interpreting by myColClasses
## It won't be foolproof as not all samples exist for every timepoint
## summary versus anova p-values

## Is this generic enough for Shannon diversity or MDS axes?
functionlm <- function(df, responseCols, explainString, rareCut) {
  ## passThreshold <- colSums(df[,responseCols] != 0) >= (nrow(df)/rareCut)
  ## sapply(df[,abundanceCols], function(x) sum(x != 0) >= nrow(df)/rareCut)
  ## The below may not be the right way anymore, but NOTE HERE
  ## keptColNames <- colnames(df)[responseCols[sapply(df[,responseCols], function(x) sum(x != 0) >= nrow(df)/rareCut)]]
  keptColNames <- list(as.data.frame(df[,responseCols]))
  ## print(keptColNames)
  ## Making a sub dataframe might work as well
  ## explainVar <- unlist(strsplit( gsub("([\\+|\\*]+)","~\\1~",gsub(" ", "",explainString)), "~" ))
  ## funlm <- lm(as.formula(paste0("bug ~ ", expRel)), data = df)
  ## funlm <- lm(as.formula(paste0(keptColNames[1]," ~ ", expRel)), data = df)
  
  manylms <- lapply(keptColNames, function(x) lm(as.formula(paste0(x, " ~ ", expRel)), data = df))
  ## testpvals <- lapply(manylms, function (x) coefficients(summary(x))[-1,4])
  ## print("Ran lm successfully")
  ## print(summary(manylms[[1]]))
  ## print(str(anova(manylms[[1]])))
  ## print(length(expRel))
  
  testpvals <- lapply(manylms, function (x) anova(x)$"Pr(>F)")
  ## print(length(testpvals))
  ## print("Ran anova successfully")
  ## perhaps one such solution using anova
  ## na.omit(anova(manylms[[1]])$"Pr(>F)"[-1])
  ## outpVals <- coefficients(summary(funlm))[-1,4]
  pvaldf <- do.call(rbind, testpvals)
  ## print(pvaldf[1,])
  ## print("The do.call worked")
  if (length(testpvals) == 0) {
    return(0)
  }
  ## print(class(pvaldf))
  pvaldf <- as.matrix(pvaldf[,!apply(is.na(pvaldf), 2, any)])
  ## print("The pvaldf was made")
  ## pvaldf <- pvaldf[,-ncol(pvaldf)]
  ## print(class(pvaldf))
  rownames(pvaldf) <- keptColNames
  ## print("rownames worked")
  colnames(pvaldf) <- rownames(anova(manylms[[1]]))[-nrow(anova(manylms[[1]]))]
  ## print("colnames worked")
  ## print(colnames(pvaldf))
  return(pvaldf)
  ## Need a detector for factors
  ## myT[,explainVar[1]]
  ## with may be a solution to contextualize the variable names
}
## pValThreshold <- 0.05
adjustpvals <- function(pvaldf, pValThreshold, helpfulName) {
  ## What does this do by default for NA values
  if(class(pvaldf) == "data.frame") {
    adjdf <- as.data.frame(apply(pvaldf, 2, function(x) p.adjust(x, method="BH")))
  } else {
    return(0)
  }
  ## print(c("adj rownames", rownames(adjdf)))
  ## print(c("original adj says", colnames(adjdf)))
  ## print(class(adjdf))
  ## print(dim(adjdf))
  intName <- sapply(colnames(adjdf), function (x) paste0("adjusted_",x))
  colnames(adjdf) <- intName
  ## print(c( "adj says", colnames(adjdf)))
  mergedf<- cbind(pvaldf, adjdf)
  ## print(dim(mergedf))
  ## print(mergedf[,1])
  ## print(c( "mergedf rownames", rownames(adjdf)))
  ## colnames(mergedf) <- c(colnames(pvaldf), colnames(adjdf))
  ## print(c("mergedf", colnames(mergedf))
  ## Include sigpicker
  ## print("post cbind()")
  sigdf <- as.data.frame(adjdf[apply(adjdf, 1, function(x) any(x < pValThreshold)),])
  ## print("post sigdf creation")
  ## print(c(colnames(adjdf)))
  ## print(colnames(sigdf))
  ## print(dim(sigdf))
  ## print(class(sigdf))
  colnames(sigdf) <- colnames(adjdf)
  ## print(c("sig colnames", colnames(sigdf)))
  mergesigdf <- merge(pvaldf, sigdf, by = 0)
  ## print(mergesigdf[,1])
  rownames(mergesigdf) <- mergesigdf[,1]
  colnames(mergesigdf) <- colnames(mergedf)
  ## I think this is causing the problem
  ## print(c("mergesigdf", colnames(mergesigdf)))
  ## print(pvaldf[1,])
  ## print(mergesigdf[1,])
  ## print(sigdf[1,])
  mergesigdf <- mergesigdf[,-1]
  colnames(mergesigdf) <- c(colnames(pvaldf), colnames(sigdf))
  ## print(c("new mergesigdf", colnames(mergesigdf)))
  ##print(c("mergesigdf", colnames(mergesigdf)))
  write.table(mergedf, file=paste(helpfulName, "_L", taxa, ".txt",sep=""), sep="\t",row.names=TRUE, col.names = NA)
  write.table(mergesigdf, file=paste(helpfulName, "_SigHitsat_L", taxa, ".txt",sep=""), sep="\t", row.names=TRUE, col.names = NA)
  ## print("Finished adjusting")
  dflist <- list(mergedf, mergesigdf)
  ## Also need to write these out.
}
## and now it is a data frame
## pvaldf <- mergesigdf
## dataType <- "qiimeProbably"
## taxa <- 2
## helpfulName <- "whoKnows"
## maybe some flag for sigdf as the histograms are meaningless then
agnoPlotter <- function(pvaldf, df, expRel, helpfulName){
  print(c(locVec,dim(pvaldf)))
  ## print(helpfulName)
  print(expRel)
  print(colnames(pvaldf))
  ## Any trimming of the data.frame must be done beforehand
  ## Same with any processing on the data within the data.frame
  ## df$Sample <- factor(df$Sample, levels = locVec)
  if( dim(pvaldf)[1] > 0) {
    pdf( paste(dataType, "_L_", taxa, "_", expRel, "_", helpfulName, "_boxplots.pdf", sep=""))
    index <- 1
    par(mfrow=c(2,2),
        oma = c(0, 3, 2, 0) + 0.1,
        mar = c(3, 2, 3, 1) + 0.1,
        mgp = c(3, 1.25, 0))
    
    ## print(colnames(pvaldf))
    maxUn <- min(grep("adjusted", colnames(pvaldf))) - 1
    ## minUn <- max(grep("^wilcoxon", colnames(pvaldf))) + 1
    minUn <- 1
    ## print(c(minUn, maxUn))
    ##print("NOT inside the nameIter loop")
    
    for (nameIter in rownames(pvaldf)) {
      ## here are the column names
      expVec <- gsub(" ", "", unlist(strsplit(expRel, split="[\\+|\\*]")))
      thisbug <- df[, nameIter]
      myFrame <- cbind(thisbug, df[expVec])
      plotIndex <- 0
      marDone <- FALSE
      ## print("Inside the nameIter loop")
      for (pIter in colnames(pvaldf)[minUn:maxUn]){
        pIndex <- which(colnames(pvaldf) == pIter) + maxUn - 1
        ## This won't work for the nonparametric tests
        ## Those can be added in manually afterwards
        ## print("Inside the pIter loop")
        if (length(grep(":", pIter)) > 0){
          ## Need to reprocess to drop the subterms
          compos <- expVec[which(sapply(expVec, function(x) grep(x, pIter)) > 0, arr.ind = TRUE)]
          ## compos <- unlist(strsplit(pIter, split=":"))
          ## varCompos <- lapply(compos, function(x) eval(parse(text=x)))
          varCompos <- myFrame[compos]
          ## Right now just solve it for the two-way interaction case
          if(class(varCompos[[1]]) == "character" & class(varCompos[[2]]) == "character"){
            boxplot(thisbug ~ as.factor(varCompos[[1]])*as.factor(varCompos[[2]]), main = paste(pIter, sprintf("%.3e", pvaldf[nameIter, pIndex])), las=1, outline=FALSE)
            stripchart(thisbug ~ as.factor(varCompos[[1]])*as.factor(varCompos[[2]]),
                       vertical = TRUE, method="jitter", pch = 16, add=TRUE)
          } else if(class(varCompos[[1]]) == "character" & class(varCompos[[2]]) == "numeric"){
            plot(varCompos[[2]], thisbug, col=as.factor(varCompos[[1]]), pch=16, main = paste(pIter, sprintf("%.3e", pvaldf[nameIter, pIndex])))
            ## Mixed case
            ## Plotting the lines might be nice here.
          } else if(class(varCompos[[2]]) == "character" & class(varCompos[[1]]) == "numeric"){
            plot(varCompos[[1]], thisbug, col=as.factor(varCompos[[2]]), pch=16, main = paste(pIter, sprintf("%.3e", pvaldf[nameIter, pIndex])))
          } else {
            plot(varCompos[[1]], thisbug, cex=abs(scale(varCompos[[2]])), pch=16, main = paste(pIter, sprintf("%.3e", pvaldf[nameIter, pIndex])))
            ## Both numeric; This is the difficult case
            ## Size by the other continuous variable
            ## Lines of a lot of them could be plotted as well.
          }
        } else {
          ## Will need to strip the non
          ## You will have to be careful with redundant naming conventions
          grabName <- expVec[which(sapply(expVec, function(x) grep(x, pIter)) > 0, arr.ind = TRUE)]
          ## pVar <- eval(parse(text=grabName))
          pVar <- myFrame[,grabName]
          
          if (class(pVar) == "numeric"){
            plot(pVar, thisbug, main = paste(pIter, format(pvaldf[nameIter, pIndex], digits=3)), pch=16)
          } else {
            boxplot(thisbug ~ as.factor(pVar), main = paste(pIter, format(pvaldf[nameIter, pIndex], digits=3)), data=myFrame, las=1, outline=FALSE)
            stripchart(thisbug ~ as.factor(pVar), data=myFrame,
                       vertical = TRUE, pch = 16, add=TRUE, method="jitter")
          }
          ## print("Plotting successful")
        }
        ## This needs to be modulo number of terms in model?
        ## mtext isn't printing for sig case
        
        if(pIndex %% 4 == 0){
          mtext(nameIter, outer=TRUE, cex = 0.7)
          marDone <- TRUE
          if(length(grep("ng.L.", nameIter)) > 0) {
            mtext("Concentration (ng/L)", outer=TRUE, side=2, line=1)
          } else if (dataType == "qiimeclosedmerge"){
            mtext("Log normalized abundance", outer=TRUE, side=2, line=1)
          } else if (dataType == "MetaPhlAn") {
            mtext("MetaPhlAn Relative Abundance", outer=TRUE, side=2, line=1)
          } else if (dataType == "HUMAnN") {
            if (dataType == "pathcoverage") {
              mtext("Coverage Percentage", outer=TRUE, side=2, line=1)
            } else {
              mtext("RPKs", outer=TRUE, side=2, line=1)
            }
          }
        }
        if (marDone == FALSE) {
          mtext(nameIter, outer=TRUE, cex = 0.7)
          if(length(grep("ng.L.", nameIter)) > 0) {
            mtext("Concentration (ng/L)", outer=TRUE, side=2, line=1)
          } else if (dataType == "qiimeclosedmerge"){
            mtext("Log normalized abundance", outer=TRUE, side=2, line=1)
          } else if (dataType == "MetaPhlAn") {
            mtext("MetaPhlAn Relative Abundance", outer=TRUE, side=2, line=1)
          } else if (dataType == "HUMAnN") {
            if (dataType == "pathcoverage") {
              mtext("Coverage Percentage", outer=TRUE, side=2, line=1)
            } else {
              mtext("RPKs", outer=TRUE, side=2, line=1)
            }
          }
          
        }
        plotIndex <- plotIndex + 1
      }
      ## Needs a fix here
      ## while(pIndex %% 4 != 0){
      ##     pIndex <- pIndex + 1
      ##     ## Blank plot from https://stat.ethz.ch/pipermail/r-help/2004-February/045948.html
      ##     plot(1, type="n", axes=FALSE, xlab="", ylab="")
      ## }
      ## Will need to plot blanks to get it to 4 plots
      ## So this messes up non-4
      while(plotIndex %% 4 != 0) {
        plot.new()
        plotIndex <- plotIndex + 1
      }
    }
    ## print(c(maxUn, nameIter))
    for(histIter in 1:maxUn) {
      hist(pvaldf[,histIter], breaks=20, main=colnames(pvaldf)[histIter])
    }
    dev.off()
  }
}

functionCompare <- function(df, locVec, responseCols, explainString, rareCut, pValThreshold, moreHelp=""){
  ## print("Made it into functionCompare successfully")
  subT.this <- df[df$Sample %in% locVec,]
  for (locItem in locVec) {
    if(length(subT.this[subT.this$Sample == locItem,]$Sample) == 0) {
      return(0)
    }
  }
  ## Attempt force ordering to be kept
  ## print(unique(factor(subT.this$Sample)))
  subT.this$Sample <- factor(subT.this$Sample, levels = locVec)
  ## print(unique(factor(subT.this$Sample)))
  ## print("Releveled factor successfully")
  pValz <- functionlm(subT.this, responseCols, explainString, rareCut)
  ## print("Ran functionlm successfully")
  if (pValz == 0) {
    print(paste("Does not work for", locVec))
    return(0)
  }
  adjList <- adjustpvals(pValz, pValThreshold, paste0(moreHelp, "_", dataType, "_", explainString, paste0(locVec, collapse="_"), collapse="_"))
  ## print("Adjusted p-values successfully")
  agnoPlotter(adjList[[1]], subT.this, explainString, paste0(moreHelp, "_", paste0(locVec, collapse="_")))
  ## print("Plotted all graphs successfully")
  ## The sig histograms will not be helpful.
  agnoPlotter(adjList[[2]], subT.this, explainString, paste0(moreHelp, "_", paste0("SIG_", paste0(locVec,  collapse="_"), collapse="_")))
  ## print("Plotted sig graphs successfully")
}

## dataType <- "qiimeclosedmerge"
## dataType <- "MetaPhlAn"
## dataType <- "HUMANn"
dataType <- "shortBRED"
baseDir <- "/Users/mbrown67/Documents/Fodor/Datasets/Wastewater"
## I think the modelDescription is largely meaningless now
modelDescription <- "StreamUpVDownTimepointDepth"
if(dataType == "qiimeclosedmerge") {
  dataDir <- paste(baseDir, "qiime/closed_reference_merged/countByTaxonomy", sep="/")
  isQIIME <- TRUE
} else if (dataType == "shortBRED") {
  dataDir <- paste(baseDir, "shortBRED/quantifyFiles", sep="/")
  isQIIME <- FALSE
} else if (dataType == "MetaPhlAn") {
  dataDir <- paste(baseDir, "MetaPhlAn/abundancetables/tables_allfourtps", sep="/")
  isQIIME <- FALSE
} else if (dataType == "HUMANn") {
  dataDir <- paste(baseDir, "HUMANn/UpdatedData", sep="/")
  isQIIME <- FALSE
}

## helpfulName <- paste0(helpfulName, modelDescription)
metadataDir <- paste(baseDir, "metadata", sep="/")
processedDir <- paste(dataDir, "processed", sep="/")
analysisDir <- paste(dataDir, "analysis", sep="/")

if (dataType == "shortBRED" | dataType == "HUMANn") {
  taxonomicLevels <- c(2)
} else {
  taxonomicLevels <- c(2) # c(2:7)
}

for (taxa in taxonomicLevels) {
  ## Load data
  
  setwd(processedDir)
  if (dataType == "shortBRED") {
    fileName <- paste("shortBREDwMetadata.txt")
  } else if (dataType == "qiimeclosedmerge") {
    fileName <- paste(dataType, "_LogNormwithMetadata_L_",taxa,".txt", sep="")
  } else if (dataType == "MetaPhlAn") {
    if( taxa == 2) {
      fileName <- paste(dataType, "_LogNormwithMetadata_","phylum",".txt", sep="")
    } else if (taxa == 3) {
      fileName <- paste(dataType, "_LogNormwithMetadata_","class",".txt", sep="")
    } else if( taxa == 4) {
      fileName <- paste(dataType, "_LogNormwithMetadata_","order",".txt", sep="")
    } else if (taxa == 5) {
      fileName <- paste(dataType, "_LogNormwithMetadata_","family",".txt", sep="")
    } else if( taxa == 6) {
      fileName <- paste(dataType, "_LogNormwithMetadata_","genus",".txt", sep="")
    } else if (taxa == 7) {
      fileName <- paste(dataType, "_LogNormwithMetadata_","species",".txt", sep="")
    }
  } else {
    ## HUMAnN processing
  }
  
  myT <- read.table(fileName, sep="\t", header=TRUE, comment.char="@")
  ## This now depends on the datatype
  ## This will change when the metadata gets reprocessed
  if (dataType == "shortBRED") {
    myColClasses <- c(rep("character", 6), rep("numeric", 3), "character", rep("numeric", 3), rep("character", 24), "numeric", rep("character", 10), rep("numeric", ncol(myT) - 48))
  } else if (dataType == "qiimeclosedmerge") {
    ## Make sure that this is correct
    myColClasses <- c(rep("character", 6), rep("numeric", 3), "character", rep("numeric", 3), rep("character", 23), "numeric", rep("character", 9), rep("numeric", ncol(myT) - 46))
  } else if (dataType == "MetaPhlAn") {
    myColClasses <- c(rep("character", 6), rep("numeric", 3), "character", rep("numeric", 3), rep("character", 23), "numeric", rep("character", 10), rep("numeric", ncol(myT) - 47))
  }
  
  ## ShortBRED there are 46
  ## MetaPhlAn there are 47
  myT <- read.table(fileName, sep = "\t", header=TRUE, colClasses = myColClasses, comment.char="@")
  myT[myT == "n/a"] <- NA
  if (dataType != "qiimeclosedmerge") {
    endMetadataIndex <- which(colnames(myT) == "SampleProject")
  } else {
    endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
  }
  
  ## myT$Location <- as.factor(myT$Location)
  myT$Timepoint <- as.factor(myT$Timepoint)
  ## myT$Sample <- as.factor(myT$Sample)
  
  ## Relevel the Timepoint variable
  levels(myT$Timepoint)[levels(myT$Timepoint) == "1"] <- "Late\nWinter"
  levels(myT$Timepoint)[levels(myT$Timepoint) == "2"] <- "Early\nSpring"
  levels(myT$Timepoint)[levels(myT$Timepoint) == "3"] <- "Mid\nSpring"
  levels(myT$Timepoint)[levels(myT$Timepoint) == "4"] <- "Mid\nSummer"
  ## myT$Timepoint <- gsub("1", "Late Winter", myT$Timepoint)
  ## myT$Timepoint <- gsub("2", "Early Spring", myT$Timepoint)
  ## myT$Timepoint <- gsub("3", "Mid Spring", myT$Timepoint)
  ## myT$Timepoint <- gsub("4", "Mid Summer", myT$Timepoint)
  
  ## Will later drop Ertapenem and Amoxicillin from consideration because of problems with their detection via calibration and suspected low quantities anyhow.
  if (imputeValues == 1){
    myT$Ertapenem..ng.L. <- gsub("ND", 0, myT$Ertapenem..ng.L.)
    myT$Ertapenem..ng.L. <- gsub("<LOQ", 0, myT$Ertapenem..ng.L.)
    myT$Ertapenem..ng.L. <- as.numeric(myT$Ertapenem..ng.L.)
    myT$Amoxicillin..ng.L. <- gsub("ND", 0, myT$Amoxicillin..ng.L.)
    myT$Amoxicillin..ng.L. <- gsub("<LOQ", 0, myT$Amoxicillin..ng.L.)
    myT$Amoxicillin..ng.L. <- as.numeric(myT$Amoxicillin..ng.L.)
    myT$Ciprofloxacin..ng.L. <- gsub("ND", 0, myT$Ciprofloxacin..ng.L.)
    myT$Ciprofloxacin..ng.L. <- gsub("<LOQ", 0, myT$Ciprofloxacin..ng.L.)
    myT$Ciprofloxacin..ng.L. <- as.numeric(myT$Ciprofloxacin..ng.L.)
    myT$Doxycycline..ng.L. <- gsub("ND", 0, myT$Doxycycline..ng.L.)
    myT$Doxycycline..ng.L. <- gsub("<LOQ", 0, myT$Doxycycline..ng.L.)
    myT$Doxycycline..ng.L. <- as.numeric(myT$Doxycycline..ng.L.)
    myT$Azithromycin..ng.L. <- gsub("ND", 0, myT$Azithromycin..ng.L.)
    myT$Azithromycin..ng.L. <- gsub("<LOQ", 0, myT$Azithromycin..ng.L.)
    myT$Azithromycin..ng.L. <- as.numeric(myT$Azithromycin..ng.L.)
    myT$Clindamycin..ng.L. <- gsub("ND", 0, myT$Clindamycin..ng.L.)
    myT$Clindamycin..ng.L. <- gsub("<LOQ", 0, myT$Clindamycin..ng.L.)
    myT$Clindamycin..ng.L. <- as.numeric(myT$Clindamycin..ng.L.)
    myT$Sulfamethoxazole..ng.L. <- gsub("ND", 0, myT$Sulfamethoxazole..ng.L.)
    myT$Sulfamethoxazole..ng.L. <- gsub("<LOQ", 0, myT$Sulfamethoxazole..ng.L.)
    myT$Sulfamethoxazole..ng.L. <- as.numeric(myT$Sulfamethoxazole..ng.L.)
    myT$Cephalexin..ng.L. <- gsub("ND", 0, myT$Cephalexin..ng.L.)
    myT$Cephalexin..ng.L. <- gsub("<LOQ", 0, myT$Cephalexin..ng.L.)
    myT$Cephalexin..ng.L. <- as.numeric(myT$Cephalexin..ng.L.)
    myT$Trimethoprim..ng.L. <- gsub("ND", 0, myT$Trimethoprim..ng.L.)
    myT$Trimethoprim..ng.L. <- gsub("<LOQ", 0, myT$Trimethoprim..ng.L.)
    myT$Trimethoprim..ng.L. <- as.numeric(myT$Trimethoprim..ng.L.)
    myT$Levofloxacin..ng.L. <- gsub("ND", 0, myT$Levofloxacin..ng.L.)
    myT$Levofloxacin..ng.L. <- gsub("<LOQ", 0, myT$Levofloxacin..ng.L.)
    myT$Levofloxacin..ng.L. <- as.numeric(myT$Levofloxacin..ng.L.)
  } 
  else if (imputeValues == 2){
    msCorrect <- read.table(paste(metadataDir,"LODLOQcorrection.txt",sep="/"), sep="\t", header=TRUE)
    myT$Ertapenem..ng.L. <- gsub("ND", msCorrect$Ertapenem[1], myT$Ertapenem..ng.L.)
    myT$Ertapenem..ng.L. <- gsub("<LOQ", msCorrect$Ertapenem[2], myT$Ertapenem..ng.L.)
    myT$Ertapenem..ng.L. <- as.numeric(myT$Ertapenem..ng.L.)
    myT$Amoxicillin..ng.L. <- gsub("ND", msCorrect$Amoxicillin[1], myT$Amoxicillin..ng.L.)
    myT$Amoxicillin..ng.L. <- gsub("<LOQ", msCorrect$Amoxicillin[2], myT$Amoxicillin..ng.L.)
    myT$Amoxicillin..ng.L. <- as.numeric(myT$Amoxicillin..ng.L.)
    myT$Ciprofloxacin..ng.L. <- gsub("ND", msCorrect$Ciprofloxacin[1], myT$Ciprofloxacin..ng.L.)
    myT$Ciprofloxacin..ng.L. <- gsub("<LOQ", msCorrect$Ciprofloxacin[2], myT$Ciprofloxacin..ng.L.)
    myT$Ciprofloxacin..ng.L. <- as.numeric(myT$Ciprofloxacin..ng.L.)
    myT$Doxycycline..ng.L. <- gsub("ND", msCorrect$Doxycycline[1], myT$Doxycycline..ng.L.)
    myT$Doxycycline..ng.L. <- gsub("<LOQ", msCorrect$Doxycycline[2], myT$Doxycycline..ng.L.)
    myT$Doxycycline..ng.L. <- as.numeric(myT$Doxycycline..ng.L.)
    myT$Azithromycin..ng.L. <- gsub("ND", msCorrect$Azithromycin[1], myT$Azithromycin..ng.L.)
    myT$Azithromycin..ng.L. <- gsub("<LOQ", msCorrect$Azithromycin[2], myT$Azithromycin..ng.L.)
    myT$Azithromycin..ng.L. <- as.numeric(myT$Azithromycin..ng.L.)
    myT$Clindamycin..ng.L. <- gsub("ND", msCorrect$Clindamycin[1], myT$Clindamycin..ng.L.)
    myT$Clindamycin..ng.L. <- gsub("<LOQ", msCorrect$Clindamycin[2], myT$Clindamycin..ng.L.)
    myT$Clindamycin..ng.L. <- as.numeric(myT$Clindamycin..ng.L.)
    myT$Sulfamethoxazole..ng.L. <- gsub("ND", msCorrect$Sulfamethoxazole[1], myT$Sulfamethoxazole..ng.L.)
    myT$Sulfamethoxazole..ng.L. <- gsub("<LOQ", msCorrect$Sulfamethoxazole[2], myT$Sulfamethoxazole..ng.L.)
    myT$Sulfamethoxazole..ng.L. <- as.numeric(myT$Sulfamethoxazole..ng.L.)
    myT$Cephalexin..ng.L. <- gsub("ND", msCorrect$Cephalexin[1], myT$Cephalexin..ng.L.)
    myT$Cephalexin..ng.L. <- gsub("<LOQ", msCorrect$Cephalexin[2], myT$Cephalexin..ng.L.)
    myT$Cephalexin..ng.L. <- as.numeric(myT$Cephalexin..ng.L.)
    myT$Trimethoprim..ng.L. <- gsub("ND", msCorrect$Trimethoprim[1], myT$Trimethoprim..ng.L.)
    myT$Trimethoprim..ng.L. <- gsub("<LOQ", msCorrect$Trimethoprim[2], myT$Trimethoprim..ng.L.)
    myT$Trimethoprim..ng.L. <- as.numeric(myT$Trimethoprim..ng.L.)
    myT$Levofloxacin..ng.L. <- gsub("ND", msCorrect$Levofloxacin[1], myT$Levofloxacin..ng.L.)
    myT$Levofloxacin..ng.L. <- gsub("<LOQ", msCorrect$Levofloxacin[2], myT$Levofloxacin..ng.L.)
    myT$Levofloxacin..ng.L. <- as.numeric(myT$Levofloxacin..ng.L.)
  } else if (imputeValues == 0){
    ## Replace with NAs. Have the modeling skip those
    myT[myT == "ND"] <- NA
    myT[myT == "<LOQ"] <- NA
  }
  
  ## Select highest depth replicate
  savemyT <- myT[FALSE,]
  for (eachVal in unique(myT$Sample.ID)) {
    myTperSample <- myT[myT$Sample.ID == eachVal,]
    if (isQIIME == TRUE){
      rep<-myTperSample[which(as.numeric(myTperSample$sampleDepth) == max(as.numeric(myTperSample$sampleDepth))),]
    }
    else {
      rep<-myTperSample[which(as.numeric(myTperSample$InputReads) == max(as.numeric(myTperSample$InputReads))),]
    }
    savemyT <- rbind(savemyT, rep)
  }
  myT <- savemyT
  
  ## Fixing some oddly formatted values
  ## This will change the metadataend Information
  justFirst <- unlist(lapply(strsplit(as.character(myT$Sample.Date.Time), " "), "[[", 1))
  justFirst<-gsub("/16", "/2016", justFirst)
  justFirst <- as.numeric(as.POSIXlt(strptime(justFirst, "%m/%d/%Y")))
  myT$Timestamp <- justFirst
  timeofDay <- unlist(lapply(strsplit(as.character(myT$Sample.Date.Time), " "), "[[", 3))
  timeofDay <- format(strptime(timeofDay, format="%I:%M%p"), format = "%H:%M")
  myT$timeofDay <- as.numeric(unlist(lapply(strsplit(timeofDay, ":"),"[[",1)))*60 + as.numeric(unlist(lapply(strsplit(timeofDay, ":"),"[[",2)))
  
  myT$Humidity <- as.numeric(gsub("%", "", myT$Humidity))/100
  
  ## The Depth part of this will have to be fixed
  ## This should be repositioned
  rareCut <- 4
  pValThreshold <- 0.05
  ## responseCols <- c((endMetadataIndex+1):(ncol(myT)-2))
  ## responseCols <- c(19:26)
  ## imputeValues <- FALSE
  ## responseCols <- c(19:26)
  
  setwd(analysisDir)
  for(locVec in runCompare) {
    for(antibioticCol in c(19:26)){
    ## responseCols <- c(19:26)
    subT <- myT
    ## subT$Antibiotic <- subT[,antibioticCol]
    ## subT <- subT[complete.cases(subT$Antibiotic),]
    subT <- subT[complete.cases(subT[,antibioticCol]),]
    ## responseCols <- c((endMetadataIndex+1):(ncol(myT)-3))
    ## use to be any
    if (all(locVec %in% elsewhere) == FALSE) {
      subT <- subT[subT$Location == "Mallard Creek" | subT$Location == "Sugar Creek",]
    }
    
    expRel <- c("Location + Sample + Timepoint")
    if (all.equal(c("UP", "DOWN"), locVec) == TRUE) {
      if (length(subT[subT$Sample == "UP A" | subT$Sample == "UP B",]$Sample) > 0) {
        subT[subT$Sample == "UP A" | subT$Sample == "UP B",]$Sample <- "UP"
    } else {
        next
    }
    if (length(subT[subT$Sample == "DS A" | subT$Sample == "DS B",]$Sample) > 0){
      subT[subT$Sample == "DS A" | subT$Sample == "DS B",]$Sample <- "DOWN"
    } else {
      next
    }
    }
    if (all.equal(c("FCE", "DOWN"), locVec) == TRUE) {
      if (length(subT[subT$Sample == "DS A" | subT$Sample == "DS B",]$Sample) > 0){
        subT[subT$Sample == "DS A" | subT$Sample == "DS B",]$Sample <- "DOWN"
      }
      else {
        next
      }
      if (length(subT[subT$Sample == "FCE",]$Sample == 0)) {
        next
      }
    }
    if (all.equal(c("UP", "RURAL"), locVec) == TRUE) {
      ## print("It gets here")
      ## subT <- myT
      if (length(subT[subT$Sample == "UP A" | subT$Sample == "UP B",]$Sample) > 0){
        subT[subT$Sample == "UP A" | subT$Sample == "UP B",]$Sample <- "UP"
      } else {
        next
      }
      if (length(subT[subT$Sample == "MNT A" | subT$Sample == "MNT B" | subT$Sample == "UW A" | subT$Sample == "UW B",]$Sample) > 0) {
      subT[subT$Sample == "MNT A" | subT$Sample == "MNT B" | subT$Sample == "UW A" | subT$Sample == "UW B",]$Sample <- "RURAL"
    } else {
      next
    }
    }
    subT$Location <- as.factor(subT$Location)
    subT$Sample <- as.factor(subT$Sample)
    
    functionCompare(subT, locVec, antibioticCol, expRel, rareCut, pValThreshold, paste0(colnames(myT)[antibioticCol],"Non-ImputedAntibiotics"))
    }
  }
}