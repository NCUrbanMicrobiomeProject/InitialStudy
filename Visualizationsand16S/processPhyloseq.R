rm(list=ls())

library("phyloseq")
library("ggplot2")
library("vegan")
library("gridExtra")
library("grid")
library("reshape")
## taxaCols <- c((which(colnames(myT) == "depthAtLevel")+1):(ncol(myT)-2))

#baseDir <- "/Users/mbrown67/Documents/Fodor/Datasets/Wastewater"
if (.Platform$OS.type == "unix") {
  baseDir <- "/Users/mbrown67/Google Drive/Datasets/Wastewater/"  
} else {
  baseDir <- "F:/Google Drive/Datasets/Wastewater"}
#baseDir <- "F:/Google Drive/Datasets/Wastewater"
dataType <- "MetaPhlAn"
metadataDir <- paste(baseDir, "metadata", sep="/")
dataDir <- paste(baseDir, dataType, sep="/")
processedDir <- paste(dataDir, "processed", sep="/")
analysisDir <- paste(dataDir, "analysis", sep="/")
baseFileName <- "MetaPhlAn_LogNormwithMetadata_"
imputeValues <- TRUE

ifelse(dataType != "QIIME", isQIIME <- FALSE, isQIIME <- TRUE)

## Cite this properly
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

grid_arrange_shared_legend_better <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs #Just works for the one legend
  ##Doing this in the general case is hard...
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  disassembledPlots <- lapply(plots, ggplot_build)
  #print(class(disassembledPlots[[1]]))
  #testo$data[[2]]$size <- 1
  #plots <- lapply(disassembledPlots, ggplot_gtable)
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none", line = element_line(size = rel(0.60)),
                                           axis.title.x = element_blank(), axis.title.y = element_blank(),
                                           plot.title = element_blank(),
                                           axis.text.x = element_text(size = rel(1.3))))
  # gl <- lapply(plots, function(x) ifelse(match(x, plots) %% ncol == 1, 
  #                                        x + theme(legend.position="none", 
  #                                           axis.title.x = element_blank(), 
  #                                           plot.title = element_blank(), 
  #                                           axis.text.x = element_text(size = rel(0.6))), 
  #                                        x + theme(legend.position="none", 
  #                                           axis.title.x = element_blank(), axis.title.y = element_blank(), 
  #                                           plot.title = element_blank(), 
  #                                           axis.text.x = element_text(size = rel(0.6)))))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

## Convert the merged data to phyloseq format
## Assumes OTU table for QIIME data
## Assumes per taxonomic level information for MetaPhlAn data
makeNice <- function(myT, taxaCols, thisLevel) {
  otunames <- colnames(myT)[taxaCols]
  otumat <- data.matrix(myT[,taxaCols])
  rownames(otumat) <- myT[,1]
  ## This won't be useful for MetaPhlAn
  if (isQIIME == TRUE) {
    taxamat <- strsplit(colnames(otumat), split="k__|.p__|.c__|.o__|.f__|.g__|.s__")
    ## Need a correction for the length 7 lists
    taxamat <- lapply(taxamat, function(x) if(length(x) == 7){append(x, "")} else {x<-x})
    taxadf <- t(as.data.frame(taxamat))
    taxadf <- taxadf[,-1]
    ## convert it to a dataframe
    colnames(taxadf) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    rownames(taxadf) <- otunames
  } else {
    taxamat <- as.matrix(colnames(otumat), ncol = 1)
    rownames(taxamat) <- taxamat[,1]
    colnames(taxamat) <- thisLevel
    taxadf <- taxamat
    ## taxadf <- as.data.frame(taxamat)
  }
  OTU = otu_table(otumat, taxa_are_rows = FALSE)
  TAX = tax_table(taxadf)
  physeq = phyloseq(OTU, TAX)
  ## plot_bar(physeq, fill = "Phylum")
  
  sampledata <- sample_data(data.frame(
    myT[,c(2:(min(taxaCols) - 1), max(taxaCols) + 1, max(taxaCols) + 2)],
    row.names = sample_names(physeq),
    stringsasFactors = TRUE
  ))
  physeq1 = merge_phyloseq(physeq, sampledata)
  ## This might be a good place for it.
}

eachBetaDiversity <- function(physeqobj, maxAx, taxaLevel) {
  ## Calculates and outputs the plots and values of the ordination
  ## It supports a capscale like approach
  ## CAP works with a formula, but I don't fully understand that yet
  ## Maybe should have some kind of taxonomic level labeling convention
  ## Change from "CAP" to "PCoA"
  returnList <- list()
  physeqobj.bray <- ordinate(physeqobj, "PCoA", "bray", ~1)
  ## May need to also output percentages
  percentVariance <- physeqobj.bray$CA$eig/sum(physeqobj.bray$CA$eig)
  
  write.table(percentVariance, sep = "\t", file=paste0(taxaLevel, "_PCoABray", "_percentVariance.txt"))
  write.table(physeqobj.bray$CA$u, sep="\t", file=paste0(taxaLevel, "_PCoABray", "_pcoa.txt"))
  write.table(physeqobj.bray$CA$eig, file=paste0(taxaLevel, "_PCoABray", "_eigenValues.txt"), sep="\t")
  
  returnList[[1]] <- physeqobj.bray
  pdf(paste0(taxaLevel, "_PCoABrayOrdination.pdf"))
  ## p1 = plot_ordination(physeqobj, physeqobj.bray, type="taxa", color="Phylum", title="taxa")
  ## print(p1)
  for(i in 1:maxAx) {
    for(j in 2:maxAx) {
      if( i == j | i > j) {
        break
      }
      ## Will need some final tweaking here
      ## The color might be botched; Can we easily convert this?
      ## He doesn't use this to plot anyways...
      plot_ordination(physeqobj, physeqobj.bray, type="taxa", color="Phylum", title=taxaLevel, axes = i:j)
    }
  }
  dev.off()
  
  if(isQIIME == TRUE) {
    ## Still need to try to handle the percent variances here
    physeqobj.weun <- ordinate(physeqobj, "PCoA", "unifrac", weighted=TRUE)
    returnList[[2]] <- physeqobj.weun
    ## May need to also output percentages
    write.table(physeqobj.weun$CA$u, sep="\t", file=paste0(taxaLevel, "_capscaleWeightedUnifrac", "_pcoa.txt"))
    write.table(physeqobj.weun$CA$eig, file=paste0(taxaLevel, "_capscaleWeightedUnifrac", "_eigenValues.txt"), sep="\t")
    
    pdf(paste0(taxaLevel, "_capscaleWeightedUnifracOrdination.pdf"))
    ## p1 = plot_ordination(physeqobj, physeqobj.bray, type="taxa", color="Phylum", title="taxa")
    ## print(p1)
    for(i in 1:maxAx) {
      for(j in 2:maxAx) {
        if( i == j | i > j) {
          break
        }
        ## Will need some final tweaking here
        plot_ordination(physeqobj, physeqobj.weun, type="taxa", color="Phylum", title=taxaLevel, axes = i:j)
      }
    }
    dev.off()
    
    physeqobj.unun <- ordinate(physeqobj, "PCoA", "unifrac", weighted=FALSE)
    returnList[[3]] <- physeqobj.unun
    write.table(physeqobj.unun$CA$u, sep="\t", file=paste0(taxaLevel, "_capscaleUnweightedUnifrac", "_pcoa.txt"))
    write.table(physeqobj.unun$CA$eig, file=paste0(taxaLevel, "_capscaleUnweightedUnifrac", "_eigenValues.txt"), sep="\t")
    
    pdf("capscaleUnweightedUnifracOrdination.pdf")
    ## p1 = plot_ordination(physeqobj, physeqobj.bray, type="taxa", color="Phylum", title="taxa")
    ## print(p1)
    for(i in 1:maxAx) {
      for(j in 2:maxAx) {
        if( i == j | i > j) {
          break
        }
        ## Will need some final tweaking here
        plot_ordination(physeqobj, physeqobj.unun, type="taxa", color="Phylum", title="taxa", axes = i:j)
      }
    }
    dev.off()
  }
  ## Moving over to the alpha diversity metrics
  ## Can't use phyloseq for this because they expect integer counts
  return(returnList)
}

taxonomicLevels <- c("phylum", "class", "order", "family", "genus", "species")

## The for loop would start here if we want to process for multiple levels
iterTaxa <- 1
thisLevel <- taxonomicLevels[iterTaxa]

setwd(processedDir)
fileName <- paste0(baseFileName, thisLevel, ".txt")
myT <- read.table(fileName, header=TRUE, comment.char = "@")

if (dataType == "shortBRED") {
  myColClasses <- c(rep("character", 6), rep("numeric", 3), "character", rep("numeric", 3), rep("character", 24), "numeric", rep("character", 10), rep("numeric", ncol(myT) - 48))
} else if (dataType == "qiimeopenmerge") {
  myColClasses <- c(rep("character", 6), rep("numeric", 3), "character", rep("numeric", 3), rep("character", 23), "numeric", rep("character", 9), rep("numeric", ncol(myT) - 46))
} else if (dataType == "MetaPhlAn") {
  myColClasses <- c(rep("character", 6), rep("numeric", 3), "character", rep("numeric", 3), rep("character", 23), "numeric", rep("character", 10), rep("numeric", ncol(myT) - 47))
}

myT <- read.table(fileName, sep = "\t", header=TRUE, colClasses = myColClasses, comment.char="@")
myT[myT == "n/a"] <- NA

## This is true for MetaPhlAn
endMetadataIndex <- which(colnames(myT) == "SampleProject")

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

## We may not always want to allow for imputed values
## Drop Ertapenem and Amoxicillin from consideration because of problems with their detection via calibration and suspected low quantities anyhow.

if (imputeValues == TRUE){
  msCorrect <- read.table(paste(metadataDir,"LODLOQcorrection.txt",sep="/"), sep="\t", header=TRUE)
  ## Will need to run models with and without this effect
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
} else {
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
## The last two columns will become metadata
justFirst <- unlist(lapply(strsplit(as.character(myT$Sample.Date.Time), " "), "[[", 1))
justFirst<-gsub("/16", "/2016", justFirst)
justFirst <- as.numeric(as.POSIXlt(strptime(justFirst, "%m/%d/%Y")))
myT$Timestamp <- justFirst
timeofDay <- unlist(lapply(strsplit(as.character(myT$Sample.Date.Time), " "), "[[", 3))
timeofDay <- format(strptime(timeofDay, format="%I:%M%p"), format = "%H:%M")
myT$timeofDay <- as.numeric(unlist(lapply(strsplit(timeofDay, ":"),"[[",1)))*60 + as.numeric(unlist(lapply(strsplit(timeofDay, ":"),"[[",2)))

myT$Humidity <- as.numeric(gsub("%", "", myT$Humidity))/100

## Properly import these columns
## impute the antibiotic values
## Select the highest depth replicate
## build those other values

setwd(analysisDir)

## Call to makeNice
## Provide it with the taxaCols only
taxaCols <- (endMetadataIndex + 1):(ncol(myT) - 2)
niceMyT <- makeNice(myT, taxaCols, thisLevel)
betaVals <- eachBetaDiversity(niceMyT, 4, "Phylum")
## Chao1 is not applicable to non-count data
alphaShannonDiversity <- apply(myT[, taxaCols], 1, diversity)
alphaT <- cbind(myT, alphaShannonDiversity)
alphaT$Sample <- as.factor(alphaT$Sample)
alphaT$Sample <- factor(alphaT$Sample, levels = c("MNT A", "MNT B", "UW A", "UW B", "UP A", "UP B", "HOSP", "RES", "INF", "PCI", "PCE", "ATE", "FCE", "UV", "DS A", "DS B"))
## Need to remember to report percentages of variance explained
betaT <- cbind(myT, betaVals[[1]]$vectors[,1:3])
betaT$Sample <- as.factor(betaT$Sample)
betaT$Sample <- factor(betaT$Sample, levels = c("MNT A", "MNT B", "UW A", "UW B", "UP A", "UP B", "HOSP", "RES", "INF", "PCI", "PCE", "ATE", "FCE", "UV", "DS A", "DS B"))
alphaT$Type<-factor(ifelse(alphaT$Sample %in% c("MNT A", "MNT B", "UW A", "UW B", "UP A", "UP B"), "Pristine", 
                    ifelse(alphaT$Sample %in% c("HOSP", "RES", "INF", "PCE"), "Pre-\nTreatment", 
                           ifelse(alphaT$Sample %in% c("PCI", "ATE", "UV", "FCE"), "Treatment", "Released"))))
alphaT$Type <- factor(alphaT$Type, levels(alphaT$Type)[c(2, 1, 4, 3)])

betaT$Type<-factor(ifelse(betaT$Sample %in% c("MNT A", "MNT B", "UW A", "UW B", "UP A", "UP B"), "Pristine", 
                           ifelse(betaT$Sample %in% c("HOSP", "RES", "PCE", "INF"), "Pre-\nTreatment", 
                                  ifelse(betaT$Sample %in% c("PCI", "ATE", "UV", "FCE"), "Treatment", "Released"))))
betaT$Type <- factor(betaT$Type, levels(betaT$Type)[c(2, 1, 4, 3)])

alphaT$LocationType <- factor(ifelse(alphaT$Location %in% c("Uwharrie Control", "Mountain Control"), "Rural", 
                                     ifelse(alphaT$Location == "Sugar Creek", "Sugar\nCreek", "Mallard\nCreek")))
alphaT$LocationType <- factor(alphaT$LocationType, levels(alphaT$LocationType)[c(3, 2, 1)])

betaT$LocationType <- factor(ifelse(betaT$Location %in% c("Uwharrie Control", "Mountain Control"), "Rural", 
                                     ifelse(betaT$Location == "Sugar Creek", "Sugar\nCreek", "Mallard\nCreek")))
betaT$LocationType <- factor(betaT$LocationType, levels(betaT$LocationType)[c(3, 2, 1)])
                              

## Do the plotting like Kevin did?

## alphaShannonPlot <- ggplot(alphaT, aes(x=Timepoint, y=alphaShannonDiversity, color=Sample, shape=Location)) +
##   labs(title="Shannon Index Alpha Diversity", x="Timepoint", y="Shannon Index") +
##     geom_jitter(alpha=1, position = position_jitter(width = .2)) +
##     theme(plot.title = element_text(hjust = 0.5))

alphaShannonPlot <- ggplot(alphaT, aes(group=Timepoint, x=Timepoint, y=alphaShannonDiversity, shape=Location, color=Sample))+
  geom_boxplot(outlier.shape = NA, notch = FALSE, show.legend = FALSE)+
  labs(title=paste0("Shannon Diversity Index at Level ", thisLevel), x="Timepoint", y="Shannon Diversity Index")+
  geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("AlphaDiversityShannonbyTimepointatLevel_", thisLevel, ".pdf"), plot = alphaShannonPlot)

alphaShannonPlot2 <- ggplot(alphaT, aes(group=LocationType, x=LocationType, y=alphaShannonDiversity, shape=Location, color=Sample))+
  geom_boxplot(outlier.shape = NA, notch = FALSE, show.legend = FALSE)+
  labs(title=paste0("Shannon Diversity Index at Level ", thisLevel), x="Location", y="Shannon Diversity Index")+
  geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("AlphaDiversityShannonbyLocationatLevel_", thisLevel, ".pdf"), plot = alphaShannonPlot2)

alphaShannonPlot2_SM <- ggplot(alphaT, aes(group=LocationType, x=LocationType, y=alphaShannonDiversity, shape=Location, color=Sample))+
       geom_boxplot(outlier.shape = NA, notch = FALSE, show.legend = FALSE)+xlim("Sugar\nCreek", "Mallard\nCreek")+
       labs(title=paste0("Shannon Diversity Index at Level ", thisLevel), x="Location", y="Shannon Diversity Index")+
       geom_jitter(alpha=1, position = position_jitter(width = .2)) +
       theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("AlphaDiversityShannonbyLocation_SMonly_atLevel_", thisLevel, ".pdf"), plot = alphaShannonPlot2_SM)

alphaShannonPlot3 <- ggplot(alphaT, aes(group=Type, x=Type, y=alphaShannonDiversity, shape=Location, color=Sample))+
  geom_boxplot(outlier.shape = NA, notch = FALSE, show.legend = FALSE)+
  labs(title=paste0("Shannon Diversity Index at Level ", thisLevel), x="Sample Type", y="Shannon Diversity Index")+
  geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("AlphaDiversityShannonbyTypeatLevel_", thisLevel, ".pdf"), plot = alphaShannonPlot3)
## Check with Kevin to see if he likes this or not.

percentVariance <- betaVals[[1]]$values$Eigenvalues/sum(betaVals[[1]]$values$Eigenvalues)

plotPC1byTimepoint <- ggplot(betaT, aes(x=Timepoint, group=Timepoint, y=Axis.1, color=Sample, shape=Location))+
  geom_boxplot(outlier.shape = NA, notch = FALSE, show.legend = FALSE)+
  labs(title=paste0("Bray-Curtis Beta Diversity at Level ", thisLevel), y=paste0("PC1 ", format(percentVariance[1]*100, digits = 3),"%")) +
  geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("BetaDiversityBrayCurtisbyTimepointat_PC1_Level_", thisLevel, ".pdf"), plot  = plotPC1byTimepoint)

plotPC1byLocation <- ggplot(betaT, aes(x=LocationType, group=LocationType, y=Axis.1, color=Sample, shape=Location))+
  geom_boxplot(outlier.shape = NA, notch = FALSE, show.legend = FALSE)+
  labs(title=paste0("Bray-Curtis Beta Diversity at Level ", thisLevel), y=paste0("PC1 ", format(percentVariance[1]*100, digits = 3),"%")) +
  geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("BetaDiversityBrayCurtisbyLocationat_PC1_Level_", thisLevel, ".pdf"), plot  = plotPC1byLocation)

plotPC1byType <- ggplot(betaT, aes(x=Type, group=Type, y=Axis.1, color=Sample, shape=Location))+
  geom_boxplot(outlier.shape = NA, notch = FALSE, show.legend = FALSE)+
  labs(title=paste0("Bray-Curtis Beta Diversity at Level ", thisLevel), y=paste0("PC1 ", format(percentVariance[1]*100, digits = 3),"%")) +
  geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("BetaDiversityBrayCurtisbyTypeat_PC1_Level_", thisLevel, ".pdf"), plot  = plotPC1byType)

comboplot <- grid_arrange_shared_legend_better(alphaShannonPlot, alphaShannonPlot2, alphaShannonPlot3, 
                                        plotPC1byTimepoint, plotPC1byLocation, plotPC1byType,
                                        ncol=3, nrow=2, position="right")

ggsave(paste0("ComboPlot_", thisLevel, ".pdf"), plot = comboplot, scale = 2)


plotPC12 <- ggplot(betaT, aes(x=Axis.1, y=Axis.2, color=Sample, shape=Location)) +
  labs(title=paste0("Bray-Curtis Beta Diversity at Level ", thisLevel), 
       x=paste0("PC1 ", format(percentVariance[1]*100, digits = 3),"%"), 
       y=paste0("PC2 ", format(percentVariance[2]*100, digits = 3),"%")) +
  #geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("BetaDiversityBrayCurtisat_PC1_PC2_Level_", thisLevel, ".pdf"), plot  = plotPC12)

plotPC12byType <- ggplot(betaT, aes(x=Axis.1, y=Axis.2, color=Type, shape=Location)) +
  labs(title=paste0("Bray-Curtis Beta Diversity at Level ", thisLevel), 
       x=paste0("PC1 ", format(percentVariance[1]*100, digits = 3),"%"), 
       y=paste0("PC2 ", format(percentVariance[2]*100, digits = 3),"%")) +
  #geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("BetaDiversityBrayCurtisbyTypeat_PC1_PC2_Level_", thisLevel, ".pdf"), plot  = plotPC12byType)

plotPC12byTimepoint <- ggplot(betaT, aes(x=Axis.1, y=Axis.2, color=Timepoint, shape=Location)) +
  labs(title=paste0("Bray-Curtis Beta Diversity at Level ", thisLevel), 
       x=paste0("PC1 ", format(percentVariance[1]*100, digits = 3),"%"), 
       y=paste0("PC2 ", format(percentVariance[2]*100, digits = 3),"%")) +
  #geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("BetaDiversityBrayCurtisbyTimepointat_PC1_PC2_Level_", thisLevel, ".pdf"), plot  = plotPC12byTimepoint)



plotPC12byLocationType <- ggplot(betaT, aes(x=Axis.1, y=Axis.2, color=LocationType, shape=Location)) +
  labs(title=paste0("Bray-Curtis Beta Diversity at Level ", thisLevel), 
       x=paste0("PC1 ", format(percentVariance[1]*100, digits = 3),"%"), 
       y=paste0("PC2 ", format(percentVariance[2]*100, digits = 3),"%")) +
  #geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("BetaDiversityBrayCurtisbyLocationTypeat_PC1_PC2_Level_", thisLevel, ".pdf"), plot  = plotPC12byLocationType)

## The below doesn't make any sense really
update_geom_defaults("point", list(size = 2))
comboplotPCoA <- grid_arrange_shared_legend_better(plotPC12byTimepoint, plotPC12byLocationType,
                                                   plotPC12byType, plotPC12,
                                                   ncol=2, nrow=2, position="right")

ggsave(paste0("ComboPlotPCoA_", thisLevel, ".pdf"), plot = comboplotPCoA, scale = 2)

plotPC13 <- ggplot(betaT, aes(x=Axis.1, y=Axis.3, color=Sample, shape=Location))+
  labs(title=paste0("Bray-Curtis Beta Diversity at Level ", thisLevel), x=paste0("PC1 ", format(percentVariance[1]*100, digits = 3),"%"), y=paste0("PC3 ", format(percentVariance[3]*100, digits = 3),"%")) +
  #geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("BetaDiversityBrayCurtisat_PC1_PC3_Level_", thisLevel, ".pdf"), plot  = plotPC13)

plotPC23 <- ggplot(betaT, aes(x=Axis.2, y=Axis.3, color=Sample, shape=Location))+
  labs(title=paste0("Bray-Curtis Beta Diversity at Level ", thisLevel),  x=paste0("PC2 ", format(percentVariance[2]*100, digits = 3), "%"), y=paste0("PC3 ", format(percentVariance[3]*100, digits = 3),"%")) +
  #geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("BetaDiversityBrayCurtisat_PC2_PC3_Level_", thisLevel, ".pdf"), plot  = plotPC23)

indexer <- c("Sample", "Location", "Timepoint")
for (indexGuy in indexer) {
  ip <- ggplot(betaT, aes_string("Axis.1", "Axis.2", color=indexGuy)) + geom_point()
  plot(ip)
}

preplot <- ggplot(betaT, aes(Axis.1, Axis.2))
for (indexGuy in indexer) {
  plot(preplot + geom_point(mapping=aes_string(color=indexGuy)))
}

ggplot(betaT, aes(Axis.1, Axis.2, color=Sample, shape=Location, group=interaction(Sample,Location))) + geom_point() + geom_line()
ggplot(betaT, aes(Axis.1, Axis.2, color=Sample, shape=Location, group=interaction(Sample,Location))) + geom_point() + geom_line() + geom_text(aes(label=Timepoint))

ggplot(betaT, aes(Axis.1, Axis.2, color=Sample, shape=Location, group=interaction(Sample,Location))) + geom_point() + geom_path() + geom_text(aes(label=as.numeric(Timepoint)))

simpleAlpha <- alphaT[,c(1, 48:69)]
simpleAlpha <- melt(alphaT[,c(1, 48:69)])
ggplot(simpleAlpha, aes(value, fill=variable)) + geom_histogram()

ggplot(simpleAlpha, aes(Row.names, value)) + geom_point() + facet_wrap(~variable, scales="free_y")

# 15-color palette from http://mkweb.bcgsc.ca/biovis2012/

safeColors <- c("0 0 0", "0 73 73", "0 146 146", "255 109 182", "255 182 219", 
                  "73 0 146", "0 109 219", "182 109 255", "109 182 255", "182 219 255",
                  "146 0 0", "146 73 0", "219 109 0", "36 255 36", "255 255 109",
                "220 220 220")
# Color converter
hexColors <- sapply(strsplit(safeColors, " "), function(x)
  rgb(x[1], x[2], x[3], maxColorValue=255))

plotPC12recolor <- ggplot(betaT, aes(x=Axis.1, y=Axis.2, color=Sample, shape=Location)) +
  labs(title=paste0("Bray-Curtis Beta Diversity at Level ", thisLevel), 
       x=paste0("PC1 ", format(percentVariance[1]*100, digits = 3),"%"), 
       y=paste0("PC2 ", format(percentVariance[2]*100, digits = 3),"%")) +
  #geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  geom_point() +
  scale_color_manual(values = hexColors) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("BetaDiversityBrayCurtisat_PC1_PC2_Level_Recolor_", thisLevel, ".pdf"), plot  = plotPC12recolor, width=7, height = 7)

alphaShannonPlotrecolor <- ggplot(alphaT, aes(group=Timepoint, x=Timepoint, y=alphaShannonDiversity, shape=Location, color=Sample))+
  geom_boxplot(outlier.shape = NA, notch = FALSE, show.legend = FALSE)+
  labs(title=paste0("Shannon Diversity Index at Level ", thisLevel), x="Timepoint", y="Shannon Diversity Index")+
  geom_jitter(alpha=1, position = position_jitter(width = .2)) +
  scale_color_manual(values = hexColors) + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("AlphaDiversityShannonbyTimepointatLevel_Recolor_", thisLevel, ".pdf"), plot = alphaShannonPlotrecolor, width=7, height=7)