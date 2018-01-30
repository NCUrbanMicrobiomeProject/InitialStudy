# Antibiotics like ShortBRED
#
# 16S Visualization like the ShortBRED one
## Sig Pattern visualizer
## ShortBRED, WGS, 16S, Antibiotics
rm(list=ls())

library(tidyverse)
library(ggplot2)

sigThreshold <- 0.05

## Get the required file

if (.Platform$OS.type == "unix") {
  baseDir <- "/Users/mbrown67/Google Drive/Datasets/Wastewater/"  
} else {
  baseDir <- "F:/Google Drive/Datasets/Wastewater"}


#fileName <- "/Users/mbrown67/Google Drive/Urban Environmental Genomics Project/Manuscript 1 Final Comparison Reports/shortbread_models/Consolidated_Markers_Format/pooledModelResults.txt"

# Just worrying about the ShortBRED aspect

setwd(baseDir)
#fileName <- "./shortBRED/quantifyFiles/analysis/shortBREDResults.txt"
fileName <- "./qiime/open_reference_CSSnormalized_merged/analysis/BestHits_L7.txt"
fileName <- "./shortBRED/quantifyFiles/analysis/ImputedNDLOQ_massSpec.txt"
myT <- read.table(fileName, header=TRUE, sep = '\t')

# myTSwap <- 
#print(dim(myT))

#Dropping INF vs FCE
#print(myT$Sites.Compared)
myT <- myT[myT$Sites.Compared != "INF_FCE",]
myT <- myT[myT$Sites.Compared != "UP_RURAL",]
myT <- myT[myT$Sites.Compared != "DS A_DS B",]
myT <- droplevels.data.frame(myT)
#print(dim(myT))
## Reread in metadata to fix the mangled gene names
#fileName <- file("./shortBRED/quantifyFiles/processed/shortBREDwMetadata.txt", "r")
#realColNames <- unlist(strsplit(readLines(fileName, n = 1), '\\"\t\\"'))
#myData <- read.table("./shortBRED/quantifyFiles/processed/shortBREDwMetadata.txt", header=TRUE, sep='\t')

#savemyData <- myData[FALSE,]
#for (eachVal in unique(myData$Sample.ID)) {
#  myDataperSample <- myData[myData$Sample.ID == eachVal,]
#  rep<-myDataperSample[which(as.numeric(myDataperSample$InputReads) == max(as.numeric(myDataperSample$InputReads))),]
#  savemyData <- rbind(savemyData, rep)
#}
#myData <- savemyData

# Drop the "other" locations from consideration
# avgFecalStress <- vector()
# avgFecalControl <- vector()
# avgFecalMale <- vector()
# avgFecalFemale <- vector()
# fecalSexDiff <- vector()
# fecalGroupDiff <- vector()
# 
# subT <- myT[myT$Source == "feces",]
# index <- 1
# for(nameIter in rownames(fecalData)){
#   avgFecalControl[index] <- mean(subT[subT$Group == "Control",][nameIter][,1])
#   avgFecalStress[index] <- mean(subT[subT$Group == "Experimental",][nameIter][,1])
#   fecalGroupDiff[index] <- sign(avgFecalControl[index] - avgFecalStress[index])
#   avgFecalMale[index] <- mean(subT[subT$Sex == "male",][nameIter][,1])
#   avgFecalFemale[index] <- mean(subT[subT$Sex == "female",][nameIter][,1])
#   fecalSexDiff[index] <- sign(avgFecalMale[index] - avgFecalFemale[index])
#   index = index + 1
# }
# fecalData <- cbind(fecalData, avgFecalControl, avgFecalStress, fecalGroupDiff, avgFecalMale, avgFecalFemale, fecalSexDiff)


#close(fileName)
## The next step loses the organism name
## myTARO <- unlist(lapply(strsplit(as.character(myT$names), split = "\\."), function(x) ifelse(length(x) > 3, x[grep("ARO_", x)], x)))
#myTgi <- unlist(lapply(strsplit(as.character(myT$Antibiotic.Resistant.Gene), split = "\\."), function(x) ifelse(length(x) > 3, x[2], x)))
## testo <- unlist(lapply(strsplit(as.character(myT$names), split = "\\."), function(x) ifelse(length(x) > 1, paste(x[grep("ARO_", x)], strsplit(x[length(x)], "__"))[[1]], x)))
## testomatic <- realColNames[grep(paste(myTARO, collapse="|"), realColNames)]
#betterNameList <- list()
#for(key in myTgi) {
#  betterNameList[[key]] <- realColNames[grep(key, realColNames)]
## betterNameList[[key]] <- betterNameList[[key]][1]
## But some of these have multiple organisms associated with the same name...
#}
#myT$Antibiotic.Resistant.Gene <- unlist(betterNameList[myTgi])
## map the names over

## Replace NAs with a value that won't be recognized as a p-value possibility
#print(levels(myT$Sites.Compared))
#myT[is.na(myT)] <- 2
## Reorder "Comparison" factor levels
## myT$Comparison <- gsub("_", "\n", myT$Comparison)
tweakSiteNames <- gsub("_", "\n", myT$Sites.Compared)
tweakSiteNames <- gsub("DS\nA", "DS_A", tweakSiteNames)
tweakSiteNames <- gsub("DS\nB", "DS_B", tweakSiteNames)
tweakSiteNames <- gsub("UP\nA", "UP_A", tweakSiteNames)
tweakSiteNames <- gsub("UP\nB", "UP_B", tweakSiteNames)
myT$Sites.Compared <- as.factor(tweakSiteNames)
myT$Sites.Compared = factor(myT$Sites.Compared, levels(myT$Sites.Compared)[c(5, 8, 6, 7, 1, 3, 4, 2)])
## Add on the color vector
#adjCols <- grep("adjusted", colnames(myT))
adjCols <- c(3,4,5)
## This doesn't grab rows with NA columns correctly

## Clean up the text display
## The ShortBRED paper just has the GeneName
## We don't need to have this match the outputs from the other supplemental figures
## SHV or TEM

sigdf <- as.data.frame(myT[apply(myT[,adjCols], 1, function(x) any(x < sigThreshold)),])
## We lose two of the comparisons here, so refer to that later in some caption creation for this figure

isSig <- sigdf[adjCols] < sigThreshold
colorVector <- vector()
## This can definitely be streamlined
for(iterz in 1:nrow(sigdf)){
  ## Changes in waterway alone are significant
  if(all(isSig[iterz,] == c(1, 0, 0)))
    colorVector[iterz] <- "red"
  ## Sample site Alone
  if(all(isSig[iterz,] == c(0, 1, 0)))
    colorVector[iterz] <- "green"
  ## Time alone
  if(all(isSig[iterz,] == c(0, 0, 1)))
    colorVector[iterz] <- "blue"
  ## Waterway and Sample site
  if(all(isSig[iterz,] == c(1, 1, 0)))
    colorVector[iterz] <- "yellow"
  ## Waterway and timepoint
  if(all(isSig[iterz,] == c(1, 0, 1)))
    colorVector[iterz] <- "magenta"
  ## Changes in Sample site and time significantly
  if(all(isSig[iterz,] == c(0, 1, 1)))
    colorVector[iterz] <- "cyan"
  if(all(isSig[iterz,] == c(1, 1, 1)))
    colorVector[iterz] <- "black"
}
colorVector <- as.factor(colorVector)
colorVector = factor(colorVector, levels(colorVector)[c(6, 4, 2, 7, 5, 3, 1)])
sigwithColor <- cbind(sigdf, colorVector)

## Cleaning up the names
#testNames <- as.character(sigwithColor$Antibiotic.Resistant.Gene)
#testNames <- gsub("SHV_", "SHV-", testNames)
#testNames <- gsub("TEM_", "TEM-", testNames)
## noARO <- sapply(strsplit(testNames, split="ARO_[0-9]*\\."), "[[",2)
#geneOrgList <- strsplit(testNames, split = "\\|")
#geneOrg <- unlist(lapply(geneOrgList, function(x) ifelse(length(x) > 1, x[length(x)], x)))
#geneOrg <- gsub('__)', '")', geneOrg, fixed = TRUE)
#geneOrg <- gsub('_"', '', geneOrg, fixed = TRUE)
#geneOrg <- gsub('_$', '', geneOrg)
#geneOrg <- gsub('_)', "')", geneOrg, fixed = TRUE)
## geneOrg <- sapply(strsplit(testNames, split="ARO_[0-9]*\\."), "[[", 4)

#justGene <- sapply(strsplit(geneOrg, split="__"), "[[",1)
## That messes up a lot, so it can't be used as casually as I had previously thought
#sigwithColor$justGene <- justGene


## Actual Visualization
## rownames are gene names
## column names are comparisons
## legend describes the color



testplot <- ggplot(sigwithColor, aes(Sites.Compared, Antibiotic.Concentration..ng.L.)) +
  geom_tile(aes(fill = colorVector)) +
  scale_fill_manual(values = levels(sigwithColor$colorVector),
                    #values = c("red", "green", "blue", "yellow", "magenta", "cyan", "black"),
                    name = "Significant\nDifferences\nwith Term(s)",
                    labels = c("Timepoint", "Sample Site", "Waterway", "Waterway and\n Sample Site", "Waterway and\nTimepoint", "Sample Site\nand Timepoint", "Sample Site\nWaterway\nand Timepoint")) +
  xlab("Sites Compared") +
  ylab("Antibiotic") +
  ggtitle("Antibiotics with Significant Differences Between Sites") +
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 8), plot.title = element_text(hjust = 0.5), legend.key.height = unit(3, "line"))

ggsave(paste0("AntibioticssiggenechangesUpdate.pdf"), width =12, height = 6, limitsize=FALSE, plot = testplot)


