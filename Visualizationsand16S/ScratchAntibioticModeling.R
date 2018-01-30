# Antibiotic Modeling from Scratch

# Visualization of the amount of ND and <LOQ values

rm(list=ls())
## library("phyloseq")
## options(warn=2)
library("ggplot2")
elsewhere <- c("MNT A", "MNT B", "UW A", "UW B", "RURAL")
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
imputeValues <- 2 #-1 for no imputation

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

## dataType <- "qiimeclosedmerge"
## dataType <- "MetaPhlAn"
## dataType <- "HUMANn"
dataType <- "metadata"
isQIIME <- FALSE
if (.Platform$OS.type == "unix") {
  baseDir <- "/Users/mbrown67/Google Drive/Programs/GibasExperiment/UEGP"  
} else {
  baseDir <- "F:/Google Drive/Programs/GibasExperiment/UEGP"}

## I think the modelDescription is largely meaningless now

## helpfulName <- paste0(helpfulName, modelDescription)
metadataDir <- paste(baseDir, "metadata", sep="/")
#processedDir <- paste(dataDir, "processed", sep="/")
analysisDir <- paste(baseDir, "results", sep="/")

## Load data
setwd(baseDir)
fileName <- "MetaMasterAllinOneWGSSfixShortBREDfixHUMANn.csv"
myT <- read.csv(fileName, header=TRUE, comment.char="@", stringsAsFactors = FALSE)

## Dropping the duplicates
myT <- myT[myT$Lane != 7,]


## This now depends on the datatype
## This will change when the metadata gets reprocessed

#myColClasses <- c(rep("character", 5), rep("numeric", 5), "character", rep("numeric", 3), rep("character", 24), "numeric", rep("character", 10))

## ShortBRED there are 46
## MetaPhlAn there are 47
#myT <- read.csv(fileName, header=TRUE, colClasses = myColClasses, comment.char="@")

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
} else if (imputeValues == 2){
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
} else {
  myT[myT == "ND"] <- -1
  myT[myT == "<LOQ"] <- 0
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

## Need to set the order of the sample levels somewhere.

setwd(analysisDir)

#Antibiotics are 16-25

antibioticX <- myT[,16:25]
antibioticX[antibioticX == "ND"] <- -1
antibioticX[antibioticX == "<LOQ"] <- -10
antibioticX <- apply(antibioticX, 2, as.numeric)

#antibioticCol <- 19
#locVec <- runCompare[[2]]
#ggplot(subset(myT, Sample %in% locVec)) +
#  geom_jitter(aes(Timepoint, as.numeric(get(colnames(myT)[antibioticCol])), color=Sample, shape=Location), width=0.25, height=0) +
#  ylab(colnames(myT)[antibioticCol]) #+
  #scale_color_manual(values = c(runCompare[[2]][1] = "red", runCompare[[2]][2] = "blue", "ND" = "black"))

for(locVec in runCompare) {
  for(antibioticCol in c(18:25)){
    shouldRun <- 1
    subT <- myT
    
    if (any(locVec %in% elsewhere) == FALSE) {
      subT <- subT[subT$Location == "Mallard Creek" | subT$Location == "Sugar Creek",]
    }
    
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
    subT <- subT[subT$Sample %in% locVec,]
    acPlot <- ggplot(subT) +
      geom_jitter(aes(Timepoint, as.numeric(get(colnames(subT)[antibioticCol])), color=Sample, shape=Location), width=0.25, height=0) +
      ylab(colnames(subT)[antibioticCol]) #+
    #scale_color_manual(name = "Missing", values = c('-1' = "black", 
    #                          '-10' = "black"))
    ggsave(paste0(strsplit(colnames(subT)[antibioticCol], split=1)[1], "_", paste0(locVec, collapse="_"), ".pdf"), plot = acPlot)
  }
}
#     subT <- subT[complete.cases(subT[,antibioticCol]),]
#     for (locItem in locVec) {
#       if(length(subT[subT$Sample == locItem,]$Sample) == 0) {
#         shouldRun <- 0
#       }
#     }
#     
#     if (shouldRun == 0) {
#       next
#     }
#     
#     subT$Location <- as.factor(subT$Location)
#     subT$Sample <- as.factor(subT$Sample)
#     location <- subT$Location
#     sample <- subT$Sample
#     timepoint <- subT$Timepoint
#     tempdf <- cbind(as.numeric(subT[,antibioticCol]), location, sample, timepoint)
#     if (all.equal(c("UP", "RURAL"), locVec) == TRUE) {
#       myLm <- lm(as.numeric(subT[,antibioticCol]) ~ sample)
#       writeVals <- cbind(colnames(subT)[antibioticCol], t(as.data.frame(anova(myLm)$"Pr(>F)"[1])))
#       write.table(writeVals, file=paste0("Impute_", imputeValues, "_", paste0(locVec, collapse="_"),"_",strsplit(colnames(subT)[antibioticCol], split="\\.")[[1]][1],".txt"), sep="\t", col.names = c("Antibiotic", rownames(anova(myLm))[1]), row.names = FALSE)
#       ## print(paste0(locVec, antibioticCol))
#       pdf(paste0("Impute_", imputeValues, "_", paste0(locVec, collapse="_"),"_", strsplit(colnames(subT)[antibioticCol], split="\\.")[[1]][1], ".pdf"))
#       par(mfrow=c(2,2),
#           oma = c(0, 3, 2, 0) + 0.1,
#           mar = c(3, 2, 3, 1) + 0.1,
#           mgp = c(3, 1.25, 0))
#       
#       boxplot(tempdf[,1]~sample, data = tempdf, main=paste0("Location ", format(anova(myLm)$"Pr(>F)"[1], digits=3)), xaxt="n", outline=FALSE)
#       axis(1, 1:2, levels(sample))
#       stripchart(tempdf[,1]~sample, vertical = TRUE, method="jitter", pch = 16, add=TRUE)
#       mtext(paste("Impute_", imputeValues, "_", paste(locVec, collapse="_"), strsplit(colnames(subT)[antibioticCol], split="\\.")[[1]][1], sep = " "), outer=TRUE, line=0)
#       mtext("Concentration (ng/L)", outer=TRUE, side=2, line=1)
#       plot.new()
#       dev.off()
#     } else {
#       myLm <- lm(as.numeric(subT[,antibioticCol]) ~ location + sample + timepoint)
#       
#       writeVals <- cbind(colnames(subT)[antibioticCol], t(as.data.frame(anova(myLm)$"Pr(>F)"[1:3])))
#       write.table(writeVals, file=paste0("Impute_", imputeValues, "_", paste0(locVec, collapse="_"),"_",strsplit(colnames(subT)[antibioticCol], split="\\.")[[1]][1],".txt"), sep="\t", col.names = c("Antibiotic", rownames(anova(myLm))[1:3]), row.names = FALSE)
#       ## print(paste0(locVec, antibioticCol))
#       pdf(paste0("Impute_", imputeValues, "_", paste0(locVec, collapse="_"),"_", strsplit(colnames(subT)[antibioticCol], split="\\.")[[1]][1], ".pdf"))
#       par(mfrow=c(2,2),
#           oma = c(0, 3, 2, 0) + 0.1,
#           mar = c(3, 2, 3, 1) + 0.1,
#           mgp = c(3, 1.25, 0))
#       
#       boxplot(tempdf[,1]~location, data = tempdf, main=paste0("Location ", format(anova(myLm)$"Pr(>F)"[1], digits=3)), xaxt="n", outline=FALSE)
#       axis(1, 1:2, levels(location))
#       stripchart(tempdf[,1]~location, vertical = TRUE, method="jitter", pch = 16, add=TRUE)
#       boxplot(tempdf[,1]~sample, data = tempdf, main=paste0("Sample ", format(anova(myLm)$"Pr(>F)"[2], digits=3)), xaxt="n", outline=FALSE)
#       axis(1, 1:2, levels(sample))
#       stripchart(tempdf[,1]~sample, vertical = TRUE, method="jitter", pch = 16, add=TRUE)
#       boxplot(tempdf[,1]~timepoint, data = tempdf, main=paste0("Timepoint ", format(anova(myLm)$"Pr(>F)"[3], digits=3)), xaxt="n", outline=FALSE)
#       axis(1, 1:4, levels(timepoint))
#       stripchart(tempdf[,1]~timepoint, vertical = TRUE, method="jitter", pch = 16, add=TRUE)
#       mtext(paste("Impute_", imputeValues, "_", paste(locVec, collapse="_"), strsplit(colnames(subT)[antibioticCol], split="\\.")[[1]][1], sep = " "), outer=TRUE, line=0)
#       mtext("Concentration (ng/L)", outer=TRUE, side=2, line=1)
#       plot.new()
#       dev.off()
#     }
#   }
# }