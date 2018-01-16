setwd('/Users/James/Desktop/Summer_2017_Waste_Water_Project/Consolidated/')
rm(list=ls())
myT=read.table('shortBREDwMetadataOnlyDeepestReads.txt',sep='\t',header = TRUE,comment.char = "@")
library(stringi)

myT=myT[which((myT$Sample=="PCE")|(myT$Sample=="ATE")),]

pValuesLocation <- vector()
pValuesPCE_ATE <- vector()
pValuesTimepoint <- vector()
names <- vector()
indexes <- vector()
absolute_indexes=vector()
index=1

for( i in 47:ncol(myT)){
  if( sum( myT[,i] >0 ) > nrow(myT) /4 ) 
  {
    PCE_ATE=factor(myT$Sample,c("PCE","ATE"))
    myLm <- lm( myT[,i] ~ myT$Location + PCE_ATE + factor(myT$Timepoint)  )		
    myAnova <- anova(myLm)
    pValuesLocation[index] <- myAnova$"Pr(>F)"[1]
    pValuesPCE_ATE[index] <- myAnova$"Pr(>F)"[2]
    pValuesTimepoint[index] <- myAnova$"Pr(>F)"[3]
    names[index] <- names(myT)[i]
    indexes[index] <- index
    absolute_indexes[index]=i
    index <- index + 1
  }
}

dFrame <- data.frame(names, indexes,absolute_indexes,pValuesPCE_ATE, pValuesLocation,pValuesTimepoint)
dFrame$pValuesPCE_ATEAdjusted<- p.adjust( dFrame$pValuesPCE_ATE, method = "BH")
dFrame$pValuesLocationAdjusted<- p.adjust( dFrame$pValuesLocation, method = "BH")
dFrame$pValuesTimepointAdjusted<- p.adjust( dFrame$pValuesTimepoint, method = "BH")
dFrame <- dFrame [order(dFrame$pValuesPCE_ATEAdjusted),]

write.table(dFrame, file=paste("shortbred_models_Updated_Format_PCE_ATE", ".txt",sep=""), sep="\t",row.names=FALSE)

pdf("shortbredPlots_Updated_Format_PCE_ATE.pdf")
for (j in 1:nrow(dFrame)){
  par(mfrow=c(2,2),oma = c(0, 3, 2, 0) + 0.1,mar = c(3, 2, 3, 1) + 0.1,mgp = c(3, 1.25, 0))
  absolute_index=dFrame$absolute_indexes[j]
  index=dFrame$indexes[j]
  bug <- myT[,absolute_index]
  time <- factor(myT$Timepoint)
  location <- factor(myT$Location)
  myFrame <- data.frame( bug, time, location, PCE_ATE )
  boxplot(myT[,absolute_index] ~ PCE_ATE, main = paste("Source ",format(dFrame$pValuesPCE_ATEAdjusted[j])), names=c("PCE","ATE"))
  stripchart(bug ~ PCE_ATE,method="jitter",data = myFrame,vertical = TRUE, pch = 20, add=TRUE ) 
  boxplot(myT[,absolute_index] ~ location, main= paste( "location ", format(dFrame$pValuesLocationAdjusted[j])),ylab="")
  stripchart(bug ~ location,method="jitter", data = myFrame,vertical = TRUE, pch = 20, add=TRUE ) 
  boxplot(myT[,absolute_index] ~ time, main= paste( "time ", format(dFrame$pValuesTimepointAdjusted[j])),ylab = "",xaxt="n")
  stripchart(bug ~ time,method="jitter",data = myFrame,vertical = TRUE, pch = 20, add=TRUE ) 
  axis(1, at=c(1,2,3,4), labels=c("Late Winter","Early Spring","Mid Spring","Mid Summer"),cex.axis=.6) 
  plot.new()
  if(grepl("gi.",names[index])){
    a=strsplit(strsplit(names[index],split="\\__")[[1]][1],split="\\.")
    aLabel=paste("Gene Accession:",a[[1]][2])
    if(substr(aLabel,nchar(aLabel)-1,nchar(aLabel)-1)=="_"){
    aLabel=stri_replace_last_regex(str=aLabel,pattern="_",replacement = "\\.")
    } 
    bLabel=paste("Gene Name:",a[[1]][4])
    bLabel=stri_replace_all_regex(str=bLabel,pattern = "_",replacement=" ")
    cLabel=stri_split_fixed(str = names[index],pattern = "__",n=2)[[1]][2]
    cLabel=stri_replace_all_regex(str=cLabel,pattern="str__",replacement="strain ")
    cLabel=stri_replace_all_regex(str=cLabel,pattern="__",replacement = " ")
    cLabel=stri_sub(cLabel, 1, -2)
    cLabel=stri_replace_all_regex(str=cLabel,pattern="_",replacement = " ")
    cLabel=stri_replace_all_regex(str=cLabel,pattern="\\.",replacement = "-")
    cLabel=paste("Source Organism:",cLabel)
    mtext(text = c(aLabel,bLabel,cLabel),side=3,line=c(1,0,-1),outer = TRUE,cex=.8)
    }
  else{
    aLabel = names[index]
    mtext(text = aLabel,side=3,line=0,outer = TRUE)
  }
  mtext(text = "ShortBred Protein Family Representative Marker Sequence (RPKM)",side=2,line=1,outer = TRUE)
}

hist(pValuesPCE_ATE)
hist(pValuesLocation)
hist(pValuesTimepoint)
plot.new()
hist(dFrame$pValuesPCE_ATEAdjusted,main='pValuesPCE_ATEAdjusted')
hist(dFrame$pValuesLocationAdjusted,main='pValuesLocationAdjusted')
hist(dFrame$pValuesTimepointAdjusted,main='pValuesTimepointAdjusted')

dev.off()