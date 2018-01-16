setwd('/Users/James/Desktop/Summer_2017_Waste_Water_Project/Consolidated/')
rm(list=ls())
myT=read.table('shortBREDwMetadataOnlyDeepestReads.txt',sep='\t',header = TRUE,comment.char = "@")
library(stringi)

myT=myT[which((myT$Sample=="ATE")|(myT$Sample=="FCE")),]

pValuesLocation <- vector()
pValuesATE_FCE <- vector()
pValuesTimepoint <- vector()
names <- vector()
indexes <- vector()
absolute_indexes=vector()
index=1

for( i in 47:ncol(myT)){
  if( sum( myT[,i] >0 ) > nrow(myT) /4 ) 
  {
    ATE_FCE=factor(myT$Sample,c("ATE","FCE"))
    myLm <- lm( myT[,i] ~ myT$Location + ATE_FCE + factor(myT$Timepoint)  )		
    myAnova <- anova(myLm)
    pValuesLocation[index] <- myAnova$"Pr(>F)"[1]
    pValuesATE_FCE[index] <- myAnova$"Pr(>F)"[2]
    pValuesTimepoint[index] <- myAnova$"Pr(>F)"[3]
    names[index] <- names(myT)[i]
    indexes[index] <- index
    absolute_indexes[index]=i
    index <- index + 1
  }
}

dFrame <- data.frame(names, indexes,absolute_indexes,pValuesATE_FCE, pValuesLocation,pValuesTimepoint)
dFrame$pValuesATE_FCEAdjusted<- p.adjust( dFrame$pValuesATE_FCE, method = "BH")
dFrame$pValuesLocationAdjusted<- p.adjust( dFrame$pValuesLocation, method = "BH")
dFrame$pValuesTimepointAdjusted<- p.adjust( dFrame$pValuesTimepoint, method = "BH")
dFrame <- dFrame [order(dFrame$pValuesATE_FCEAdjusted),]

write.table(dFrame, file=paste("shortbred_models_Updated_Format_ATE_FCE", ".txt",sep=""), sep="\t",row.names=FALSE)

pdf("shortbredPlots_Updated_Format_ATE_FCE.pdf")
for (j in 1:nrow(dFrame)){
  par(mfrow=c(2,2),oma = c(0, 3, 2, 0) + 0.1,mar = c(3, 2, 3, 1) + 0.1,mgp = c(3, 1.25, 0))
  absolute_index=dFrame$absolute_indexes[j]
  index=dFrame$indexes[j]
  bug <- myT[,absolute_index]
  time <- factor(myT$Timepoint)
  location <- factor(myT$Location)
  myFrame <- data.frame( bug, time, location, ATE_FCE )
  boxplot(myT[,absolute_index] ~ ATE_FCE, main = paste("Source ",format(dFrame$pValuesATE_FCEAdjusted[j])), names=c("ATE","FCE"))
  stripchart(bug ~ ATE_FCE,method="jitter",data = myFrame,vertical = TRUE, pch = 20, add=TRUE ) 
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

hist(pValuesATE_FCE)
hist(pValuesLocation)
hist(pValuesTimepoint)
plot.new()
hist(dFrame$pValuesATE_FCEAdjusted,main='pValuesATE_FCEAdjusted')
hist(dFrame$pValuesLocationAdjusted,main='pValuesLocationAdjusted')
hist(dFrame$pValuesTimepointAdjusted,main='pValuesTimepointAdjusted')

dev.off()