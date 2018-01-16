setwd('/Users/James/Desktop/Summer_2017_Waste_Water_Project/Consolidated/')
rm(list=ls())
myT=read.table('shortBREDwMetadataOnlyDeepestReads.txt',sep='\t',header = TRUE,comment.char = "@")
myT <- myT[ myT$Location == "Sugar Creek" | myT$Location == "Mallard Creek" | myT$Location == "Mountain Control" | myT$Location == "Uwharrie Control",]
myT <- myT[ myT$Sample=="UP A" | myT$Sample=="UP B" | myT$Location == "Mountain Control" | myT$Location == "Uwharrie Control",]

library(stringi)

pValuesupRural <- vector()
names <- vector()
index <- 1
indexes <- vector()
absolute_indexes=vector()

for( i in 47:ncol(myT)){
  if( sum( myT[,i] >0 ) > nrow(myT) /4 ) 
  {
    upRural <- factor(ifelse( myT$Sample=="UP A" | myT$Sample=="UP B" , "UP", "Rural"))
    upRural <- factor(upRural, c("UP", "Rural"))
    myLm <- lm( myT[,i] ~ upRural)			
    myAnova <- anova(myLm)
    pValuesupRural[index] <- myAnova$"Pr(>F)"[1]
    names[index] <- names(myT)[i]
    indexes[index] <- index
    absolute_indexes[index]=i
    index <- index + 1
  }
}

dFrame <- data.frame(names, indexes, absolute_indexes,pValuesupRural)
dFrame$pValuesupRuralAdjusted<- p.adjust( dFrame$pValuesupRural, method = "BH")
dFrame <- dFrame [order(dFrame$pValuesupRuralAdjusted),]

write.table(dFrame, file=paste("shortbred_models_Updated_Format_upRural", ".txt",sep=""), sep="\t",row.names=FALSE)

pdf("shortbredPlots_Updated_Format_upRural.pdf")
for (j in 1:nrow(dFrame)){
  par(mfrow=c(2,2),oma = c(0, 3, 2, 0) + 0.1,mar = c(3, 2, 3, 1) + 0.1,mgp = c(3, 1.25, 0))
  absolute_index=dFrame$absolute_indexes[j]
  index=dFrame$indexes[j]
  bug <- myT[,absolute_index]
  myFrame <- data.frame( bug,upRural)
  boxplot(myT[,absolute_index] ~ upRural, main = paste("Source ",format(dFrame$pValuesupRuralAdjusted[j])), names=c("UP","Rural"))
  stripchart(bug ~ upRural,method="jitter",data = myFrame,vertical = TRUE, pch = 20, add=TRUE ) 
  plot.new()
  plot.new()
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
hist(pValuesupRural)
hist(dFrame$pValuesupRuralAdjusted,main='pValuesupRuralAdjusted')

dev.off()