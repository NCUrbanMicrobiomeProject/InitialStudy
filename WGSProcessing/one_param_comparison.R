rm(list=ls())
setwd("/Users/anju/Desktop/UEGP_stats")
myT <- read.table("/Users/anju/Desktop/UEGP_stats/MetaPhlAn_LogNormwithMetadata_species1.txt", header=TRUE, sep="\t")
### Replacing dots with an unscore in column headers:
#names(myT)<-gsub(".", "_", colnames(myT), fixed = TRUE)
library(stringr)
myT$Timepoint<-as.character(myT$Timepoint)
myT$Timepoint<-str_replace(myT$Timepoint,"1","Late Winter")
myT$Timepoint<-str_replace(myT$Timepoint,"2","Early Spring")
myT$Timepoint<-str_replace(myT$Timepoint,"3","Mid Spring")
myT$Timepoint<-str_replace(myT$Timepoint,"4","Mid Summer")

### Selecting the replicate with highest depth (Mathew's snippet):

savemyT <- myT[FALSE,]
## Iterates over the unique Sample.IDs
for (eachVal in unique(myT$Sample.ID)) {
  ## Selects just those samples
  myTperSample <- myT[myT$Sample.ID == eachVal,]
  ## Grabs the row with the maximum value for depth
  rep<-myTperSample[which(as.numeric(myTperSample$InputReads) == max(as.numeric(myTperSample$InputReads))),]
  ## Builds the "save" dataframe one row at a time
  savemyT <- rbind(savemyT, rep)
}
## Overwrite my T
myT <- savemyT

x<-"HOSP"
y<-"RES"

###Simple linear model (Dr Fodor's code):
myT <- myT[ myT$Location == "Sugar Creek" | myT$Location == "Mallard Creek" ,]
myT <- myT[ myT$Sample==x | myT$Sample==y ,]
pValuesLocation <- vector()
pValuesSource <- vector()
pValuesTimepoint <- vector()
names <- vector()
index <- 1
indexes <- vector()
pdf(paste("metaphlan_species_Plots_" , x , y, ".pdf", sep=""))
par(mfrow=c(2,2),
    oma = c(0, 3, 2, 0) + 0.1,
    mar = c(4, 3, 3, 1) + 0.1,
    mgp = c(3, 1.5, 0),
    las = 1)

for( i in 48:ncol(myT))
  if( sum( myT[,i] >0 ) > nrow(myT) /4 ) 
  {
    xy <- factor(ifelse( myT$Sample== x, x, y))
    xy<- factor(xy,c(x,y))
    myT$Timepoint<-factor(myT$Timepoint, c("Late Winter", "Early Spring", "Mid Spring", "Mid Summer"))
    
    myT$Location<-factor(myT$Location)
    myLm <- lm( myT[,i] ~ myT$Location + xy + factor(myT$Timepoint)  )		
    myAnova <- anova(myLm)
    pValuesLocation[index] <- myAnova$"Pr(>F)"[1]
    pValuesSource[index] <- myAnova$"Pr(>F)"[2]
    pValuesTimepoint[index] <- myAnova$"Pr(>F)"[3]
    names[index] <- names(myT)[i]
    indexes[index] <- index
    
    bug <- myT[,i]
    time <- myT$Timepoint
    location <- myT$Location
    
    myFrame <- data.frame( bug, time, location, xy)
    
    index <- index + 1
    
  }
dFrame <- data.frame(names, indexes, pValuesSource, pValuesLocation, pValuesTimepoint)

dFrame$pValuesSourceAdjusted<- p.adjust( dFrame$pValuesSource, method = "BH")
dFrame$pValuesLocationAdjusted<- p.adjust( dFrame$pValuesLocation, method = "BH")
dFrame$pValuesTimepointAdjusted<- p.adjust( dFrame$pValuesTimepoint, method = "BH")

# Filter the df to select significant adjusted pvalues
dFrame2<-dFrame[dFrame$pValuesSourceAdjusted<0.05 |dFrame$pValuesLocationAdjusted<0.05 | dFrame$pValuesTimepointAdjusted<0.05 ,]
dFrame2 <- dFrame2 [order(dFrame2$pValuesSourceAdjusted),]

# Plots for bugs that are significantly different in the two enviornments

for( i in 48:ncol(myT)) {
  
  for(n in 1:nrow(dFrame2)){
    if(names(myT)[i] == dFrame2$names[n]){
      xy <- factor(ifelse( myT$Sample== x, x, y))
      xy<- factor(xy,c(x,y))
      bug <- myT[,i]
      
      boxplot(  bug ~ xy , main = paste("Source", format(dFrame2$pValuesSourceAdjusted[n],digits=3)), outline = FALSE)
      
      stripchart(bug~ xy, method = "jitter", vertical = TRUE, pch = 16, add=TRUE )		
      plot(  bug ~ location, main= paste( "Location", format(dFrame2$pValuesLocationAdjusted[n],digits=3)), outline = FALSE, ylab ="", xlab="")
      
      stripchart(bug~ location, method = "jitter", vertical = TRUE, pch = 16, add=TRUE )		
      
      
      plot(bug  ~ time, main= paste("Timepoint", format(dFrame2$pValuesTimepointAdjusted[n],digits=3)), xaxt="n", outline = FALSE,ylab ="", xlab="" )
      axis(side = 1, at = time , labels = gsub(" ", "\n", time))
      
      stripchart(bug~ time, method = "jitter", vertical = TRUE, pch = 16, add=TRUE )
      plot(1, type="n", axes=F, xlab="", ylab="")  
      
      mtext(names(myT)[i], outer=TRUE, cex = 1)
      mtext("Metaphlan2 Relative Abundance", outer=TRUE, side=2, line=1,las=3,  cex = 1)
    } 
  }
}  


write.table(dFrame, file=paste("Metaphlan_species_models_", x, y, ".txt", sep=""), sep="\t",row.names=FALSE)

hist(dFrame$pValuesSource)
hist(dFrame$pValuesLocation)
hist(dFrame$pValuesTimepoint)
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5,paste( "Date:", Sys.Date(), "R_version:", R.version$minor))

dev.off()