rm(list=ls())
#Read in libraries
library(vegan)
library(beeswarm)

####
#***Install if you dont have packages!
####


####
#1) read in the data
####

#first the actual bacteria
otu=read.csv(file = "OTUs_clean.csv")[,-1]
#then the taxa
tax=read.csv(file="TAX_clean.csv")
#then their treatment group
meta=read.csv(file = "META_clean.csv")

####
#2) subset to specific week
####

week_indx=which(meta$SAMPLEWEEK=="Week 02")

otu_week=otu[week_indx,]
meta_week=meta[week_indx,]

zeroBacIndx=which(colSums(otu_week)==0)
otu_week=otu_week[,-zeroBacIndx]

#####
#
#3) Alpha diversity
#
#####

#three different meassurements
RICH=estimateR(otu_week )[1,]
SHANNON=diversity(x = otu_week,index = "shannon")
SIMPSON=diversity(x = otu_week,index = "simpson")


#plot
beeswarm(RICH~meta_week$OUA)
#test
t.test(RICH~meta_week$OUA)

####
#***repeat plot and test with the two other diversity metrics!
####

#####
#
#4) Beta diversity
#
#####

#visualization

#calculating nMDS
nmds=metaMDS(otu_week)

#plotting scores
#color by group
plot(scores(nmds, display = "sites"), col=(meta_week$OUA+1), pch=16)

#add a visual aid
ordispider(nmds,groups =  meta_week$OUA, col = (meta_week$OUA+1)[2:1])

#test for multivariate difference
adonis2(otu_week~meta_week$OUA)



#####
#
#3) Investigating individual genera
#
#####
i=1
testDF=data.frame(taxa=colnames(otu_week),pvals=-1)
for(i in 1:NCOL(otu_week)) {
  TTEST=t.test(log10(otu_week[,i]+1)~meta_week$OUA)
  testDF$pvals[i]=TTEST$p.value
  
}

testDF=testDF[order(testDF$pvals),]

head(testDF)

#plot
beeswarm(otu_week$Bacteria_Desulfobacterota_Desulfovibrionia_Desulfovibrionales_Desulfovibrionaceae_Desulfovibrio_NA~meta_week$OUA)
#test
t.test(otu_week$Bacteria_Desulfobacterota_Desulfovibrionia_Desulfovibrionales_Desulfovibrionaceae_Desulfovibrio_NA~meta_week$OUA)

#plot with log transform
beeswarm(log10(1+otu_week$Bacteria_Desulfobacterota_Desulfovibrionia_Desulfovibrionales_Desulfovibrionaceae_Desulfovibrio_NA)~meta_week$OUA)
#test
t.test(log10(1+otu_week$Bacteria_Desulfobacterota_Desulfovibrionia_Desulfovibrionales_Desulfovibrionaceae_Desulfovibrio_NA)~meta_week$OUA)

####
#***repeat plot and test with other interesting bacteria!
####

