library(ggpubr)
library(ggplot2)
library(Hmisc)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(purrr)
library(readxl)

sheets <- excel_sheets('InputPCOA_MicrobiomeSource.xlsx')
sample_sheets <- sheets[grepl("Sheet1|Sheet2|Sheet3|Sheet4|Sheet5|Sheet6|Sheet7|Sheet8|Sheet9", sheets)]
list_all <- lapply(sample_sheets , function(x) read_excel(path = 'InputPCOA_MicrobiomeSource.xlsx', sheet = x))
merged=Reduce(function(...)merge(...,by="Species",all=T),list_all)
merged[is.na(merged)] <- 0
input3<-merged[,-1]
rownames(input3)<-merged$Species
otu=otu_table(input3,taxa_are_rows = TRUE)
metadata=import_qiime_sample_data("MetadataPCOA_MicrobiomeSource.txt")
ps=merge_phyloseq(otu,metadata)

ps2<-subset_samples(ps,Specific %in% c("Environmental_control","stool","Ancient","subgingival_plaque","supragingival_plaque","Soil","AncientPublished2","AncientPublished","Negative","Soil_Italy","right_retroauricular_crease","left_retroauricular_crease","Neanderthal"))
bray=ordinate(ps2,method="PCoA",distance="bray")

theme_set(theme_pubr())
p12<-plot_ordination(ps2,bray,color = "Description",axes=1:2)+geom_point(size=7,alpha=0.75)+theme_classic()+
geom_hline(yintercept=0,linetype="dashed",color="black")+geom_vline(xintercept=0,linetype="dashed",color="black")
p23<-plot_ordination(ps2,bray,color = "Description",axes=2:3)+geom_point(size=7,alpha=0.75)+theme_classic()+
geom_hline(yintercept=0,linetype="dashed",color="black")+geom_vline(xintercept=0,linetype="dashed",color="black")
figure <- ggarrange(p12,p23,ncol = 1, nrow = 2)
jpeg("PCoA_MicrobiomeSource.jpg",width=2500,height=2500,res=330)  
figure
dev.off()

metadata2 <- read.table("MetadataSourceTracker.txt",h=T,row.names=1,check=F,comment='') #a new metadata file with SourceSink column
otus2=t(otu)
alpha1 <- alpha2 <- 0.001
otus2<-ceiling(otus2)
train.ix <- which(metadata2$SourceSink=='source')
test.ix <- which(metadata2$SourceSink=='sink')
envs <- metadata2$Description
source("../Script/SourceTracker/sourcetracker-master/src/SourceTracker.r")
st <- sourcetracker(otus2[train.ix,], envs[train.ix])
results <- predict(st,otus2[test.ix], alpha1=alpha1, alpha2=alpha2)
