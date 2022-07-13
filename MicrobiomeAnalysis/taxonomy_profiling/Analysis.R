library(phyloseq)
library(ggplot2)
library(dplyr)
library(readxl)
library(microViz)

input<-read.table("B-out.txt")
input3<-input[,-1]
rownames(input3)<-input$Taxa
otu=otu_table(input3,taxa_are_rows = TRUE)
metadata <- read.delim("MappingPhyloseq.txt", header=TRUE, sep = "\t")
metadatadf <- data.frame(metadata)
sample <- sample_data(metadatadf)
rownames(sample)<-samples_df$SampleID
tax<-tax_table(as.matrix(rownames(otu)))
rownames(tax)<-rownames(otu)
ps=merge_phyloseq(otu,sample,tax)
id.contam<-c("P-T3S","GIL12","Fer.T65","Fer.T69")
ps<-subset_samples(psprova,!Name %in% id.contam)
ps = filter_taxa(ps, function(x) sum(x > 3) > (0.02*length(x)), TRUE) 
#Cluster
ps %>%
    tax_transform(rank = "unique", trans = "identity") %>%
    dist_calc(dist = "aitchison") %>%
    ord_calc(
        method = "PCoA"
    ) %>% 
    ord_plot(
        axes = c(1, 2),
        colour = "Cluster", fill = "Cluster",
        shape = "Period", alpha = 0.5,
        size = 5
    ) + 
    
    stat_chull(
        ggplot2::aes(colour = Cluster)
    ) + scale_colour_manual(values=c("#00FF66","#0066FF","#CC00FF","#CCFF00","#FF0000"))
    
#Gradient reads
pp<-ps %>%
    tax_transform(rank = "unique", trans = "identity") %>%
    dist_calc(dist = "aitchison") %>%
    ord_calc(
        method = "PCoA"
    ) %>%  
    ord_plot(
        axes = c(1, 2),
        fill = "Cleaned.reads",
        alpha = 0.5,
        size = 5)

pp1<-pp$data
ggplot(pp1,aes(MDS1,MDS2))+geom_point(aes(colour=as.numeric(Cleaned.reads)),size=5)+scale_colour_gradientn(colours = terrain.colors(10))+theme_minimal()
    
#Extraction batch
ps %>%
    tax_transform(rank = "unique", trans = "identity") %>%
    dist_calc(dist = "aitchison") %>%
    ord_calc(
        method = "PCoA"
    ) %>% 
    ord_plot(
        axes = c(1, 2),
        colour = "Extraction_batch", fill = "Extraction_batch",
        alpha = 0.5,
        size = 5)
        
#Library_batch
ps %>%
    tax_transform(rank = "unique", trans = "identity") %>%
    dist_calc(dist = "aitchison") %>%
    ord_calc(
        method = "PCoA"
    ) %>% 
    ord_plot(
        axes = c(1, 2),
        colour = "Library_batch", fill = "Library_batch",
        alpha = 0.5,
        size = 5)
  
    
#Network
library(NetCoMi)    
    
net_aitchison <- netConstruct(psprova2,
                              measure = "aitchison",
                              zeroMethod = "multRepl",
                              sparsMethod = "knn", 
                              kNeighbor = 6,
                              verbose = 3)

props_aitchison <- netAnalyze(net_aitchison,
                              clustMethod = "hierarchical",
                              clustPar = list(method = "average", k = 5),
                              hubPar = "eigenvector") 
plot(props_aitchison, 
    nodeColor = "cluster", 
    nodeSize = "fix",
    hubTransp = 40,
    edgeTranspLow = 60,
    charToRm = "00000",
    labelLength = 6L,
    mar = c(3, 1, 1, 1),cexNodes = 2)
    
# Plot Significant species
library(dplyr)
ps1=transform_sample_counts(ps,function(x) 100*x/sum(x))

tb <- psmelt(ps2) %>%
    group_by(Cluster,Name,ta1)
    
tb$Period <- factor(tb$Period,levels = c("PA", "EN", "MN","FMN", "LN","CA"))

tf<-subset(tb,tb$OTU== c("Tannerella_forsythia"))
pg<-subset(tb,tb$OTU== c("Porphyromonas_gingivalis"))
td<-subset(tb,tb$OTU== c("Treponema_denticola"))
rc<-rbind(tf,pg,td)



cg<-subset(tb,tb$OTU== c("Campylobacter_gracilis"))
cr<-subset(tb,tb$OTU== c("Campylobacter_rectus"))
cs<-subset(tb,tb$OTU== c("Campylobacter_showae"))
pri<-subset(tb,tb$OTU== c("Prevotella_intermedia"))

oc<-rbind(cg,cr,cs,pri)


sm<-subset(tb,tb$OTU== c("Streptococcus_mitis"))
so<-subset(tb,tb$OTU== c("Streptococcus_oralis"))
ss<-subset(tb,tb$OTU== c("Streptococcus_sanguinis"))
sg<-subset(tb,tb$OTU== c("Streptococcus_gordonii"))
yc<-rbind(sm,so,ss,sg)


cg<-subset(tb,tb$OTU== c("Capnocytophaga_sputigena"))
ce<-subset(tb,tb$OTU== c("Capnocytophaga_endodontalis"))
ch<-subset(tb,tb$OTU== c("Capnocytophaga_haemolytica"))
co<-subset(tb,tb$OTU== c("Capnocytophaga_ochracea"))
ec<-subset(tb,tb$OTU== c("Eikenella_corrodens"))

gc<-rbind(cg,ce,ch,co,ec)


ah<-subset(tb,tb$OTU== c("Actinomyces_howellii"))
ao<-subset(tb,tb$OTU== c("Actinomyces_oris"))
an<-subset(tb,tb$OTU== c("Actinomyces_naeslundii"))
ai<-subset(tb,tb$OTU== c("Actinomyces_israelii"))
a414<-subset(tb,tb$OTU== c("Actinomyces_sp._oral_taxon_414"))
asi<-subset(tb,tb$OTU== c("Actinomyces_slackii"))
av<-subset(tb,tb$OTU== c("Actinomyces_viscosus"))
a171<-subset(tb,tb$OTU== c("Actinomyces_sp._oral_taxon_171"))
a169<-subset(tb,tb$OTU== c("Actinomyces_sp._oral_taxon_169"))

ac<-rbind(ah,ao,an,av,a171,a169,ai,a414,asi)
ac2<-ac %>% mutate(group = if_else(OTU == "Actinomyces_sp._oral_taxon_414", "gr1", "gr2"))


#OTHERS
ad<-subset(tb,tb$OTU== c("Abiotrophia_defectiva"))
aa<-subset(tb,tb$OTU== c("Aggregatibacter_aphrophilus"))
bp<-subset(tb,tb$OTU== c("Burkholderia_pseudomallei"))
co<-subset(tb,tb$OTU== c("Capnocytophaga_sp._oral_taxon_323"))
do<-subset(tb,tb$OTU== c("Desulfomicrobium_orale"))
fa<-subset(tb,tb$OTU== c("Filifactor_alocis"))
kp<-subset(tb,tb$OTU== c("Klebsiella_pneumoniae"))
lmi<-subset(tb,tb$OTU== c("Lautropia_mirabilis"))
ne<-subset(tb,tb$OTU== c("Neisseria_elongata"))
oo<-subset(tb,tb$OTU== c("Olsenella_sp._oral_taxon_807"))
ot<-subset(tb,tb$OTU== c("Ottowia_sp._oral_taxon_894"))
pm<-subset(tb,tb$OTU== c("Parvimonas_micra"))
ra<-subset(tb,tb$OTU== c("Rothia_aeria"))
so<-subset(tb,tb$OTU== c("Schaalia_odontolytica"))
sa<-subset(tb,tb$OTU== c("Streptococcus_anginosus"))
ss<-subset(tb,tb$OTU== c("Streptococcus_salivarius"))
sto<-subset(tb,tb$OTU== c("Streptococcus_sp._oral_taxon_061"))
stro<-subset(tb,tb$OTU== c("Streptococcus_sp._oral_taxon_431"))

flora<-rbind(ad,aa,bp,co,do,fa,kp,lmi,ne,oo,ot,pm,ra,so,sa,ss,sto,stro)
other.plot<-ggplot(flora,aes(Period,Abundance)) +geom_boxplot(outlier.shape = NA) + theme_classic()+ylab("Relative Abundance")+xlab("Period")+ggtitle("Other")+ theme(plot.title = element_text(hjust = 0.5))+geom_jitter(shape=16, position=position_jitter(0.2),aes(colour=OTU))+facet_wrap(~OTU, scales = "free")+theme(legend.position = "none")

rc.plot2<-ggplot(rc,aes(Period,Abundance,fill=OTU)) +geom_boxplot(outlier.shape = NA) + theme_classic()+stat_summary(fun=mean, geom="line", aes(group=OTU,colour=OTU))+scale_fill_brewer(palette="Reds")+scale_color_brewer(palette="Reds")+ylab("Relative Abundance")+xlab("Period")+ggtitle("Red Complex")+ theme(plot.title = element_text(hjust = 0.5))+ geom_jitter(shape=16, position=position_jitter(0.2),aes(colour=OTU))

oc.plot2<-ggplot(oc,aes(Period,Abundance,fill=OTU)) +geom_boxplot(outlier.shape = NA) + theme_classic()+stat_summary(fun=mean, geom="line", aes(group=OTU,colour=OTU))+scale_fill_brewer(palette="Oranges")+scale_color_brewer(palette="Oranges")+ylab("Relative Abundance")+xlab("Period")+ggtitle("Orange Complex")+ theme(plot.title = element_text(hjust = 0.5))+ geom_jitter(shape=16, position=position_jitter(0.2),aes(colour=OTU))

yc.plot<-ggplot(yc,aes(Period,Abundance,fill=OTU)) +geom_boxplot(outlier.shape = NA) + theme_classic()+stat_summary(fun=mean, geom="line", aes(group=OTU,colour=OTU))+scale_fill_manual(values=c("#B8860B","#FFD700","#FFFF00","#DAA520"))+scale_color_manual(values=c("#B8860B","#FFD700","#FFFF00","#DAA520"))+ylab("Relative Abundance")+xlab("Period")+ggtitle("Yellow Complex")+ theme(plot.title = element_text(hjust = 0.5))+geom_jitter(shape=16, position=position_jitter(0.2),aes(colour=OTU))

gc.plot<-ggplot(gc,aes(Period,Abundance,fill=OTU)) +geom_boxplot(outlier.shape = NA) + theme_classic()+stat_summary(fun=mean, geom="line", aes(group=OTU,colour=OTU))+scale_fill_brewer(palette="Greens")+scale_color_brewer(palette="Greens")+coord_cartesian(ylim = c(0.00,3))+ylab("Relative Abundance")+xlab("Period")+ggtitle("Green Complex")+ theme(plot.title = element_text(hjust = 0.5))+geom_jitter(shape=16, position=position_jitter(0.2),aes(colour=OTU))

ac.plot<-ggplot(ac2,aes(Period,Abundance,fill=OTU)) +geom_boxplot(outlier.shape = NA) + theme_classic()+stat_summary(fun=mean, geom="line", aes(group=OTU,colour=OTU))+scale_fill_brewer(palette="Blues")+scale_color_brewer(palette="Blues")+facet_wrap(~group,ncol=1,scales = "free")+theme(strip.background = element_blank(),strip.text.x = element_blank())+ylab("Relative Abundance")+xlab("Period")+ggtitle("Purple Complex")+ theme(plot.title = element_text(hjust = 0.5))+geom_jitter(shape=16, position=position_jitter(0.2),aes(colour=OTU))

