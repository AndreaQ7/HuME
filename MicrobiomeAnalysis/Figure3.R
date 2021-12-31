#Graphical rappresentation of significant species obtained through Kruskal-Wallis analysis (corrected by Benjamini & Hochberg)

library(ggpubr)
library(reshape2)
library(dplyr)

inputkw<-read.table("InputKWadjust.txt",header=T)

inputkw$period<-gsub("Paleolithic","PA",inputkw$period)
inputkw$period<-gsub("Early_Neolithic","EN",inputkw$period)
inputkw$period<-gsub("Final-Middle_Neolithic","FMN",inputkw$period)
inputkw$period<-gsub("Middle_Neolithic","MN",inputkw$period)
inputkw$period<-gsub("Late_Neolithic","LN",inputkw$period)
inputkw$period<-gsub("Copper_Age","CA",inputkw$period)


inputkw$period <- factor(inputkw$period,levels = c("PA", "EN", "MN","FMN", "LN","CA"))
 
inputkw.melt<-melt(inputkw,value.name="ID")

#Create one plot for each complex

tf<-subset(inputkw.melt,inputkw.melt$variable== c("Tannerella_forsythia"))
pg<-subset(inputkw.melt,inputkw.melt$variable== c("Porphyromonas_gingivalis"))
td<-subset(inputkw.melt,inputkw.melt$variable== c("Treponema_denticola"))
rc<-rbind(tf,pg,td) #red complex assembled


rc.plot<-ggplot(rc,aes(period,ID.1,fill=variable)) +geom_boxplot() + theme_classic()+stat_summary(fun=mean, geom="line", aes(group=variable,colour=variable))+scale_fill_brewer(palette="Reds")+coord_cartesian(ylim = c(0.00,0.03))+scale_color_brewer(palette="Reds")+ylab("Relative Abundance")+xlab("Period")+ggtitle("Red Complex")+ theme(plot.title = element_text(hjust = 0.5))

cg<-subset(inputkw.melt,inputkw.melt$variable== c("Campylobacter_gracilis"))
cr<-subset(inputkw.melt,inputkw.melt$variable== c("Campylobacter_rectus"))
cs<-subset(inputkw.melt,inputkw.melt$variable== c("Campylobacter_showae"))
fn<-subset(inputkw.melt,inputkw.melt$variable== c("Fusobacterium_nucleatum"))
pri<-subset(inputkw.melt,inputkw.melt$variable== c("Prevotella_intermedia"))
sc<-subset(inputkw.melt,inputkw.melt$variable== c("Streptococcus_constellatus"))

oc<-rbind(cg,cr,cs,fn,pri,sc)   #orange complex assembled

oc.plot<-ggplot(oc,aes(period,ID.1,fill=variable)) +geom_boxplot() + theme_classic()+stat_summary(fun=mean, geom="line", aes(group=variable,colour=variable))+scale_fill_brewer(palette="Oranges")+coord_cartesian(ylim = c(0.00,0.01))+scale_color_brewer(palette="Oranges")+ylab("Relative Abundance")+xlab("Period")+ggtitle("Orange Complex")+ theme(plot.title = element_text(hjust = 0.5))

sm<-subset(inputkw.melt,inputkw.melt$variable== c("Streptococcus_mitis"))
so<-subset(inputkw.melt,inputkw.melt$variable== c("Streptococcus_oralis"))
ss<-subset(inputkw.melt,inputkw.melt$variable== c("Streptococcus_sanguinis"))
sg<-subset(inputkw.melt,inputkw.melt$variable== c("Streptococcus_gordonii"))


yc<-rbind(sm,so,ss,sg)  #yellow complex assembled
yc.plot<-ggplot(yc,aes(period,ID.1,fill=variable)) +geom_boxplot() + theme_classic()+stat_summary(fun=mean, geom="line", aes(group=variable,colour=variable))+scale_fill_brewer(palette="YlOrRd")+coord_cartesian(ylim = c(0.00,0.2))+scale_color_brewer(palette="YlOrRd")+ylab("Relative Abundance")+xlab("Period")+ggtitle("Yellow Complex")+ theme(plot.title = element_text(hjust = 0.5))

cg<-subset(inputkw.melt,inputkw.melt$variable== c("Capnocytophaga_gingivalis"))
ce<-subset(inputkw.melt,inputkw.melt$variable== c("Capnocytophaga_endodontalis"))
ch<-subset(inputkw.melt,inputkw.melt$variable== c("Capnocytophaga_haemolytica"))
co<-subset(inputkw.melt,inputkw.melt$variable== c("Capnocytophaga_ochracea"))
ec<-subset(inputkw.melt,inputkw.melt$variable== c("Eikenella_corrodens"))

gc<-rbind(cg,ce,ch,co,ec)   #green complex assembled
gc.plot<-ggplot(gc,aes(period,ID.1,fill=variable)) +geom_boxplot() + theme_classic()+stat_summary(fun=mean, geom="line", aes(group=variable,colour=variable))+scale_fill_brewer(palette="Greens")+coord_cartesian(ylim = c(0.00,0.01))+scale_color_brewer(palette="Greens")+ylab("Relative Abundance")+xlab("Period")+ggtitle("Green Complex")+ theme(plot.title = element_text(hjust = 0.5))

ah<-subset(inputkw.melt,inputkw.melt$variable== c("Actinomyces_howellii"))
ao<-subset(inputkw.melt,inputkw.melt$variable== c("Actinomyces_oris"))
an<-subset(inputkw.melt,inputkw.melt$variable== c("Actinomyces_naeslundii"))
av<-subset(inputkw.melt,inputkw.melt$variable== c("Actinomyces_viscosus"))
a171<-subset(inputkw.melt,inputkw.melt$variable== c("Actinomyces_sp._oral_taxon_171"))
a169<-subset(inputkw.melt,inputkw.melt$variable== c("Actinomyces_sp._oral_taxon_169"))
ai<-subset(inputkw.melt,inputkw.melt$variable== c("Actinomyces_israelii"))
a414<-subset(inputkw.melt,inputkw.melt$variable== c("Actinomyces_sp._oral_taxon_414"))
asi<-subset(inputkw.melt,inputkw.melt$variable== c("Actinomyces_slackii"))

ac<-rbind(ah,ao,an,av,a171,a169,ai,a414,asi)
ac2<-ac %>% mutate(group = if_else(variable == "Actinomyces_sp._oral_taxon_414", "gr1", "gr2"))     #Purple complex assembled

ac.plot<- ggplot(ac2,aes(period,ID.1,fill=variable)) +geom_boxplot() + theme_classic()+stat_summary(fun=mean, geom="line", aes(group=variable,colour=variable))+scale_fill_brewer(palette="Blues")+scale_color_brewer(palette="Blues")+facet_wrap(~group,ncol=1,scales = "free")+theme(strip.background = element_blank(),strip.text.x = element_blank())+ylab("Relative Abundance")+xlab("Period")+ggtitle("Purple Complex")+ theme(plot.title = element_text(hjust = 0.5))

theme_set(theme_pubr())
jpeg("ComplexNEW.jpg",width = 5000,height = 3800, res=300)
ggarrange(rc.plot, oc.plot, yc.plot,gc.plot,ac.plot,labels = c("A", "B", "C","D","E"), ncol = 2, nrow = 3)
dev.off()
