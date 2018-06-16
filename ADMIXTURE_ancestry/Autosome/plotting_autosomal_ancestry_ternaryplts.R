library(plyr)
library(ggplot2)
library(ggtern)
library(gridExtra)

#read individual id info for admixture results
fam<-read.table('1kg_nam_flipped_pruned.fam',header=F,stringsAsFactors = F)
fam<-fam[,c(1,2)]
colnames(fam)<-c("FID","IID")
#read file linking IDs to pop groups
pop<-read.table('../../Data_tables/pop_1kg.txt',header=T,sep="\t",stringsAsFactors = F)

#link fam file to pop designations
fam<-join(fam,pop,by="IID")

#fill in pop designations for native americans
nat.index<-which(is.na(fam$pop))
fam$pop[nat.index]<-sapply(fam$IID[nat.index],function(x){unlist(strsplit(x,split="_"))[1]})

#read global ancestry from ADMIXTURE Q3
q3<-read.table('1kg_nam_flipped.3.Q',header=F)
colnames(q3)<-c("C1","C2","C3")
q3$IID<-fam$IID
q3$pop<-fam$pop
#continental grouping (African, European, Native American, or admixed) for plotting
q3$cont_group<-NA
q3$cont_group[which(q3$pop=="YRI")]<-"AFR"
q3$cont_group[which(q3$pop=="CEU")]<-"EUR"
q3$cont_group[which(q3$pop%in%c("MAYAN","AYMARAN","QUECHUAN","NAHUAN"))]<-"NAT"
q3$cont_group[which(is.na(q3$cont_group)=="TRUE")]<-"Admixed"

#read MT haplogroup info
haplo<-read.table('../../Data_tables/AMR_mt.haplo',header=T,sep="\t")
colnames(haplo)[1]<-c("IID")
haplo<-haplo[,c(1,4)]
q3<-join(q3,haplo,by="IID")
q3<-q3[-which(q3$Major_Haplogroup=="M"),]
q3$region_haplogroup<-NA
q3$region_haplogroup[which(q3$Major_Haplogroup%in%c("A","B","C","D"))]<-"NAT"
q3$region_haplogroup[which(q3$Major_Haplogroup%in%c("K","J","U","H","T","V","W","I"))]<-"EUR"
q3$region_haplogroup[which(q3$Major_Haplogroup%in%c("L"))]<-"AFR"


fig_2a<-ggtern()+
  geom_point(data=q3[which(q3$cont_group=="Admixed"),],aes(C1,C2,C3),color="grey",alpha=0.6)+
  geom_point(data=q3[which(q3$cont_group%in%c("NAT","EUR","AFR")),],aes(C1,C2,C3,color=cont_group),size=6,alpha=0.2)+
  theme_bw()+
  theme(legend.text=element_text(size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=16),
        plot.title=element_text(size=20))+
  labs(color="Continental\ngroup",title="A.")

fig_2b<-ggtern()+
  geom_point(data=q3[which(q3$cont_group=="Admixed"),],aes(C1,C2,C3,color=pop),alpha=0.6)+
  theme_bw()+
  theme(legend.text=element_text(size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=16),
        plot.title=element_text(size=20))+
  labs(color="Population",x="EUR",y="AFR",z="NAT",title="B.")

fig_2c<-ggtern()+
  geom_point(data=q3,aes(C1,C2,C3,color=Major_Haplogroup),alpha=0.6)+
  theme_bw()+
  theme(legend.text=element_text(size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=16),
        plot.title=element_text(size=20))+
  labs(color="mtDNA\nhaplogroup",x="EUR",y="AFR",z="NAT",title="C.")

fig_2d<-ggtern()+
  geom_point(data=q3,aes(C1,C2,C3,color=region_haplogroup),alpha=0.6)+
  theme_bw()+
  theme(legend.text=element_text(size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=16),
        plot.title=element_text(size=20))+
  labs(color="mtDNA\nhaplogroup",x="EUR",y="AFR",z="NAT",title="D.")

fig_2.combined<-ggtern::grid.arrange(fig_2a,fig_2b,fig_2c,fig_2d,nrow=2,padding=unit(0.01,"line"))

ggsave("ancestry_ternary_plots.pdf",fig_2.combined,height=13,width=13)
