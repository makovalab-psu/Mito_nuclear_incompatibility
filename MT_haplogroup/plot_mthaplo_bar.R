#plot breakdown of hpalogroup from each region
library(ggplot2)
library(plyr)
library(reshape2)

#setwd("~/Documents/mtproj_files/mitonuclear_project/mtnuc_organization/MT_haplogroup/")
#read mtDNA haplogroup info
haplo<-read.table("../Data_tables/AMR_mt.haplo",header=T,sep="\t")
colnames(haplo)[1]<-"IID"
#two individuals carry an Asian haplogroup (M). Remove from analyses
haplo.red<-haplo[-which(haplo$Major_Haplogroup=="M"),]
haplo.red<-droplevels(haplo.red)
pop<-read.table("../Data_tables//pop_1kg.txt",header=T,sep="\t")
haplo.red<-join(haplo.red,pop,by="IID")
dhaplo.red<-ddply(haplo.red,.(pop),function(x){table(x$Major_Haplogroup)/sum(table(x$Major_Haplogroup))})
mdhaplo.red<-melt(dhaplo.red,id.vars=c("pop"))
colnames(mdhaplo.red)<-c("pop","major_haplogroup","freq")
mdhaplo.red$region_haplogroup<-NA
mdhaplo.red$region_haplogroup[which(mdhaplo.red$major_haplogroup%in%c("A","B","C","D"))]<-"Native American"
mdhaplo.red$region_haplogroup[which(mdhaplo.red$major_haplogroup%in%c("L"))]<-"African"
mdhaplo.red$region_haplogroup[which(mdhaplo.red$major_haplogroup%in%c("K","J","U","H","T","V","W","I"))]<-"European"

col1<-c("#238b45","#41ae76","#66c2a4","#99d8c9","#8c2d04","#d94801","#f16913","#fd8d3c","#fdae6b","#fdd0a2","#feedde","#fff5eb","#8c6bb1")
mdhaplo.red$major_haplogroup<-factor(mdhaplo.red$major_haplogroup,levels=c("A","B","C","D","H","I","J","K","T","U","V","W","L"))

mthaplo.bar<-ggplot()+
  geom_bar(data=mdhaplo.red,aes(region_haplogroup,freq,fill=major_haplogroup),stat="identity")+
  facet_wrap(~pop)+
  scale_fill_manual(values=col1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))+
  labs(x="Admixed population",y="Haplogroup frequency",fill="MT haplogroup")

ggsave("mt_haplo_barplot.pdf",mthaplo.bar,width=7,height=5)
