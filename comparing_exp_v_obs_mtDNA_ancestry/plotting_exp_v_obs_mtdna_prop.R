#this script plots the expected distribution of females from each ancestral population and the observed proportions of mtDNA

#load necessary libraries
library(ggplot2)
library(plyr)
library(reshape2)

#setwd("/Users/Azaidi/Documents/mtproj_files/mitonuclear_project/mtnuc_organization/sex_biased_admixture")


#load boostrapped m/f contributions
mf.boot<-read.table('../Sex_bias_estimation/sex_spec_ancestry_contribs_bootstap.txt',header=T)
#we only need f contributions - normalize by 0.5
mf.boot<-mf.boot[which(mf.boot$sex=="f"),]
mf.boot$norm.proportion<-mf.boot$norm.proportion*2
dmfest.boot<-ddply(mf.boot,.(anc.pop,admx.pop),summarize,mean.prop=mean(norm.proportion),lower=quantile(norm.proportion,probs=0.025),upper=quantile(norm.proportion,probs=0.975))


#read mtDNA haplogroup info
haplo<-read.table('../Data_tables/AMR_mt.haplo',header=T,sep="\t")
colnames(haplo)[1]<-"IID"

haplo$region_haplogroup<-NA
haplo$region_haplogroup[which(haplo$Major_Haplogroup%in%c("A","B","C","D"))]<-"Native_American"
haplo$region_haplogroup[which(haplo$Major_Haplogroup%in%c("L"))]<-"African"
haplo$region_haplogroup[which(haplo$Major_Haplogroup%in%c("K","J","U","H","T","V","W"))]<-"European"
#two individuals carry an Asian haplogroup (M). Remove from analyses
haplo.red<-haplo[which(is.na(haplo$region_haplogroup)=="FALSE"),]

#read population and sex info for each individual
pop.file<-read.table("../Data_tables//pop_1kg.txt",header=T,sep="\t",stringsAsFactors = F)
haplo.red<-join(haplo.red,pop.file,by="IID")
colnames(haplo.red)[7]<-"admx.pop"
#calculate frequency of mt haplogroup in each population
mtprop.pop<-ddply(haplo.red,.(admx.pop),function(x){African_mt=nrow(x[which(x$region_haplogroup=="African"),])/nrow(x)
nat_mt=nrow(x[which(x$region_haplogroup=="Native_American"),])/nrow(x)
Euro_mt=nrow(x[which(x$region_haplogroup=="European"),])/nrow(x)
return(data.frame(African=African_mt,European=Euro_mt,Native_American=nat_mt))})

#convert to long format for plotting
mmtprop.pop<-melt(mtprop.pop,id.vars = c("admx.pop"))
colnames(mmtprop.pop)[c(2,3)]<-c("anc.pop","mtprop")

mt.plt<-ggplot()+
  geom_boxplot(data=dmfest.boot,aes(x=anc.pop,ymin=lower,ymax=upper,lower=lower,upper=upper,middle=upper,fill=anc.pop),color=NA,alpha=0.5,position=position_dodge(),stat="identity",show.legend = T)+
  geom_errorbar(data=mmtprop.pop,aes(x=anc.pop,ymin=mtprop,ymax=mtprop),stat="identity",show.legend = F)+
  facet_wrap(~admx.pop,scales="free")+
  theme_bw()+
  scale_fill_manual(values=c("#8da0cb","#fc8d62","#66c2a5"))+
  #scale_color_manual(values=c("#8da0cb","#fc8d62","#66c2a5"),guide=F)+
  labs(x="Ancestry",y="Expected/Observed MT haplogroup proportion",fill="Ancestry")+
  theme(axis.text.x = element_blank())

ggsave("mt_exp_vs_obs.pdf",mt.plt,height=5,width=7)
