library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)

#load mtdna copynumber
mt.copyn<-read.table("mtdna_copyn_allsamples_06012018.txt",header=T,sep="\t")

#remove all samples with mtcopyn less than 250
mt.copyn<-mt.copyn[which(mt.copyn$med.mtcopyn>250),]

#read pop file
pop.file<-read.table("../../../Data_tables/pop_1kg.txt",header=T)

#merge global ancestry with mtcopyn info
mt.copyn.anc<-join(mt.copyn[,c(1,3)],pop.file,by="IID")
mt.copyn.anc$pop<-factor(mt.copyn.anc$pop,levels=c("CEU","YRI","PEL","ACB","ASW","CLM","MXL","PUR"))

plt.mtcopy<-ggplot(mt.copyn.anc)+
  geom_boxplot(aes(pop,med.mtcopyn),outlier.shape = NA)+
  geom_point(aes(pop,med.mtcopyn),alpha=0.5,position=position_jitter(height=0))+
  geom_hline(yintercept=905.29,color="red",linetype="dashed",size=2)+
  theme_bw()+
  theme(axis.title=element_text(size=16))+
  labs(x="Population",y="mtDNA copy number")

ggsave("mtcopy_population.pdf",plt.mtcopy,height=10,width=10)
  
