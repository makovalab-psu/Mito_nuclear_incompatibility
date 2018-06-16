library(ggplot2)
library(plyr)

#setwd("~/Documents/mtproj_files/mitonuclear_project/mtnuc_organization/Copy_number/copyn_analysis_06012018/")

ceu.copyn<-read.table("../Low_coverage/CEU/CEU.header.mtdna_no",header=T,sep="\t")
yri.copyn<-read.table("../Low_coverage/YRI/YRI.header.mtdna_no",header=T,sep="\t")
admx.copyn<-read.table("../Low_coverage/Americans_lowcov.mtdna_no/Americans_lowcov.mtdna_no",header=F,sep="\t")
colnames(admx.copyn)<-colnames(ceu.copyn)


ceu.copyn$pop<-"CEU"
yri.copyn$pop<-"YRI"

admx.pop<-read.table("../../../Data_tables/pop_1kg.txt",header=T,sep="\t")
admx.copyn<-join(admx.copyn,admx.pop[,c("IID","pop")])

all.copyn<-rbind(ceu.copyn,yri.copyn,admx.copyn)
#remove samples with less than 250 mtdna_nopyn - these are blood samples


#parent.copyn<-rbind(ceu.copyn,yri.copyn)
#parent.copyn<-parent.copyn[which(parent.copyn$chr!="MT"),]

# #plot mtdna_copyn
# ggplot(parent.copyn,aes(chr,mtdna_no,fill=pop))+
#   geom_boxplot(outlier.shape=NA)+
#   geom_point(aes(color=pop),position=position_jitterdodge(),alpha=0.5)+
#   theme_bw()
# 
# #plot coverage
# ggplot(parent.copyn,aes(chr,coverage,fill=pop))+
#   geom_boxplot(outlier.shape=NA)+
#   geom_point(aes(color=pop),position=position_jitterdodge(),alpha=0.5)+
#   theme_bw()

#calculate median copy number across chromosomes (excluding MT) for each individual
dall.copyn<-ddply(all.copyn[which(all.copyn$chr!="MT"),],.(IID,pop),summarize,med.mtcopyn=median(mtdna_no))

#write to file
write.table(dall.copyn,"mtdna_copyn_allsamples_06012018.txt",sep="\t",col.names=T,row.names=F,quote=F)
