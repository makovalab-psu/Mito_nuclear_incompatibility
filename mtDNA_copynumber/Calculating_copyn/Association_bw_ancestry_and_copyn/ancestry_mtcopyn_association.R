library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)


#load mtdna copynumber
mt.copyn<-read.table("mtdna_copyn_allsamples_06012018.txt",header=T,sep="\t")

#load info about DNA source
DNAsrc<-read.table("../../../Data_tables/1kg_DNAsource.txt",header=T,sep="\t")
colnames(DNAsrc)<-c("IID","pop","ebv.cvrg","dna.source","has.blood")

#add source info
mt.copyn.src<-join(mt.copyn,DNAsrc[,c(1,3,4,5)],by="IID")

#plot distribution of copy number and color DNA source
plt<-ggplot(mt.copyn.src,aes(med.mtcopyn))+
  geom_density()+
  geom_point(aes(med.mtcopyn,0.0005,color=dna.source),position = position_jitter(width=0,height=0.0005),alpha=0.7)+
  theme_bw()+
  scale_color_manual(values=c("#e41a1c","#377eb8","grey"))+
  labs(x="MT-DNA copy number",y="Density",color="DNA source")

ggsave('lowcov_mtdna_copyn_distribution_06012018.pdf',plt,height=5,width=10)

#remove all samples with mtcopyn less than 250
mt.copyn<-mt.copyn[which(mt.copyn$med.mtcopyn>250),]

#read global ancestry for these people
qfile<-read.table("../../../ADMIXTURE_ancestry/Autosome/1kg_nam_flipped.3.Q",header=F)
#read fam file corresponding to the qfile
fam<-read.table("../../../ADMIXTURE_ancestry/Autosome/1kg_nam_flipped.fam",header=F)
fam<-fam[,c(1:2)]
colnames(fam)<-c("FID","IID")
global.anc<-cbind(fam,qfile)

#read pop file
pop.file<-read.table("../../../Data_tables/pop_1kg.txt",header=T)
#add pop info to global.anc file
global.anc<-join(global.anc,pop.file,by="IID")

#C2 is African, C1 is European, and C3 is Native American
colnames(global.anc)[c(3:5)]<-c("European","African","Native_American")

#merge global ancestry with mtcopyn info
mt.copyn.anc<-join(mt.copyn[,c(1,3)],global.anc,by="IID")

#load mthaplogroup info
#load haplogroup information 
mt.haplo<-read.table("../../../Data_tables/AMR_mt.haplo",header=T,sep="\t")
colnames(mt.haplo)<-c("IID","Sex","Population","major_haplogroup")
#group haplo into regional categories
mt.haplo$region_haplogroup<-NA
mt.haplo$region_haplogroup[which(mt.haplo$major_haplogroup%in%c("A","B","C","D"))]<-"Native_American"
mt.haplo$region_haplogroup[which(mt.haplo$major_haplogroup%in%c("L"))]<-"African"
mt.haplo$region_haplogroup[which(mt.haplo$major_haplogroup%in%c("K","J","U","H","I","T","V","W"))]<-"European"
#two individuals carry a South Asian haplogroup (M). Remove from analyses
mt.haplo.red<-mt.haplo[which(is.na(mt.haplo$region_haplogroup)=="FALSE"),]

#add mt.haplogroup info to global ancestry info
global.anc.mt<-merge(global.anc,mt.haplo.red[,c("IID","major_haplogroup","region_haplogroup")],by="IID")

#melt data.frame
m.global.anc.mt<-melt(global.anc.mt,id.vars=c("IID","FID","sex","pop","major_haplogroup","region_haplogroup"))
colnames(m.global.anc.mt)[c(7,8)]<-c("ancestry","fraction")

#add MT copy number info
merged.dat.full<-merge(mt.copyn,m.global.anc.mt,by="IID")
merged.dat.full$ancestry<-as.character(merged.dat.full$ancestry)

#estimate combined slope of effect of mtDNA and ancestry from a different mt haplogroup in admixed samples
gen.dissimilar<-merged.dat.full[which(merged.dat.full$region_haplogroup==merged.dat.full$ancestry),]
gen.dissimilar<-gen.dissimilar[-which(gen.dissimilar$pop.x%in%c("CEU","YRI")),]
gen.dissimilar$opp.fraction<-1-gen.dissimilar$fraction

gen.dissimilar$scaled.mtcopyn<-scale(gen.dissimilar$med.mtcopyn)
gen.dissimilar$scaled.opp.fraction<-scale(gen.dissimilar$opp.fraction)
summary(lm(data=gen.dissimilar,scaled.mtcopyn~scaled.opp.fraction))

#plot distribution of copy number in parental populations

gen.plt.parent<-ggplot(merged.dat.full[which(merged.dat.full$pop.x%in%c("CEU","YRI")),],aes(pop.x,med.mtcopyn))+
  geom_boxplot(aes(fill=pop.x),alpha=0.8)+
  geom_point(alpha=0.6,aes(x=pop.x,y=med.mtcopyn),position=position_jitter(width=0.3))+
  theme_bw()+
  geom_hline(yintercept=900,color="red",linetype="dashed",size=1)+
  labs(x="Parental groups",y="mtDNA copy number")+
  ylim(c(300,1750))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),legend.position="none")+
  scale_fill_manual(values = c("#d95f02","#7570b3"))

gen.dissimilar.plt.admx<-ggplot(gen.dissimilar,aes(opp.fraction,med.mtcopyn))+
  geom_point(alpha=0.6)+
  stat_smooth(method="lm",se=F)+
  theme_bw()+
  geom_hline(yintercept=900,color="red",linetype="dashed",size=1)+
  labs(x="Mito-nuclear Discordance",y="MTDNA copy number")+
  ylim(c(300,1750))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14),axis.text.y=element_blank(),axis.title.y=element_blank())

plt.combined<-grid.arrange(gen.plt.parent,gen.dissimilar.plt.admx,nrow=1)

ggsave("mitonuclear_copynumber_assoc_parents.pdf",plt.combined,height=7,width=10)

#plot ancestry vs mtcopy number broken down by ancestry and mtDNA haplogroup
plt.anc.v.mtcopy<-ggplot(merged.dat.full)+
  geom_point(aes(fraction,med.mtcopyn),alpha=0.6)+
  facet_grid(region_haplogroup~ancestry)+
  theme_bw()+
  geom_smooth(aes(fraction,med.mtcopyn),method="lm",se=F)+
  theme(axis.text.x=element_text(angle=90))+
  labs(x="Autosomal ancestry fraction",y="MTDNA copy number")
