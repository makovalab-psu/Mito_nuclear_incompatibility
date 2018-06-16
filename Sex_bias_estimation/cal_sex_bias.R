library(plyr)
library(reshape2)
library(ggplot2)

#setwd("/Users/Azaidi/Documents/mtproj_files/mitonuclear_project/mtnuc_organization/sex_biased_admixture")

#read pop file
pop<-read.table('../Data_tables/pop_1kg.txt',header=T,stringsAsFactors = F)
colnames(pop)<-c("IID","Sex","Population")

#define admixed populations
admx.pops<-c("ACB","ASW","CLM","MXL","PEL","PUR")

#read autosomal ancestry
#glanc.autosome<-read.table("~/Documents/mtproj_files/mitonuclear_project/local_ancestry_1kgcalls/viterbi/clean_viterbi/AMR_04182018.autosome.glanc",header=T)
qfile.autosome<-read.table("~/Documents/mtproj_files/mitonuclear_project/mtnuc_organization/admixture_results/Autosome/1kg_nam_flipped.3.Q",header=F)
colnames(qfile.autosome)<-c("European","African","Native_American")
auto.fam<-read.table('~/Documents/mtproj_files/mitonuclear_project/mtnuc_organization/admixture_results/Autosome/1kg_nam_flipped.fam',header=F)
auto.fam<-auto.fam[,c(1:2)]
colnames(auto.fam)<-c("FID","IID")
qfile.autosome<-cbind(auto.fam,qfile.autosome)
qfile.autosome<-merge(qfile.autosome,pop,by="IID")
#glanc.autosome<-merge(glanc.autosome,pop,by="IID")
#colnames(glanc.autosome)[c(2:4)]<-c("European","African","Native_American")

#calculate mean ancestry across all individuals
#fpop.auto<-ddply(glanc.autosome[,c('Population','European','African','Native_American')],.(Population),colwise(mean))

#calculate average autosomal ancestry for each population
f.pop.auto<-ddply(qfile.autosome[,c("European","African","Native_American","Population")],.(Population),summarize,European=mean(European),African=mean(African),Native_American=mean(Native_American))



#read X chromosome ancestry
qfile.x<-read.table("../ADMIXTURE_ancestry/X_chromosome/Americans_1kg_namsnps_X_pruned.3.Q",header=F)
colnames(qfile.x)<-c("African","Native_American","European")
x.fam<-read.table('../ADMIXTURE_ancestry/X_chromosome/Americans_1kg_namsnps_X_pruned.fam',header=F)
x.fam<-x.fam[,c(1:2)]
colnames(x.fam)<-c("FID","IID")
qfile.x<-cbind(x.fam,qfile.x)
qfile.x<-merge(qfile.x,pop,by="IID")

#define function to calculate average X chromosomal ancestry for a population
#this function will adjust for ploidy differences between males and females
cal.avg.x<-function(x){
  males<-x[which(x$Sex=="male"),]
  females<-x[which(x$Sex=="female"),]
  males.afr<-mean(males$African)
  males.eur<-mean(males$European)
  males.nat<-mean(males$Native_American)
  
  females.afr<-mean(females$African)
  females.eur<-mean(females$European)
  females.nat<-mean(females$Native_American)
  
  African<-(males.afr+2*females.afr)/3
  European<-(males.eur+2*females.eur)/3
  Native_American<-(males.nat+2*females.nat)/3
  return(data.frame(European,African,Native_American))
}

f.pop.x<-ddply(qfile.x,.(Population),cal.avg.x)

#create grid of values for each population - 3 decimals

f.male<-seq(0,0.5,by=0.001)
f.female<-seq(0,0.5,by=0.001)
mat1<-expand.grid(f=f.female,m=f.male)

f.pop<-function(anc.pop,admx.pop){
  auto.fractions<-f.pop.auto[which(f.pop.auto$Population==admx.pop),]
  x.fractions<-f.pop.x[which(f.pop.x$Population==admx.pop),]
  #anc.pop is the name of the ancestral population e (European, African, Native_American)
  #pop is the name of the current population e (PUR, CLM, ASW, ACB, MXL, PEL)
  auto.p<-as.numeric(round(auto.fractions[,anc.pop],3))
  x.p<-as.numeric(round(x.fractions[,anc.pop],3))
  #if(anc.pop=="African"){mat=afr.mat}
  #if(anc.pop=="European"){mat=eur.mat}
  #if(anc.pop=="Native_American"){mat=nat.mat}
  
  #select values from grid that match autosomal fraction
  comb.auto.fraction<-round(mat1$f+mat1$m,3)
  mat<-mat1[which(comb.auto.fraction==auto.p),]
  mat$pred.x<-(mat$m+(2*mat$f))/1.5
  mat$dev<-(mat$pred.x-x.p)^2
  row.index<-which(mat$dev==min(mat$dev))
  if(length(row.index)>1){row.index<-row.index[1]}
  mat<-cbind(anc.pop,mat[row.index,c('f','m')])
  return(mat)
}

sex.bias.pop<-function(pop){
  #pop is the name of the current population e (PUR, CLM, ASW, ACB, MXL, PEL)
  nat<-f.pop("Native_American",pop)
  eur<-f.pop("European",pop)
  afr<-f.pop("African",pop)
  dat<-rbind(nat,eur,afr)
  dat$admx.pop<-pop
  return(dat)
}

mfestimates<-sapply(admx.pops,sex.bias.pop,simplify=F,USE.NAMES=T)
mfestimates<-do.call(rbind,mfestimates)
mfestimates<-melt(mfestimates,id.vars=c("anc.pop","admx.pop"))
colnames(mfestimates)[c(3,4)]<-c("sex","proportion")

#normalize such that m = f = 0.5
mfestimates<-ddply(mfestimates,.(admx.pop,sex),function(x){
  norm.prop<-x$proportion/sum(x$proportion)*0.5
  x$norm.proportion<-norm.prop
  return(x)
}
)

#plot m and f proportions for each population
mf.plt<-ggplot(mfestimates,aes(sex,norm.proportion,fill=anc.pop))+
  geom_bar(stat="density")+
  facet_wrap(~admx.pop)+
  scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))+
  theme_bw()+
  labs(x="Sex",y="Average contribution",fill="Ancestry")

ggsave('Sex_specific_ancestry_contributions_avg.pdf',mf.plt,height=7,width=10)
write.table(mfestimates,"sex_spec_ancestry_contribs_avg.txt",row.names=F,col.names=T,quote=F,sep="\t")


#generate bootstrap intervals for m and f contributions
mfest.boot<-mfestimates
mfest.boot$bootstrap<-0
pb<-txtProgressBar(min=,max=1000,style=3)
for(i in 1:1000){
  f.pop.auto<-ddply(qfile.autosome[,c("European","African","Native_American","Population")],.(Population),function(x){
    boot.index<-sample(nrow(x),replace=T)
    x2<-x[boot.index,]
    European=mean(x2$European)
    African=mean(x2$African)
    Native_American=mean(x2$Native_American)
    return(data.frame(European,African,Native_American))})
  f.pop.x<-ddply(qfile.x,.(Population),function(x){
    boot.index<-sample(nrow(x),replace=T)
    x2<-x[boot.index,]
    return(cal.avg.x(x2))})
  
  mfestimates<-sapply(admx.pops,sex.bias.pop,simplify=F,USE.NAMES=T)
  mfestimates<-do.call(rbind,mfestimates)
  mfestimates<-melt(mfestimates,id.vars=c("anc.pop","admx.pop"))
  colnames(mfestimates)[c(3,4)]<-c("sex","proportion")
  #normalize such that m = f = 0.5
  mfestimates<-ddply(mfestimates,.(admx.pop,sex),function(x){
    norm.prop<-x$proportion/sum(x$proportion)*0.5
    x$norm.proportion<-norm.prop
    return(x)
  }
  )
  mfestimates$bootstrap<-i
  mfest.boot<-rbind(mfest.boot,mfestimates)
  setTxtProgressBar(pb,i)
  
}


mfest.boot<-mfest.boot[-which(mfest.boot$bootstrap==0),]
dmfest.boot<-ddply(mfest.boot,.(anc.pop,admx.pop,sex),summarize,mean.prop=mean(norm.proportion),lower=quantile(norm.proportion,probs=0.025),upper=quantile(norm.proportion,probs=0.975))
boot.plt<-ggplot(dmfest.boot)+
  geom_bar(aes(sex,mean.prop,fill=anc.pop,color=anc.pop),alpha=0.7,stat="identity",position=position_dodge())+
  facet_wrap(~admx.pop)+
  theme_bw()+
  geom_errorbar(aes(sex,ymin=lower,ymax=upper,color=anc.pop),position=position_dodge(width=0.9),width=0.2)+
  scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))+
  scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))+
  labs(x="Sex",y="Proportion of females/males from ancestral population")

ggsave('Sex_specific_ancestry_contributions_bootstrap.pdf',boot.plt,height=7,width=7)
write.table(mfest.boot,"sex_spec_ancestry_contribs_bootstap.txt",row.names=F,col.names=T,quote=F,sep="\t")
