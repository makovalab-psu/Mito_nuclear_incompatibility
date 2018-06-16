library(data.table)
library(ggplot2)

#write function to simulate variation in local ancestry
sim.admix<-function(a0,b0,c0,g,nind=1250){
  nchrom<-nind*2
  nchrom.a0<-a0*nchrom
  nchrom.b0<-b0*nchrom
  nchrom.c0<-c0*nchrom
  #nchrom.present=ne*2
  #nchrom0<-nchrom.present*(a0+b0+c0)
  #nchrom.a0<-a0*nchrom.present
  #nchrom.b0<-b0*nchrom.present
  #nchrom.c0<-c0*nchrom.present
  for(i in 1:g){
    fa0<-nchrom.a0/(nchrom.a0+nchrom.b0+nchrom.c0)
    fb0<-nchrom.b0/(nchrom.a0+nchrom.b0+nchrom.c0)
    fc0<-nchrom.c0/(nchrom.a0+nchrom.b0+nchrom.c0)
    
    sims<-rmultinom(nchrom,1,c(fa0,fb0,fc0))
    updated.p<-apply(sims,1,function(x){length(which(x==1))/length(x)})
    nchrom.a0<-updated.p[1]*nchrom
    nchrom.b0<-updated.p[2]*nchrom
    nchrom.c0<-updated.p[3]*nchrom
  }
  
  return(c(nchrom.a0,nchrom.b0,nchrom.c0)/nchrom)
}



#read local ancestry frequency for the population of interest
avg.lanc<-fread('lanc_popavg_CLM.txt',sep="\t",header=T)
avg.lanc$hg19chr<-paste("chr",avg.lanc$chr,sep="")

#subsample 1000 loci
avg.lanc.red<-avg.lanc[sample(nrow(avg.lanc),10000),]

output<-as.data.frame(t(replicate(10000,sim.admix(mean(avg.lanc.red$afr.lanc),mean(avg.lanc.red$eur.lanc),mean(avg.lanc.red$nat.lanc),14))))
colnames(output)<-c("African","European","Native_American")


#generate quantiles for qqplot
qqdat<-function(nq,variable1,variable2){
  #nq is number of quantiles
  p <- (1 : nq) / nq - 0.5 / nq
  qvar1<-quantile(variable1,p)
  qvar2<-quantile(variable2,p)
  return(data.frame(qvar1,qvar2))
}

qq.afr<-qqdat(100,output$African,avg.lanc.red$afr.lanc)
qq.afr$anc.pop<-"African"
qq.nat<-qqdat(100,output$Native_American,avg.lanc.red$nat.lanc)
qq.nat$anc.pop<-"Native_American"
qq.eur<-qqdat(100,output$European,avg.lanc.red$eur.lanc)
qq.eur$anc.pop<-"European"
qq.all<-rbind(qq.afr,qq.nat,qq.eur)

#plot simulated vs observed values of local ancestry
qq.plot<-ggplot(qq.all,aes(qvar1,qvar2))+
  geom_point(alpha=0.7)+
  geom_abline(intercept=0,slope=1,color="red")+
  labs(x="Simulated local ancestry",y="Observed local ancestry")+
  facet_wrap(~anc.pop,scales="free")+
  theme_bw()+
  theme(axis.text=element_text(size=14),axis.title = element_text(size=18),panel.spacing.x = unit(1,"cm"))

ggsave("qq_autosomal_CLM.pdf",qq.plot,height=5,width=10)

#sim.mt


#loading sex specific ancestry contributions for each population
mf.props<-read.table("../../Sex_bias_estimation//sex_spec_ancestry_contribs_avg.txt",header=T)
mf.clm<-mf.props[which(mf.props$admx.pop=="CLM" & mf.props$sex=="f"),]

#write function to simulate variation in mtDNA given proportion of females
sim.mt<-function(af,bf,cf,g,nind=1250){
  #mt has 1/4th the effective population size of the autosome
  #or half the number of individuals since only females pass it on
  nchrom<-nind/2
  #calculate number of mtDNA of each ancestry from the frequency of females
  nchrom.af<-af*nchrom
  nchrom.bf<-bf*nchrom
  nchrom.cf<-cf*nchrom
  #run drift simulation for g generations
  for(i in 1:g){
    af<-nchrom.af/(nchrom.af+nchrom.bf+nchrom.cf)
    bf<-nchrom.bf/(nchrom.af+nchrom.bf+nchrom.cf)
    cf<-nchrom.cf/(nchrom.af+nchrom.bf+nchrom.cf)
    
    #multinomial sampling of the three haplogroups
    sims<-rmultinom(nchrom,1,c(af,bf,cf))
    updated.p<-apply(sims,1,function(x){length(which(x==1))/length(x)})
    nchrom.af<-updated.p[1]*nchrom
    nchrom.bf<-updated.p[2]*nchrom
    nchrom.cf<-updated.p[3]*nchrom
  }
  
  return(c(nchrom.af,nchrom.bf,nchrom.cf)/nchrom)
}

#simulate
mt.sim.dat<-as.data.frame(t(replicate(10000,sim.mt(mf.clm$proportion[3],mf.clm$proportion[2],mf.clm$proportion[3],14))))
colnames(mt.sim.dat)<-c("African","European","Native_American")

#read observed mt frequencies
mt.obs.dat<-read.table("../../Data_tables//mtprop_table.txt",header=T)
mt.obs.pur<-mt.obs.dat[which(mt.obs.dat$Population=="CLM"),]

#plt simulated data against observed mt haplogroup freuency
mt.sim.dat<-melt(mt.sim.dat)
colnames(mt.sim.dat)<-c("Ancestry","sim.f")

mt.plt<-ggplot()+
  geom_violin(data=mt.sim.dat,aes(Ancestry,sim.f,fill=Ancestry))+
  geom_boxplot(data=mt.sim.dat,aes(Ancestry,sim.f,fill=Ancestry),width=0.5,outlier.shape=NA)+
  geom_point(data=mt.obs.pur,aes(Ancestry,mt),color="red")+
  theme_bw()+
  scale_fill_manual(values=c("#8da0cb","#fc8d62","#66c2a5"))+
  theme(legend.position="none")+
  labs(x="Ancestry",y="Simulated/Observed MT haplogroup frequency")

ggsave("Simulated_v_observed_MT_CLM.pdf",mt.plt,height=7,width=7)
