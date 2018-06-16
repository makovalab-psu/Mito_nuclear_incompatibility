#script to generate bootstrapped loal ancestry estimates for genic windows
#load dependencies
library(data.table)
library(plyr)

args<-commandArgs(TRUE)
pop<-args[1]
n.bootstraps<-as.numeric(args[2])
weight<-as.logical(args[3])

message("loading lanc file")

#load gene-based local ancestry estimates
lanc.mt<-fread(paste("mt_lanc_annot_",pop,".bed",sep=""),header=F,sep="\t")

colnames(lanc.mt)<-c("gene.chr","gene.start","gene.stop","gene.name","broad.category","narrow.category","snp.chr","snp.start","snp.stop","eur.lanc","afr.lanc","nat.lanc")

#remove rows where snps are missing
lanc.mt<-lanc.mt[which(lanc.mt$snp.chr!="."),]
lanc.mt$eur.lanc<-as.numeric(lanc.mt$eur.lanc)
lanc.mt$afr.lanc<-as.numeric(lanc.mt$afr.lanc)
lanc.mt$nat.lanc<-as.numeric(lanc.mt$nat.lanc)

message("generating list of genes")
#generate list of gene names
genes<-ddply(lanc.mt,.(narrow.category),function(x){return(data.frame(gene.name=unique(x$gene.name)))})
genes$gene.name<-as.character(genes$gene.name)

message("calculating expected ancestry")
avg.eur.lanc<-mean(lanc.mt$eur.lanc)
avg.afr.lanc<-mean(lanc.mt$afr.lanc)
avg.nat.lanc<-mean(lanc.mt$nat.lanc)

#write function to bootstrap for a specific function category
boot.dev<-function(functional_cat,nboot,weighted=FALSE){
  lanc.deviation<-matrix(NA,nrow=nboot,ncol=3)
  
  if(weighted==TRUE){
  pb<-txtProgressBar(min=0,max=nboot,style=3)
  for(i in 1:nboot){
  sample.genes<-sample(genes$gene.name[which(genes$narrow.category==functional_cat)],164,replace=T)
  sample.lanc<-lanc.mt[which(lanc.mt$gene.name%in%sample.genes),]
  lanc.deviation[i,1]<-mean(sample.lanc$eur.lanc)-avg.eur.lanc
  lanc.deviation[i,2]<-mean(sample.lanc$afr.lanc)-avg.afr.lanc
  lanc.deviation[i,3]<-mean(sample.lanc$nat.lanc)-avg.nat.lanc
  setTxtProgressBar(pb,i)
  }
  colnames(lanc.deviation)<-c("European","African","Native American")
  lanc.deviation<-as.data.frame(lanc.deviation)
  lanc.deviation$category<-functional_cat
  return(lanc.deviation)
  }
  
  if(weighted==FALSE){
    pb<-txtProgressBar(min=0,max=nboot,style=3)
    for(i in 1:nboot){
      sample.genes<-sample(genes$gene.name[which(genes$narrow.category==functional_cat)],164,replace=T)
      sample.lanc<-lanc.mt[which(lanc.mt$gene.name%in%sample.genes),]
      dsample.lanc<-ddply(sample.lanc,.(gene.name),summarize,eur.lanc=mean(eur.lanc),afr.lanc=mean(afr.lanc),nat.lanc=mean(nat.lanc))
      lanc.deviation[i,1]<-mean(dsample.lanc$eur.lanc)-avg.eur.lanc
      lanc.deviation[i,2]<-mean(dsample.lanc$afr.lanc)-avg.afr.lanc
      lanc.deviation[i,3]<-mean(dsample.lanc$nat.lanc)-avg.nat.lanc
      setTxtProgressBar(pb,i)
    }
    colnames(lanc.deviation)<-c("European","African","Native American")
    lanc.deviation<-as.data.frame(lanc.deviation)
    lanc.deviation$category<-functional_cat
    return(lanc.deviation)
  }
}

message("bootstrapping high-mt genes")
high.lanc<-boot.dev("high",n.bootstraps, weight)

message("bootstrapping low-mt genes")
low.lanc<-boot.dev("low",n.bootstraps, weight)

message("bootstrapping non-mt genes")
non.lanc<-boot.dev("non_mito",n.bootstraps, weight)

comb.lanc<-rbind(high.lanc,low.lanc,non.lanc)

message("writing to file")

if(weight==TRUE){
  fwrite(comb.lanc,paste(pop,"_lanc_deviation_bootstrap_weighted.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}

if(weight==FALSE){
  fwrite(comb.lanc,paste(pop,"_lanc_deviation_bootstrap_unweighted.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}


