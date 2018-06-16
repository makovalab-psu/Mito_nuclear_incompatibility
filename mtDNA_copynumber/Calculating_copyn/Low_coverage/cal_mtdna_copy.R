#!/usr/bin/env Rscript

args<-commandArgs(TRUE)

if (length(args)<2) {
  stop("At least two arguments must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==2) {


idxstats<-read.table(args[1],header=F)
read_length<-as.numeric(args[2])
colnames(idxstats)<-c("chr","bp","mapped","unmapped")
idxstats<-idxstats[which(idxstats$chr%in%c(1:22,"MT")),]

idxstats$coverage<-idxstats$mapped*read_length/idxstats$bp

#coverage<-read.table(args[1],header=F)

#colnames(coverage)<-c("chr","bp","mapped","med_coverage")

name<-unlist(strsplit(as.character(args[1]),split=".",fixed=T))[1]
print(name)
#coverage$name<-as.character(name)


mtcov<-idxstats$coverage[which(idxstats$chr=="MT")]
idxstats$mtdna_no<-with(idxstats,(mtcov/coverage)*2)
idxstats$IID<-name

write.table(idxstats,paste(name,".mtdna_no",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}
