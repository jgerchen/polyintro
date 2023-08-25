library(MASS)
args = commandArgs(trailingOnly=TRUE)
PC_file<-args[1]
max_k<-args[2]
prefix_qk<-args[3]
input_PCS<-read.table(PC_file)
nind<-dim(input_PCS)[1]
for(k in 2:max_k){
  init.admix<-matrix(0, nrow=nind, ncol=k)
  kn<-kmeans(input_PCS[,1:5], k, iter.max=10, nstart=10, algorithm="Hartigan-Wong")
  ldakn<-lda(input_PCS[,1:5], grouping=kn$cluster, CV=T)
  init.admix<-ldakn$posterior
  write.table(round(init.admix, 2), quote=F, row.names=F, col.names=F, file=paste0(prefix_qk,"_qk",k))
}
