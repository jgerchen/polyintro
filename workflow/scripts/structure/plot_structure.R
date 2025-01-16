#input_files
#multiple files
#save R
library(viridis)
args = commandArgs(trailingOnly =TRUE)
n_input_files<-(length(args)-2)/2
ploidy_list_file<-args[1]
pop_map_file<-args[2]
output_image<-args[3]

structure_out_clumpp_all<-args[4:(n_input_files+3)]
structure_plot_clumpp_all<-args[(n_input_files+4):((2*n_input_files)+3)]

ploidies<-read.table(ploidy_list_file, sep="\t", header=FALSE)
pops<-read.table(pop_map_file, sep="\t", header=FALSE, col.names=c("ind", "pop" ))
pops$ploidies<-rep(0, length(pops$ind))
for(i in 1:length(pops$ind)){
  pops$ploidies[i]<-ploidies$V2[ploidies$V1==pops$ind[i]]
}


for(structure_out_clumpp_all_i in 1:length(structure_out_clumpp_all)){
  structure_out_clumpp_all_temp<-read.table(structure_out_clumpp_all[structure_out_clumpp_all_i], header=FALSE)
  a_pop<-cbind(pops,structure_out_clumpp_all_temp[,4] ,structure_out_clumpp_all_temp[,6:ncol(structure_out_clumpp_all_temp)])
  a_pop<-a_pop[order(a_pop$pop, a_pop$ploidies),]
  pop_table=table(a_pop$pop)
  ploidy_table=table(a_pop$pop, a_pop$ploidies)
  pdf(structure_plot_clumpp_all[structure_out_clumpp_all_i] ,width=nrow(pops)*0.2, height=8)
  par(mar=c(5.1, 4.1, 1.5, 1.1), xpd=TRUE)
  barplot(t(as.matrix(a_pop[,5:ncol(a_pop)])), border=NA, space=0, col=viridis(ncol(structure_out_clumpp_all_temp)-5), xlab="", ylab="admixture coefficients", names=rep('',length(a_pop$ind)))
  abpos=0
  for(i in 1:length(pop_table)){
    ploi_pos=0
    for(u in 1:length(ploidy_table[i,])){
      if( ploidy_table[i,u]>0){
        text(abpos+ploi_pos+(ploidy_table[i,u]/2), 1.03, paste(colnames(ploidy_table)[u], "x", sep=""), xpd=NA)
        ploi_pos<-ploi_pos+ploidy_table[i,u]
        if(ploi_pos<pop_table[i]){
          abline(v=abpos+ploi_pos, col="blue", lty=3)
        }
        
      }
    }
    abpos<-abpos+pop_table[i]
    abline(v=abpos)
    text(abpos-(pop_table[i]/2), -0.075, names(pop_table)[i], srt=90, xpd=NA)
  }
  dev.off()
}

save.image(output_image)
