#binput_files
#multiple files
#save R
library(viridis)
args = commandArgs(trailingOnly =TRUE)
#ploidy_list_file<-args[1]
input_file=args[1]
pop_map_file<-args[2]
output_image<-args[3]

pops<-read.table(pop_map_file, sep="\t", header=FALSE, col.names=c("ind", "pop", "ploidies" ))
results<-read.table(input_file, header=FALSE)
k=ncol(results)-5
pops$mean_ploidy<-pops$pop
populations<-unique(pops$pop)
for(population in populations){
	mean_ploidy<-mean(pops$ploidies[pops$pop==population])
	pops$mean_ploidy[pops$pop==population]<-mean_ploidy
}

all_data<-cbind(pops, results[, 6:ncol(results)])
all_data_sorted<-all_data[order(all_data$mean_ploidy, all_data$pop, all_data$ploidies),]

pdf(output_image, width=nrow(pops)*0.2, height=8)
par(mar=c(5.1, 4.1, 1.5, 1.1), xpd=TRUE)
barplot(t(as.matrix(all_data_sorted[,5:ncol(all_data_sorted)])), border=NA, space=0, col=viridis(k), xlab="", ylab="admixture coefficients", names=rep('',length(all_data_sorted$ind)))
#  abpos=0
curr_pop<-"NNN"
curr_ploidy<--1
last_pop_boundary<-0
last_ploidy_boundary<-0
for(i in 1:nrow(all_data_sorted)){
  print(i)
	if(all_data_sorted$pop[i]!=curr_pop){
		if(curr_pop!="NNN"){
			abline(v=i-1)
			text(last_pop_boundary+((i-last_pop_boundary)/2), 1.03, labels=curr_pop, adj=0.5, xpd=NA)
			last_pop_boundary=i-1
			last_ploidy_boundary=i-1
		}
		curr_pop<-all_data_sorted$pop[i]
		curr_ploidy<-all_data_sorted$ploidies[i]
	}
	#only plot ploidy when we're not switching pops
	else{
		if(all_data_sorted$ploidies[i]!=curr_ploidy){
			if(curr_ploidy!=-1){
				abline(v=i-1, col="blue", lty=3)
				#text(last_ploidy_boundary+(i/2), 1.03, curr_ploidy, "x", sep="", xpd=NA)
			}
			curr_ploidy=all_data_sorted$ploidies[i]
			last_ploidy_boundary=i-1
		}
	}
}
#print last pop
text(last_pop_boundary+((i-last_pop_boundary)/2), 1.03, labels=curr_pop, adj=0.5, xpd=NA)
dev.off()











#for(structure_out_clumpp_all_i in 1:length(structure_out_clumpp_all)){
#  structure_out_clumpp_all_temp<-read.table(structure_out_clumpp_all[structure_out_clumpp_all_i], header=FALSE)
#  a_pop<-cbind(pops,structure_out_clumpp_all_temp[,4] ,structure_out_clumpp_all_temp[,6:ncol(structure_out_clumpp_all_temp)])
#  a_pop<-a_pop[order(a_pop$pop, a_pop$ploidies),]
#  pop_table=table(a_pop$pop)
#  ploidy_table=table(a_pop$pop, a_pop$ploidies)
#  pdf(structure_plot_clumpp_all[structure_out_clumpp_all_i] ,width=nrow(pops)*0.2, height=8)
#  par(mar=c(5.1, 4.1, 1.5, 1.1), xpd=TRUE)
#  barplot(t(as.matrix(a_pop[,5:ncol(a_pop)])), border=NA, space=0, col=viridis(ncol(structure_out_clumpp_all_temp)-5), xlab="", ylab="admixture coefficients", names=rep('',length(a_pop$ind)))
#  abpos=0
#  for(i in 1:length(pop_table)){
#    ploi_pos=0
#    for(u in 1:length(ploidy_table[i,])){
#      if( ploidy_table[i,u]>0){
#        text(abpos+ploi_pos+(ploidy_table[i,u]/2), 1.03, paste(colnames(ploidy_table)[u], "x", sep=""), xpd=NA)
#        ploi_pos<-ploi_pos+ploidy_table[i,u]
#        if(ploi_pos<pop_table[i]){
#          abline(v=abpos+ploi_pos, col="blue", lty=3)
#        }
        
#      }
#    }
#    abpos<-abpos+pop_table[i]
#    abline(v=abpos)
#    text(abpos-(pop_table[i]/2), -0.075, names(pop_table)[i], srt=90, xpd=NA)
#  }
#  dev.off()
#}

#save.image(output_image)
