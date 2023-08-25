#input_files
#multiple files
#save R

args = commandArgs(trailingOnly =TRUE)
n_input_files<-(length(args)-2)/2
ploidy_list_file<-args[1]
pop_map_file<-args[2]
output_image<-args[3]

structure_out_clumpp_all<-args[4:(n_input_files+3)]
structure_plot_clumpp_all<-args[(n_input_files+4):((2*n_input_files)+3)]

ploidy_list<-read.table(ploidy_list_file, sep="\t", header=FALSE)
pop_map<-read.table(pop_map_file, sep="\t", header=FALSE)
pop_table<-table(pop_map$V2)

#Plot admixture plots for all

plot_bars<-function(value_matrix, pop_table, ploidy_list, output_pdf)
{
  pdf(output_pdf, width=12, height=8)
  par(mar=c(5.1, 3.1, 1.5, 5.1), xpd=TRUE)
  #barplot(t(value_matrix), border = NA, space = 0,
  barplot(value_matrix, border = NA, space = 0,
          xlab = "", ylab = "Admixture coefficients", col=rainbow(nrow(value_matrix)))
  abpos=0
  for(i in 1:length(pop_table)){
    abpos<-abpos+pop_table[i]
    abline(v=abpos)
    text(abpos-(pop_table[i]/2), -0.075, names(pop_table)[i], srt=90, xpd=NA)
  }
  points(0.5:length(ploidy_list), rep(1.01, length(ploidy_list)), cex=0.5, pch=ploidy_list, xpd=NA)
  legend("topright",legend=unique(ploidy_list), pch=unique(ploidy_list), inset=c(-0.05,0), title="Ploidy")
  dev.off()
}

for(structure_out_clumpp_all_i in 1:length(structure_out_clumpp_all)){
  structure_out_clumpp_all_temp<-read.table(structure_out_clumpp_all[structure_out_clumpp_all_i], header=FALSE)
  structure_q_k_split<-tail(unlist(strsplit(structure_out_clumpp_all[structure_out_clumpp_all_i], "_")), n=1)
  structure_q_k<-head(unlist(strsplit(structure_q_k_split, ".", fixed=TRUE)), n=1)
  structure_q_k_matrix<-as.matrix(structure_out_clumpp_all_temp[,6:(6+as.numeric(structure_q_k)-1)])
  #structure_q_k_matrix<-matrix(structure_out_clumpp_all_temp[,6:(6+as.numeric(structure_q_k)-1)], ncol=as.numeric(structure_q_k))
  print(dim(structure_q_k_matrix))
  plot_bars(t(structure_q_k_matrix), pop_table, ploidy_list$V1, structure_plot_clumpp_all[structure_out_clumpp_all_i])
}

#Plot Rhat


save.image(output_image)
