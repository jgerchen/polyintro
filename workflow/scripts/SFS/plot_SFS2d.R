options(scipen=999)
#input_files
#multiple files
#save R
library(viridis)
args = commandArgs(trailingOnly =TRUE)
input_file_name<-args[1]
output_file_name<-args[2]
################Delete me!!!!
#input_file_name<-"~/Desktop/Stuff/Medicago/SFS/medicago_SFS_2d.txt"
#output_file_name<-"~/Desktop/Stuff/Medicago/SFS/medicago_SFS_2d_plot"
#############################
input_file <- file(input_file_name ,open="r")
input_lines<-readLines(input_file)

curr_matrix<-matrix()
matrices<-list()
n_comparisons<-0
for(input_line in input_lines){
  if(substring(input_line, 1,1)=="#"){
    curr_comparison<-substring(input_line, 2,nchar(input_line))
    #add old matrix to matrix list
    if(n_comparisons!=0){
      matrices[[n_comparisons]]<-curr_matrix
      names(matrices)[n_comparisons]<-curr_comparison
    }
    n_comparisons<-n_comparisons+1
    new_matrix<-TRUE
  }
  else{
    items=unlist(strsplit(input_line, "\t"))
    values=matrix(as.numeric(items[1:length(items)]), nrow = 1)
    if(new_matrix==TRUE){
      curr_matrix<-matrix(values, nrow = 1)      
      new_matrix<-FALSE
    }
    else{
      curr_matrix<-rbind(curr_matrix, values)
    }
  }
}
matrices[[n_comparisons]]<-curr_matrix
names(matrices)[n_comparisons]<-curr_comparison





#heatmap(log(matrices$MF37_MF47), Rowv=NA, Colv=NA, revC=TRUE, labRow = 0:nrow(matrices$MF37_MF47), labCol =  0:ncol(matrices$MF37_MF47), col=viridis(256), scale="none")
#for(out_matrix in names(matrices)){
#  print(out_matrix)
#  assign(out_matrix, out_matrix)
#  matrices$out_matrix
#}


for(out_matrix_i in 1:length(matrices)){
  comp_name<-names(matrices)[out_matrix_i]
  comp_matrix<-unname(matrices)[[out_matrix_i]]
  pdf(paste(output_file_name,"_", comp_name,".pdf", sep=""), width=10, height=8)
  #par(mfrow = c(length(input_lines), 2))
  layout_mat <- matrix(c(1, 1, 2, 0), nrow = 2)
  my_lay <- layout(mat = layout_mat, 
                   heights = c(3, 2),
                   widths = c(4, 1), respect =TRUE)
  #layout.show(my_lay)
  par(mar=c(5.1, 4.1, 2.5, 1.1), xpd=TRUE)
  comp_1<-strsplit(comp_name, "_")[[1]][1]
  comp_2<-strsplit(comp_name, "_")[[1]][2]
  q_arr<-c(quantile(comp_matrix, 0:254/255), max(comp_matrix))
  image(comp_matrix, col=viridis(255), breaks=q_arr, xaxt="n", yaxt="n", xlab=comp_1, ylab=comp_2)
  axis(1, at=0:(dim(comp_matrix)[1]-1)/(dim(comp_matrix)[1]-1), labels=0:(dim(comp_matrix)[1]-1))
  axis(2, at=0:(dim(comp_matrix)[2]-1)/(dim(comp_matrix)[2]-1), labels=0:(dim(comp_matrix)[2]-1))
  image(matrix(q_arr, nrow=1), col=viridis(255), breaks=q_arr, xaxt="n", yaxt="n")
  axis(2, at=c(0, (0:255/255)[seq(0, 255, 32)],1), labels=c(0, round(q_arr[seq(0, 255, 32)]), max(comp_matrix)), cex.axis=0.7)
  dev.off()
}




#par(mfrow = c(1, 2))











