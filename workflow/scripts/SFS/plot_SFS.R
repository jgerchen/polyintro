#input_files
#multiple files
#save R
args = commandArgs(trailingOnly =TRUE)
input_file_name<-args[1]
output_file_name<-args[2]


input_file <- file(input_file_name ,open="r")
input_lines<-readLines(input_file)





pdf(output_file_name, width=10, height=8*length(input_lines))
par(mfrow = c(length(input_lines), 1))
par(mar=c(5.1, 4.1, 1.5, 1.1), xpd=TRUE)
for(input_line in input_lines){
  items=unlist(strsplit(input_line, "\t"))
  pop=items[1]
  values=as.numeric(items[2:length(items)])
  barplot(values, names.arg = 0:(length(values)-1), space = 0)
  title(pop)
}

dev.off()


