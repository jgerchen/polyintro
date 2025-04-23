library(StAMPP)
library(vcfR)
library(adegenet)
library(phangorn)
library(viridis)
args = commandArgs(trailingOnly=TRUE)
vcf_input<-args[1]
#ploidy_map<-args[2]
pop_map<-args[2]
plot_heatmap_pop<-args[3]
plot_heatmap_ind<-args[4]
plot_nj_pop<-args[5]
plot_nj_ind<-args[6]
plot_neighnet_pop<-args[7]
plot_neighnet_ind<-args[8]
out_nei_ind_phylip<-args[9]

#ploidies<-read.table(ploidy_map, sep="\t", header=FALSE)
pops<-read.table(pop_map, sep="\t", header=FALSE, col.names=c("ind", "pop", "ploidies"))
#pops$ploidies<-rep(0, length(pops$ind))
#for(i in 1:length(pops$ind)){
#  pops$ploidies[i]<-ploidies$V2[ploidies$V1==pops$ind[i]]
#}


dip_colour<-"#1e90ff"
tet_colour<-"#ffa500"
tri_colour<-"#db3192"
#no 5x colour yet
pent_colour<-"#bcdb31"
hex_colour<-"#7cfc00"

pops$colours<-pops$ploidies
pops$colours[pops$colours==2]<-dip_colour
pops$colours[pops$colours==3]<-tri_colour
pops$colours[pops$colours==4]<-tet_colour
pops$colours[pops$colours==5]<-pent_colour
pops$colours[pops$colours==6]<-hex_colour


vcfR2genlight.tetra <- function (x, n.cores = 1) 
{
  bi <- is.biallelic(x)
  if (sum(!bi) > 0) {
    msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
    msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
    msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
    warning(msg)
    x <- x[bi, ]
  }
  x <- addID(x)
  CHROM <- x@fix[, "CHROM"]
  POS <- x@fix[, "POS"]
  ID <- x@fix[, "ID"]
  x <- extract.gt(x)
  x[x == "0|0"] <- 0
  x[x == "0|1"] <- 1
  x[x == "1|0"] <- 1
  x[x == "1|1"] <- 2
  x[x == "0/0"] <- 0
  x[x == "0/1"] <- 1
  x[x == "1/0"] <- 1
  x[x == "1/1"] <- 2
  x[x == "1/1/1"] <- 3
  x[x == "0/1/1"] <- 2
  x[x == "0/0/1"] <- 1
  x[x == "0/0/0"] <- 0
  x[x == "1/1/1/1"] <- 4
  x[x == "0/1/1/1"] <- 3
  x[x == "0/0/1/1"] <- 2
  x[x == "0/0/0/1"] <- 1
  x[x == "0/0/0/0"] <- 0
  x[x == "1/1/1/1/1"] <- 5
  x[x == "0/1/1/1/1"] <- 4
  x[x == "0/0/1/1/1"] <- 3
  x[x == "0/0/0/1/1"] <- 2
  x[x == "0/0/0/0/1"] <- 1
  x[x == "0/0/0/0/0"] <- 0
  x[x == "1/1/1/1/1/1"] <- 6
  x[x == "0/1/1/1/1/1"] <- 5
  x[x == "0/0/1/1/1/1"] <- 4
  x[x == "0/0/0/1/1/1"] <- 3
  x[x == "0/0/0/0/1/1"] <- 2
  x[x == "0/0/0/0/0/1"] <- 1
  x[x == "0/0/0/0/0/0"] <- 0
  if (requireNamespace("adegenet")) {
    x <- new("genlight", t(x), n.cores = n.cores)
  }
  else {
    warning("adegenet not installed")
  }
  adegenet::chromosome(x) <- CHROM
  adegenet::position(x) <- POS
  adegenet::locNames(x) <- ID
  return(x)
}

vcf <- read.vcfR(vcf_input)
aa.genlight <- vcfR2genlight.tetra(vcf)
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_") # add real SNP.names

indnames<-indNames(aa.genlight)
#popnames<-rep(NA, length(indnames))
#for(i in 1:length(indnames)){
#  popnames[i]<-pops$pop[pops$ind==indnames[i]]
#}

pop_order<-match(indnames,pops$ind)
pop(aa.genlight)<-pops$pop[pop_order]
#pop(aa.genlight)<-popnames  # add pop names
include_individuals<-indNames(aa.genlight)[is.na(pop(aa.genlight))==F]
sub.genlight<-new("genlight",(as.matrix(aa.genlight)[include_individuals,]))

new_pop_order<-match(indNames(sub.genlight),pops$ind)
pop(sub.genlight)<-pops$pop[new_pop_order]

aa.D.pop <- stamppNeisD(sub.genlight, pop = TRUE)   # Nei's 1972 genetic distance between pops
aa.D.ind <- stamppNeisD(sub.genlight, pop = FALSE)   # Nei's 1972 genetic distance between pops
colnames(aa.D.pop)<-rownames(aa.D.pop)
colnames(aa.D.ind)<-rownames(aa.D.ind)

#fsts<-stamppFst(sub.genlight)

pdf(plot_heatmap_pop)
heatmap.bp(aa.D.pop)
dev.off()

pdf(plot_heatmap_ind, width=15, height=15)
heatmap.bp(aa.D.ind)
dev.off()

pop_nj<-nj(aa.D.pop)
pdf(plot_nj_pop)
plot(pop_nj, "u", cex=0.7)
dev.off()

ind_nj<-nj(aa.D.ind)

ind_nj_tips<-ind_nj$tip.label
ind_nj_order<-match(ind_nj_tips, pops$ind)


pdf(plot_nj_ind, width=15, height=15)
plot(ind_nj, "u", cex=0.7, tip.color=pops$colours[ind_nj_order])
dev.off()

aa.splits<-neighborNet(dist(aa.D.pop))

pdf(plot_neighnet_pop)
plot(as.networx(aa.splits), edge.width=1, cex=0.7)
dev.off()

aa.ind.splits<-neighborNet(dist(aa.D.ind))

ind_networx<-as.networx(aa.ind.splits)
ind_networx_tips<-ind_networx$tip.label
ind_networx_order<-match(ind_networx_tips, pops$ind)

pdf(plot_neighnet_ind, width=20, height=12)
plot(ind_networx, tip.color=pops$colours[ind_networx_order], edge.width=0.7, cex=0.5)
dev.off()

#write.phylip()
stamppPhylip(aa.D.ind, out_nei_ind_phylip)

