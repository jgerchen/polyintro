#install.packages("StAMPP")
library(StAMPP)
library(vcfR)
library(adegenet)
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
vcf_input<-("Desktop/polyintro/popgen_workflow/workflow/scripts/cardamine_ld_pruned.vcf")

vcf <- read.vcfR(vcf_input)
aa.genlight <- vcfR2genlight.tetra(vcf)
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_") # add real SNP.names
pop(aa.genlight)<-substr(indNames(aa.genlight),1,5)  # add pop names
aa.D.pop <- stamppNeisD(aa.genlight, pop = TRUE)   # Nei's 1972 genetic distance between pops
colnames(aa.D.pop)<-rownames(aa.D.pop)
heatmap(aa.D.pop, Rowv = NA, Colv=NA)
plot(nj(aa.D.pop))


