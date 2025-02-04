
library(adegenet)
library(vcfR)

args = commandArgs(trailingOnly=TRUE)
vcf_input<-args[1]
pop_map<-args[2]
pca_PCs_out<-args[3]
pca_Eigenvalues_out<-args[4]
pca_plot<-args[5]
#Adegenet functions taken from https://github.com/mbohutinska/PolyplChapter and extended to work with triploids, pentaploids and hexaploids
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

#######################


## -------------------------------------------------------
### a patch for MUCH MUCH faster PCA calculation on genlight objects
# see https://github.com/thibautjombart/adegenet/pull/150
glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  
  if(!inherits(x, "genlight")) stop("x is not a genlight object")
  # keep the original mean / var code, as it's used further down
  # and has some NA checks..
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  # convert to full data, try to keep the NA handling as similar
  # to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)
  # all dot products at once using underlying BLAS
  # to support thousands of samples, this could be
  # replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights
  ## PERFORM THE ANALYSIS ##
  ## eigenanalysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
  ## scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  ## rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging
  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  ## GET LOADINGS ##
  ## need to decompose X^TDV into a sum of n matrices of dim p*r
  ## but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }
    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }
  ## FORMAT OUTPUT ##
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}


dip_colour<-"#1e90ff"
tet_colour<-"#ffa500"
tri_colour<-"#DB3192"
#no 5x colour yet
pent_colour<-"#bcdb31"
hex_colour<-"#7cfc00"
# ---------------------------------------------------------
#read VCF
vcf <- read.vcfR(vcf_input)
pops<-read.table(pop_map, sep="\t", header=FALSE, col.names=c("ind", "pop", "ploidy"))
aa.genlight <- vcfR2genlight.tetra(vcf)
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_") # add real SNP.names
indnames<-indNames(aa.genlight)
#popnames<-rep(NA, length(indnames))
#for(i in 1:length(indnames)){
#  popnames[i]<-pops$pop[pops$ind==indnames[i]]
#}
pop_order<-match(indnames,pops$ind)
#warnings()
pop(aa.genlight)<-pops$pop[pop_order]
#make PCA
pca.1 <- glPcaFast(aa.genlight, nf=300)
#write PCA results
write.table(pca.1$scores,file=pca_PCs_out ,quote=FALSE)
write.table(pca.1$eig,file=pca_Eigenvalues_out,quote=FALSE)
#make plot
col_vector<-pops$ploidy
col_vector[col_vector==2]<-dip_colour
col_vector[col_vector==3]<-tri_colour
col_vector[col_vector==4]<-tet_colour
col_vector[col_vector==5]<-pent_colour
col_vector[col_vector==6]<-hex_colour

pdf (pca_plot, width=12, height=8)
par(mar=c(5.1, 4.1, 0.5, 8.1), xpd=TRUE)
plot(pca.1$scores, xlab=paste(round((pca.1$eig[1]/sum(pca.1$eig)*100),1)," %"), ylab=paste(round((pca.1$eig[2]/sum(pca.1$eig)*100),1)," %"), col=col_vector[pop_order], pch=as.numeric(pop(aa.genlight)))
legend("topright",legend=unique(pop(aa.genlight)), pch=unique(as.numeric(pop(aa.genlight))), inset=c(-0.15,0), title="Pop")
legend("bottomright",legend=unique(pops$ploidy), col=unique(col_vector) , pch=19, inset=c(-0.1,0), title="Ploidy")
dev.off()



