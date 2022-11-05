#!/s/bin/R35

if(!require(data.table)){
    install.packages("data.table")
    library(data.table)
}
if(!require(BEDMatrix)){
    install.packages("BEDMatrix")
    library(BEDMatrix)
}
if(!require(optparse)){
    install.packages("optparse")
    library(optparse)
}

# Read the argument into R
options(stringsAsFactors=F)
option_list = list(
  make_option("--k", action = "store", default = NA, type = "numeric"),
  make_option("--ref_path", action = "store", default = NA, type = "character"),
  make_option("--trait_name", action = "store", default = NA, type = "character"),
  make_option("--prs_method", action = "store", default = NA, type = "character"),
  make_option("--xty_path", action = "store", default = NA, type = "character"),
  make_option("--stats_path", action = "store", default = NA, type = "character"),
  make_option("--weight_path", action = "store", default = NA, type = "character"),
  make_option("--output_path", action = "store", default = NA, type = "character")
)

opt = parse_args(OptionParser(option_list=option_list))
k <- opt$k
ref_path <- opt$ref_path
trait_name <- opt$trait_name
prs_method <- opt$prs_method
xty_path <- opt$xty_path
stats_path <- opt$stats_path
weight_path <- opt$weight_path
output_path <- opt$output_path

### functions

## calculate sumstats-based R2
get_sumR2 <- function(XtY.t, weight, var.Y, N.t, X.ref){
  weight <- t((t(weight) - colMeans(weight)) / apply(weight, 2, sd)) # idk why we cannot use colSds (from matrixStats)
  cov.Y_Y.Hat <- colSums(weight*(XtY.t/N.t))
  Y.Hat <- X.ref %*% weight
  var.Y.Hat <- apply(Y.Hat, 2, var)
  sum.R2 <- (unlist(cov.Y_Y.Hat))^2/(var.Y*unlist(var.Y.Hat))
  sum.R2[is.na(sum.R2)] <-0
  return(sum.R2)
}


## main function
main <- function(){
    snp.w <- as.data.frame(fread(paste0(weight_path,trait_name,".",prs_method,".ite1.txt"),header=T))
    suffix <- colnames(snp.w)[-c(1:4)]
    xty <- as.data.frame(fread(paste0(xty_path,trait_name,".xty.ite1.txt"),header=T))
    stats.pumas <- as.data.frame(fread(paste0(stats_path,trait_name,".forEVAL.txt"),header=T))
    
    ## add a function here to coordinate A1/A2 between XtY and SNP weights data ##
    
    ## end ##
    
    geno_ref <- BEDMatrix(ref_path)
    rs.geno <- fread(paste0(ref_path,".bim"),header=F)$V2
    geno_ref <- geno_ref[,match(xty$SNP,rs.geno)]
    pumas.cor2 <- c()
    for (j in 1:k) {
        xty <- fread(paste0(xty_path,trait_name,".xty.ite",j,".txt"),header=T)
        snp.w <- as.data.frame(fread(paste0(weight_path,trait_name,".",prs_method,".ite",j,".txt"),header=T))
        snp.weight <- snp.w[, suffix]
        snp.weight[snp.weight==Inf|snp.weight==-Inf] <- 0
        snp.weight[is.na(snp.weight)] <- 0
        pumas.cor2.tmp <- get_sumR2(XtY.t=xty$test, weight=snp.weight, var.Y=stats.pumas$var.Y, N.t=stats.pumas$N.t, X.ref=geno_ref)
        pumas.cor2 <- rbind(pumas.cor2,pumas.cor2.tmp)
    }
    colnames(pumas.cor2) <- suffix
    # store R2
    write.table(pumas.cor2,paste0(output_path,trait_name,".",prs_method,".txt"),col.names = T,row.names=F,quote=F,sep="\t")
}

main()
