#!/s/bin/R35

# PUMAS R2
library(data.table)
library(MASS)
library(parallel)
library(matrixStats)
# Read the argument into R
library(optparse)
options(stringsAsFactors=F)
option_list = list(
  make_option("--k", action = "store", default = NA, type = "numeric"),
  make_option("--partitions", action = "store", default = NA, type = "character"),
  make_option("--trait_name", action = "store", default = NA, type = "character"),
  make_option("--gwas_path", action = "store", default = NA, type = "character"),
  make_option("--ld_path", action = "store", default = NA, type = "character"),
  make_option("--output_path", action = "store", default = NA, type = "character"),
  make_option("--chr", action = "store", default = NULL, type = "numeric"),
  make_option("--multicore", action = "store_true", default = FALSE),
  make_option("--cores_number", action = "store", default = detectCores(), type = "numeric")
)

opt = parse_args(OptionParser(option_list=option_list))
k <- opt$k
partitions <- as.numeric(unlist(lapply((strsplit(opt$partitions, ',')),trimws)))
trait_name <- opt$trait_name
gwas_path <- opt$gwas_path
ld_path <- opt$ld_path
output_path <- opt$output_path
chr <- opt$chr
multicore <- opt$multicore
cores_number <- opt$cores_number

if (!multicore) {cores_number = 1}

### functions below

# subsample sumstats based on a subset of individuals from the entire GWAS data
# sample splitting is 60, 20, 10, 10 (Ntr, Nt, Nvtr, Nvt)
PUMAS.II.FUN.I <- function(FunII.GWAS, FunII.Nt, FunII.LD, FunII.LD.block, trait_name, k, chr){
  
  # parameter initialization
  FunII.N.samplesize <- FunII.GWAS$N
  FunII.Ntr <- FunII.N.samplesize - FunII.Nt
  FunII.MAF <- FunII.GWAS$MAF
  FunII.X.VAR <- 2*FunII.MAF*(1-FunII.MAF)
  FunII.SE <- FunII.GWAS$SE
  FunII.beta <- FunII.GWAS$BETA
  
  # statistics initialization
  FunII.Var.Y <- quantile(FunII.SE^2*FunII.N.samplesize*FunII.X.VAR,probs=seq(0,1,0.1))[10]
  FunII.XtY <- FunII.beta*FunII.N.samplesize*FunII.X.VAR
  XtY_mu <- (FunII.Nt/FunII.N.samplesize)*FunII.XtY

  # sampling covariance matrix for XtY_t
  project_mat <- mclapply(1:length(FunII.LD.block), function(j){
      if (FunII.LD.block[[j]][1]==FunII.LD.block[[j]][2]) {
          FunII.SE.tmp <- FunII.SE[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
          FunII.X.VAR.tmp <- FunII.X.VAR[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
          FunII.N.samplesize.tmp <- FunII.N.samplesize[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
          FunII.NminusNv.tmp <- FunII.N.samplesize.tmp - FunII.Nt
          
          epsilon_mat <- FunII.SE.tmp*sqrt(FunII.X.VAR.tmp)*sqrt(FunII.N.samplesize.tmp)
          XtY_sigma <- (median(FunII.NminusNv.tmp)*(median(FunII.N.samplesize.tmp) - median(FunII.NminusNv.tmp))/median(FunII.N.samplesize.tmp))*(epsilon_mat^2)*as.numeric(FunII.LD[[j]])*FunII.X.VAR.tmp
          return(sqrt(XtY_sigma))
      }
      FunII.SE.tmp <- FunII.SE[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
      FunII.X.VAR.tmp <- FunII.X.VAR[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
      FunII.X.VAR.mat.tmp <- sqrt(FunII.X.VAR.tmp %*% t(FunII.X.VAR.tmp))
      FunII.N.samplesize.tmp <- FunII.N.samplesize[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
      FunII.NminusNv.tmp <- FunII.N.samplesize.tmp - FunII.Nt
        
      SE_mat_1 <- matrix(rep(FunII.SE.tmp*sqrt(FunII.X.VAR.tmp)*sqrt(FunII.N.samplesize.tmp),length(FunII.SE.tmp)),length(FunII.SE.tmp),length(FunII.SE.tmp))
      SE_mat_2 <- matrix(rep(FunII.SE.tmp*sqrt(FunII.X.VAR.tmp)*sqrt(FunII.N.samplesize.tmp),each=length(FunII.SE.tmp)),length(FunII.SE.tmp),length(FunII.SE.tmp))
      epsilon_mat <- SE_mat_1
      epsilon_mat[SE_mat_1>=SE_mat_2] <- SE_mat_2[SE_mat_1>=SE_mat_2]
      XtY_sigma <- (median(FunII.NminusNv.tmp)*(median(FunII.N.samplesize.tmp) - median(FunII.NminusNv.tmp))/median(FunII.N.samplesize.tmp))*(epsilon_mat^2)*FunII.LD[[j]]*FunII.X.VAR.mat.tmp
      eigen_decom <- eigen(XtY_sigma)
      eigen_decom_value <- pmax(eigen_decom$value,0)
      project_block <- eigen_decom$vector %*% diag(sqrt(eigen_decom_value))
      return(project_block)
    }, mc.cores=cores_number)

    # sample summary statistics
    FunII.Cor.squared <- c()
    for (i in 1:k){
      if (!is.null(chr)) {
        set.seed(as.numeric(substring(chr,2)))
      } else {
        set.seed(i)
      }
      FunII.XtY.test.temp <- c()
      for (j in 1:length(FunII.LD.block)) {
        print(j)
        XtY_mu_tmp <- XtY_mu[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
        if (FunII.LD.block[[j]][1] == FunII.LD.block[[j]][2]) {
          FunII.XtY.test.temp <- c(FunII.XtY.test.temp,unlist(XtY_mu_tmp + project_mat[[j]] * rnorm(n=1,0,1)))
        } else {
          FunII.XtY.test.temp <- c(FunII.XtY.test.temp,unlist(XtY_mu_tmp + project_mat[[j]] %*% rnorm(ncol(FunII.LD[[j]]),0,1)))
        }
      }

    ## calculate training summary statistics
    FunII.XtY.train.temp <- FunII.XtY - FunII.XtY.test.temp
    FunII.beta.train.temp <- FunII.XtY.train.temp/(FunII.Ntr*FunII.X.VAR)
    FunII.SE.train.temp <- sqrt(FunII.N.samplesize/FunII.Ntr)*FunII.GWAS$SE
    FunII.Zscore.temp <- FunII.XtY.train.temp/(sqrt(FunII.N.samplesize*FunII.Ntr)*FunII.X.VAR*FunII.GWAS$SE)
    FunII.pvalue.temp <- 2*pnorm(abs(FunII.Zscore.temp),lower.tail=F)
        
    ## write new beta, SE, and p-value in the GWAS sumstats
    FunII.GWAS.tmp <- FunII.GWAS
    FunII.GWAS.tmp$`BETA` <- FunII.beta.train.temp
    FunII.GWAS.tmp$`SE` <- FunII.SE.train.temp
    FunII.GWAS.tmp$`P` <- FunII.pvalue.temp
    FunII.GWAS.tmp$`N` <- FunII.Ntr

    fwrite(FunII.GWAS.tmp,paste0(output_path,trait_name,".gwas.ite",i,chr,".txt"),col.names=T,row.names=F,sep="\t",quote=F)  
    
    ## write XtY
    fwrite(data.frame(CHR=FunII.GWAS$CHR,SNP=FunII.GWAS$SNP,A1=FunII.GWAS$A1,A2=FunII.GWAS$A2,test=FunII.XtY.test.temp),paste0(output_path,trait_name,".xty.ite",i,chr,".txt"),col.names=T,row.names=F,sep="\t",quote=F)  
  }
  
  return(FunII.Var.Y)
}


## match order of sumstats, LD genotype, and LD matrices
match_gwas_LD <- function(gwas,LD,rs,bp=NULL){

  match <- Filter(length, mclapply(1:length(LD), function(i){
    snp.overlap <- intersect(gwas$SNP,rs[[i]])
    if (length(snp.overlap)==0) {
      return()
    } else {
      ld_blk <- LD[[i]][match(snp.overlap,rs[[i]]),match(snp.overlap,rs[[i]])]
      rs_blk <- rs[[i]][match(snp.overlap,rs[[i]])]
      gwas <- gwas[match(snp.overlap,gwas$SNP),]
      ld_blk_size <- length(snp.overlap)
      return(list(ld_blk=ld_blk, rs_blk=rs_blk, gwas=gwas, ld_blk_size=ld_blk_size))
    }
  }, mc.cores=cores_number))

  new_ld_blocks <- list()
  new_rs_blocks <- list()
  ld_block_size <- list()
  gwas_matched <- c()
  start.num <- 1
  end.num <- 0
  for (i in 1:length(match)) {
    print(i)
    new_ld_blocks[[i]] <- match[[i]]$ld_blk
    new_rs_blocks[[i]] <- match[[i]]$rs_blk
    gwas_matched <- rbind(gwas_matched, match[[i]]$gwas)
    end.num <- end.num + match[[i]]$ld_blk_size
    ld_block_size[[i]] <- c(start.num, end.num)
    start.num <- end.num + 1
  }
  
  return(list(gwas_matched=gwas_matched,new_rs_blocks=new_rs_blocks,new_ld_blocks=new_ld_blocks,ld_block_size=ld_block_size))
}



# main function
main <- function(){
    gwas <- fread(paste0(gwas_path,trait_name,".txt.gz"),header=T)
    if (!is.null(chr)) {
      gwas <- gwas[gwas$`CHR` == chr]
      chr <- paste0(".", chr)
    }

    load(paste0(ld_path, "/ld_1kg.RData"))
    load(paste0(ld_path, "/rs_1kg.RData"))

    # match GWAS SNPs with LD reference
    matched_data <- match_gwas_LD(gwas=gwas,LD=LD_ref,rs=rs_ref,bp=NULL)
    rm(LD_ref)
    rm(rs_ref)
    rm(gwas)

    # PUMAS R2, k times
    pumas_tmp <- PUMAS.II.FUN.I(FunII.GWAS=matched_data$gwas_matched, FunII.Nt=floor(partitions[2]*min(matched_data$gwas_matched$N)),FunII.LD=matched_data$new_ld_blocks, FunII.LD.block=matched_data$ld_block_size, trait_name=trait_name, k=k, chr=chr)

    fwrite(data.frame(var.Y=pumas_tmp,N.t=floor(partitions[2]*min(matched_data$gwas_matched$N))),paste0(output_path,trait_name,".forEVAL",chr,".txt"),col.names=T,row.names=F,sep="\t",quote=F)
}

## execute main function
main()
