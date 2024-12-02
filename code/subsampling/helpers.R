if(!require(data.table)){
    install.packages("data.table")
    library(data.table)
}
if(!require(parallel)){
    install.packages("parallel")
    library(parallel)
}
if(!require(optparse)){
    install.packages("optparse")
    library(optparse)
}

# PUMAS SUBSAMPLING
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
    }, mc.cores = threads)

    # sample summary statistics
    FunII.Cor.squared <- c()
    subsample <- mclapply(1:k, function(i){
      if (!is.null(chr)) {
        set.seed(as.numeric(substring(chr,2)))
      } else {
        set.seed(i)
      }
      FunII.XtY.test.temp <- c()
      for (j in 1:length(FunII.LD.block)) {
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
  }, mc.cores = threads)
  
  return(FunII.Var.Y)
}

## PUMA-CUBS SUBSAMPLING
# subsample sumstats based on a subset of individuals from the entire GWAS data
# sample splitting is 60, 20, 10, 10 (Ntr, Nt, Nvtr, Nvt)
PUMAS.II.FUN.II <- function(FunII.GWAS, FunII.Nt, FunII.Nvtr, FunII.Nvt, FunII.LD, FunII.LD.block, trait_name, k, chr){
  
  # parameter initialization
  FunII.N.samplesize <- FunII.GWAS$N
  FunII.Ntr <- FunII.N.samplesize - FunII.Nt - FunII.Nvtr - FunII.Nvt
  FunII.MAF <- FunII.GWAS$MAF
  FunII.X.VAR <- 2*FunII.MAF*(1-FunII.MAF)
  FunII.SE <- FunII.GWAS$SE
  FunII.beta <- FunII.GWAS$BETA
  
  # statistics initialization
  FunII.Var.Y <- quantile(FunII.SE^2*FunII.N.samplesize*FunII.X.VAR,probs=seq(0,1,0.1))[10]
  FunII.XtY <- FunII.beta*FunII.N.samplesize*FunII.X.VAR
  XtY_mu <- (FunII.Nt/FunII.N.samplesize)*FunII.XtY
  
  # sampling covariance matrix for XtY_t
  cat("\n\nGenerating sampling covariance martrix ...\n")
  begin.time = start_time("Generating sampling covariance martrix")
  FunII.SE1 <- FunII.SE
  FunII.SE2 <- FunII.GWAS$SE * sqrt(FunII.N.samplesize/(FunII.N.samplesize - FunII.Nt)) #added
  FunII.SE3 <- FunII.GWAS$SE * sqrt(FunII.N.samplesize/(FunII.N.samplesize - FunII.Nt - FunII.Nvtr)) #added
  project_mat <- mclapply(1:length(FunII.LD.block), function(j){
      if (FunII.LD.block[[j]][1]==FunII.LD.block[[j]][2]) {
          FunII.X.VAR.tmp <- FunII.X.VAR[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]

          FunII.SE.tmp <- FunII.SE1[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
          FunII.N.samplesize.tmp <- FunII.N.samplesize[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
          FunII.NminusNv.tmp <- FunII.N.samplesize.tmp - FunII.Nt
          epsilon_mat <- FunII.SE.tmp*sqrt(FunII.X.VAR.tmp)*sqrt(FunII.N.samplesize.tmp)
          XtY_sigma <- (median(FunII.NminusNv.tmp)*(median(FunII.N.samplesize.tmp) - median(FunII.NminusNv.tmp))/median(FunII.N.samplesize.tmp))*(epsilon_mat^2)*as.numeric(FunII.LD[[j]])*FunII.X.VAR.tmp
          project_mat1 <- sqrt(XtY_sigma)

          FunII.SE.tmp <- FunII.SE2[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
          FunII.N.samplesize.tmp <- FunII.NminusNv.tmp
          FunII.NminusNv.tmp <- FunII.N.samplesize.tmp - FunII.Nvtr
          epsilon_mat <- FunII.SE.tmp*sqrt(FunII.X.VAR.tmp)*sqrt(FunII.N.samplesize.tmp)
          XtY_sigma <- (median(FunII.NminusNv.tmp)*(median(FunII.N.samplesize.tmp) - median(FunII.NminusNv.tmp))/median(FunII.N.samplesize.tmp))*(epsilon_mat^2)*as.numeric(FunII.LD[[j]])*FunII.X.VAR.tmp
          project_mat2 <- sqrt(XtY_sigma)

          FunII.SE.tmp <- FunII.SE3[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
          FunII.N.samplesize.tmp <- FunII.NminusNv.tmp
          FunII.NminusNv.tmp <- FunII.N.samplesize.tmp - FunII.Nvt
          epsilon_mat <- FunII.SE.tmp*sqrt(FunII.X.VAR.tmp)*sqrt(FunII.N.samplesize.tmp)
          XtY_sigma <- (median(FunII.NminusNv.tmp)*(median(FunII.N.samplesize.tmp) - median(FunII.NminusNv.tmp))/median(FunII.N.samplesize.tmp))*(epsilon_mat^2)*as.numeric(FunII.LD[[j]])*FunII.X.VAR.tmp
          project_mat3 <- sqrt(XtY_sigma)
        return(list(project_mat1=project_mat1, project_mat2=project_mat2, project_mat3=project_mat3))
      }
      FunII.X.VAR.tmp <- FunII.X.VAR[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
      FunII.X.VAR.mat.tmp <- sqrt(FunII.X.VAR.tmp %*% t(FunII.X.VAR.tmp))

      FunII.SE.tmp <- FunII.SE1[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
      FunII.N.samplesize.tmp <- FunII.N.samplesize[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
      FunII.NminusNv.tmp <- FunII.N.samplesize.tmp - FunII.Nt
      SE_mat_1 <- matrix(rep(FunII.SE.tmp*sqrt(FunII.X.VAR.tmp)*sqrt(FunII.N.samplesize.tmp),length(FunII.SE.tmp)),length(FunII.SE.tmp),length(FunII.SE.tmp))
      SE_mat_2 <- matrix(rep(FunII.SE.tmp*sqrt(FunII.X.VAR.tmp)*sqrt(FunII.N.samplesize.tmp),each=length(FunII.SE.tmp)),length(FunII.SE.tmp),length(FunII.SE.tmp))
      epsilon_mat <- SE_mat_1
      epsilon_mat[SE_mat_1>=SE_mat_2] <- SE_mat_2[SE_mat_1>=SE_mat_2]
      XtY_sigma <- (median(FunII.NminusNv.tmp)*(median(FunII.N.samplesize.tmp) - median(FunII.NminusNv.tmp))/median(FunII.N.samplesize.tmp))*(epsilon_mat^2)*FunII.LD[[j]]*FunII.X.VAR.mat.tmp
      #XtY_sigma_diag <- diag(diag(XtY_sigma))
      eigen_decom <- eigen(XtY_sigma)
      #print(eigen_decom$value)
      eigen_decom_value <- pmax(eigen_decom$value,0)
      project_mat1 <- eigen_decom$vector %*% diag(sqrt(eigen_decom_value))

      FunII.SE.tmp <- FunII.SE2[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
      FunII.N.samplesize.tmp <- FunII.NminusNv.tmp
      FunII.NminusNv.tmp <- FunII.N.samplesize.tmp - FunII.Nvtr
      SE_mat_1 <- matrix(rep(FunII.SE.tmp*sqrt(FunII.X.VAR.tmp)*sqrt(FunII.N.samplesize.tmp),length(FunII.SE.tmp)),length(FunII.SE.tmp),length(FunII.SE.tmp))
      SE_mat_2 <- matrix(rep(FunII.SE.tmp*sqrt(FunII.X.VAR.tmp)*sqrt(FunII.N.samplesize.tmp),each=length(FunII.SE.tmp)),length(FunII.SE.tmp),length(FunII.SE.tmp))
      epsilon_mat <- SE_mat_1
      epsilon_mat[SE_mat_1>=SE_mat_2] <- SE_mat_2[SE_mat_1>=SE_mat_2]
      XtY_sigma <- (median(FunII.NminusNv.tmp)*(median(FunII.N.samplesize.tmp) - median(FunII.NminusNv.tmp))/median(FunII.N.samplesize.tmp))*(epsilon_mat^2)*FunII.LD[[j]]*FunII.X.VAR.mat.tmp
      #XtY_sigma_diag <- diag(diag(XtY_sigma))
      eigen_decom <- eigen(XtY_sigma)
      #print(eigen_decom$value)
      eigen_decom_value <- pmax(eigen_decom$value,0)
      project_mat2 <- eigen_decom$vector %*% diag(sqrt(eigen_decom_value))

      FunII.SE.tmp <- FunII.SE3[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
      FunII.N.samplesize.tmp <- FunII.NminusNv.tmp
      FunII.NminusNv.tmp <- FunII.N.samplesize.tmp - FunII.Nvt
      SE_mat_1 <- matrix(rep(FunII.SE.tmp*sqrt(FunII.X.VAR.tmp)*sqrt(FunII.N.samplesize.tmp),length(FunII.SE.tmp)),length(FunII.SE.tmp),length(FunII.SE.tmp))
      SE_mat_2 <- matrix(rep(FunII.SE.tmp*sqrt(FunII.X.VAR.tmp)*sqrt(FunII.N.samplesize.tmp),each=length(FunII.SE.tmp)),length(FunII.SE.tmp),length(FunII.SE.tmp))
      epsilon_mat <- SE_mat_1
      epsilon_mat[SE_mat_1>=SE_mat_2] <- SE_mat_2[SE_mat_1>=SE_mat_2]
      XtY_sigma <- (median(FunII.NminusNv.tmp)*(median(FunII.N.samplesize.tmp) - median(FunII.NminusNv.tmp))/median(FunII.N.samplesize.tmp))*(epsilon_mat^2)*FunII.LD[[j]]*FunII.X.VAR.mat.tmp
      #XtY_sigma_diag <- diag(diag(XtY_sigma))
      eigen_decom <- eigen(XtY_sigma)
      #print(eigen_decom$value)
      eigen_decom_value <- pmax(eigen_decom$value,0)
      project_mat3 <- eigen_decom$vector %*% diag(sqrt(eigen_decom_value))
      return(list(project_mat1=project_mat1, project_mat2=project_mat2, project_mat3=project_mat3))
  }, mc.cores = threads)
  end_time("Generating sampling covariance martrix", begin.time)

  project_mat1 <- list()
  project_mat2 <- list()
  project_mat3 <- list()
  for (j in 1:length(project_mat)){
    project_mat1[[j]] <- project_mat[[j]]$project_mat1
    project_mat2[[j]] <- project_mat[[j]]$project_mat2
    project_mat3[[j]] <- project_mat[[j]]$project_mat3
  }
  
  FunII.Cor.squared <- c()
  cat("\n\nRunning subsampling over", k, "iterations ...\n")
  begin.time = start_time("Subsampling")
  subsampling <- mclapply(1:k, function(i){
    if (!is.null(chr)) {
      set.seed(as.numeric(substring(chr,2)))
    } else {
      set.seed(i)
    }
    FunII.XtY.test.temp <- c()
    FunII.XtY.vali.train.temp <- c()
    FunII.XtY.vali.test.temp <- c()
    for (j in 1:length(FunII.LD.block)) {
        XtY_mu_tmp <- XtY_mu[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
        if (FunII.LD.block[[j]][1] == FunII.LD.block[[j]][2]) {
          FunII.XtY.test.temp <- c(FunII.XtY.test.temp,unlist(XtY_mu_tmp + project_mat1[[j]] * rnorm(n=1,0,1)))
        } else {
          FunII.XtY.test.temp <- c(FunII.XtY.test.temp,unlist(XtY_mu_tmp + project_mat1[[j]] %*% rnorm(ncol(FunII.LD[[j]]),0,1)))
        }
    }

    FunII.XtY.temp1 <- FunII.XtY - FunII.XtY.test.temp
    XtY_mu1 <- FunII.XtY.temp1*(FunII.Nvtr)/(FunII.N.samplesize - FunII.Nt)
    for (j in 1:length(FunII.LD.block)) {
      XtY_mu_tmp1 <- XtY_mu1[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
      if (FunII.LD.block[[j]][1] == FunII.LD.block[[j]][2]) {
        FunII.XtY.vali.train.temp <- c(FunII.XtY.vali.train.temp,unlist(XtY_mu_tmp1 + project_mat2[[j]] * rnorm(n=1,0,1)))
      } else {
        FunII.XtY.vali.train.temp <- c(FunII.XtY.vali.train.temp,unlist(XtY_mu_tmp1 + project_mat2[[j]] %*% rnorm(ncol(FunII.LD[[j]]),0,1)))
      }
    }
    
    FunII.XtY.temp2 <- FunII.XtY.temp1 - FunII.XtY.vali.train.temp
    XtY_mu2 <- FunII.XtY.temp2*(FunII.Nvt)/(FunII.N.samplesize - FunII.Nt - FunII.Nvtr)
    for (j in 1:length(FunII.LD.block)) {
      XtY_mu_tmp2 <- XtY_mu2[FunII.LD.block[[j]][1]:FunII.LD.block[[j]][2]]
      if (FunII.LD.block[[j]][1] == FunII.LD.block[[j]][2]) {
        FunII.XtY.vali.test.temp <- c(FunII.XtY.vali.test.temp,unlist(XtY_mu_tmp2 + project_mat3[[j]] * rnorm(n=1,0,1)))
      } else {
        FunII.XtY.vali.test.temp <- c(FunII.XtY.vali.test.temp,unlist(XtY_mu_tmp2 + project_mat3[[j]] %*% rnorm(ncol(FunII.LD[[j]]),0,1)))
      }
    }

    ## calculate training summary statistics
    FunII.XtY.train.temp <- FunII.XtY.temp2 - FunII.XtY.vali.test.temp
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
    fwrite(FunII.GWAS.tmp,paste0(output_path,trait_name,".gwas.omnibus.ite",i,chr,".txt"),col.names=T,row.names=F,sep="\t",quote=F)
    
    ## write XtY's
    fwrite(data.frame(CHR=FunII.GWAS$CHR,SNP=FunII.GWAS$SNP,A1=FunII.GWAS$A1,A2=FunII.GWAS$A2,test=FunII.XtY.test.temp,validation_train=FunII.XtY.vali.train.temp,validation_test=FunII.XtY.vali.test.temp),paste0(output_path,trait_name,".xty.omnibus.ite",i,chr,".txt"),col.names=T,row.names=F,sep="\t",quote=F)
    cat("\nSubsampling finished for iteration:", i, "\n")
    ## if running SL models, do subsampling for the innerloop of SL
    if (ensemble == "SL" || ensemble == "all") {
        cat("\nStarting subsampling for innerloop of Superlearning for iteration:", i, "\n")
        # make sumstats for inner loop subsampling
        sl.gwas.tmp <- FunII.GWAS
        sl.gwas.tmp$SE <- sl.gwas.tmp$SE * sqrt(sl.gwas.tmp$N/FunII.Nt)
        sl.gwas.tmp$N <- FunII.Nt
        sl.gwas.tmp$BETA <- FunII.XtY.test.temp/(sl.gwas.tmp$N * 2 * (1 - sl.gwas.tmp$MAF) * sl.gwas.tmp$MAF)
        sl.gwas.tmp$P <- 2*pnorm(abs(sl.gwas.tmp$BETA/sl.gwas.tmp$SE),lower.tail=F)
        # write data with innerloop gwas data
        fwrite(sl.gwas.tmp,paste0(output_path,trait_name,".innerloop.ite",i,chr,".txt.gz"),col.names=T,row.names=F,sep="\t",quote=F)

        # run subsample method with 0.5/0.5 split of innerloop data
        innerloop_trait_name <- paste0(trait_name,".innerloop.ite",i,chr)
        PUMAS.II.FUN.I(sl.gwas.tmp, FunII.Nt=0.5*min(sl.gwas.tmp$N), FunII.LD, FunII.LD.block, innerloop_trait_name, k, chr)
        cat("\nInnerloop of Superlearning subsampling finished for iteration:", i, "\n")
    }
  }, mc.cores = threads)
  end_time("Subsampling", begin.time)
  return(FunII.Var.Y)
}


## match order of sumstats, LD genotype, and LD matrices
match_gwas_LD <- function(gwas,LD,rs,bp=NULL){

  match <- Filter(length, mclapply(1:length(LD), function(i){
    if (!is.matrix(LD[[i]]) & !is.list(LD[[i]])) {
      return()
    }
    snp.overlap <- intersect(gwas$SNP,rs[[i]])
    if (length(snp.overlap)==0) {
      return()
    } else {
      ld_blk <- LD[[i]][match(snp.overlap,rs[[i]]),match(snp.overlap,rs[[i]])]
      rs_blk <- rs[[i]][match(snp.overlap,rs[[i]])]
      gwas <- gwas[match(snp.overlap,gwas$SNP),]
      ld_blk_size <- length(snp.overlap)
      return(list(ld_blk=ld_blk, rs_blk=rs_blk, gwas=gwas, ld_blk_size=ld_blk_size, i=i))
    }
  }, mc.cores = threads))

  new_ld_blocks <- list()
  new_rs_blocks <- list()
  ld_block_size <- list()
  gwas_matched <- c()
  start.num <- 1
  end.num <- 0
  for (i in 1:length(match)) {
    new_ld_blocks[[i]] <- match[[i]]$ld_blk
    new_rs_blocks[[i]] <- match[[i]]$rs_blk
    gwas_matched <- rbind(gwas_matched, match[[i]]$gwas)
    end.num <- end.num + match[[i]]$ld_blk_size
    ld_block_size[[i]] <- c(start.num, end.num)
    start.num <- end.num + 1
  }
  
  return(list(gwas_matched=gwas_matched,new_rs_blocks=new_rs_blocks,new_ld_blocks=new_ld_blocks,ld_block_size=ld_block_size))
}