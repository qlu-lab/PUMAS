pumas.main<-function(input_path,output_path, beta_header,maf_header,se_header, sample_size,flag_plot){
  data.real <- read.table(paste0(input_path), header=T,na.strings=c("na","NA",""))
  data.real=na.omit(data.real)
  if (!is.character(beta_header)){stop(print("beta header should be a character"))}
  if (!is.character(se_header)){stop(print("se header should be a character"))}
  #TODO to create a MAF
  if (is.null(maf_header)){print("MAF processing with reference")}
  #If usr puts in a number for sample size
  #if(is.numeric(sample_size)){N.sample = sample_size}
  beta=data.real[,which(names(data.real)==beta_header)]
  MAF=data.real[,which(names(data.real)==maf_header)]
  se=data.real[,which(names(data.real)==se_header)]
  X_VAR=2*MAF*(1-MAF)
  N.sample = sample_size
  
  p.value=data.real$Pval
  p.value=p.value[order(p.value,decreasing=F)]
  
  Res_TildeRL <- FunII.TildeRL(FunII.beta=beta,FunII.SE=se,FunII.sigma=X_VAR,FunII.N.samplesize=N.sample,FunII.Nv=0.25*N.sample,FunII.rep=4, FunII.Corr=F)
  
  if(flag_plot){FunI.plot(Res_TildeRL=Res_TildeRL,p.value=p.value)}
  
  
}

FunI.plot<-function(Res_TildeRL,p.value,output_path){
  png(paste0(output_path),units="in",width = 9,height = 5,res = 300)
  plot(y=Res_TildeRL,x=log(p.value,base=10),main=paste0("T0030_pruned"),ylab="COR^2",xlab="log(P-value)")
  points(y=max(Res_TildeRL),x=log(p.value,base=10)[which(Res_TildeRL==max(Res_TildeRL))],col="red", pch=19)
  abline(v=log(p.value,base=10)[which(Res_TildeRL==max(Res_TildeRL))],col="red",lty=2)
  text(labels=paste0(p.value[which(Res_TildeRL==max(Res_TildeRL))]),y=0,x=log(p.value,base=10)[which(Res_TildeRL==max(Res_TildeRL))],pos=3,col="red")
 # legend("p-value is : ", col="red", cex=0.8)
  dev.off()
}
### ---  Funtion File ---###

#--- Prediction Error: --->
#------ the paramters should be ordered in advance
FunI.PEhat<-function(FunI.beta_hat,FunI.test.Y,FunI.test.X){
  if (! is.vector(FunI.beta_hat)){stop(print("beta_hat should be a vector!"))}
  Y_hat=matrix(NA,nrow=dim(FunI.test.X)[1],ncol=dim(FunI.test.X)[2])
  for (ii in 1:dim(FunI.test.X)[1]){
    Y_hat[ii,]=cumsum(FunI.test.X[ii,]*FunI.beta_hat)
  }
  FunII.R2=c()
  for (jj in 1:dim(Y_hat)[2]){
    FunII.R2[jj] = (cor(Y_hat[,jj],FunI.test.Y))^2
  }
  return(FunII.R2)
}


#------ the paramters should be ordered in advance, including FunI.Sigma:
#---------  *1. FunI.sigma is the true covariance matrix for (x_1,...,x_m), not necessary a diagnal matrix;
FunI.PEtilde<-function(FunI.XtY.tr, FunI.XtY.v, FunI.ntr, FunI.nv, FunI.Sigma,FunI.Var.Y){
  if (! is.vector(FunI.XtY.tr)){stop(print("X_t_Y.tr should be a vector!"))}
  if (! is.vector(FunI.XtY.v)){stop(print("X_t_Y.v should be a vector!"))}
  if (length(FunI.XtY.tr)!=length(FunI.XtY.v)){stop(print("Please check dimensions of parameters!"))}
  
  if(is.vector(FunI.Sigma)){
    FunI.beta_tilde <- (FunI.XtY.tr/FunI.ntr)/(FunI.Sigma)
    FunI.PE.term3 <- FunI.Sigma*(FunI.beta_tilde^2)
    FunI.Var.Y.Hat = cumsum(FunI.PE.term3)
    FunI.Cov.Y_Y.Hat = cumsum(FunI.beta_tilde*(FunI.XtY.v/FunI.nv))
    FunI.Cor.squared = FunI.Cov.Y_Y.Hat^2/(FunI.Var.Y.Hat*FunI.Var.Y)
  }else{
    FunI.Sigma.diag <- diag(FunI.Sigma)
    FunI.beta_tilde <- (FunI.XtY.tr/FunI.ntr)*(FunI.Sigma.diag^-1)
    FunI.PE.term2 <- -2*FunI.beta_tilde*(FunI.XtY.v/FunI.nv)
    FunI.PE.term3 <- apply(FunI.Sigma*FunI.beta_tilde,1,function(z){z*FunI.beta_tilde})
    FunI.PE <- cumsum(FunI.PE.term2) + diag(apply(apply(FunI.PE.term3,2,cumsum),1,cumsum))
  }
}

###--- Tilde CV ---###
#*** 1. This v4 version is for Repteaded Learning only;
#*** 2. if Corr=F (default), VarXY is approximated in the simplest way, which is a diagonal matrix;
#***    if Corr=T, VarXY is a m-dim matrix.
#*** 3. FunII.beta & FunII.SE are m-dim vectors;
#***    FunII.sigma is LD matrix with E(X_j^2) as diagonal element, X_j's are centered; When FunII.corr=F, FunII.sigma should be a vector;
#***    FunII.N.samplesize is total sample size for every SumStat, can be either a number or vector;
#***    FunII.Ntr is a number for how large sample size is used to train, which is smaller than max(FunII.N.samplesize)
#***    FunII.rep is a number for repeated-learning.
FunII.TildeRL<-function(FunII.beta, FunII.SE, FunII.sigma, FunII.N.samplesize, FunII.Nv, FunII.rep, FunII.Corr=F){
  if(length(FunII.N.samplesize)==1){
    FunII.N.samplesize <- rep(FunII.N.samplesize, length(FunII.beta))
  }else{
    if(length(FunII.N.samplesize)!=length(FunII.beta)){
      stop("Error: Length of FunII.N is not the same as FunII.beta!")
    }
  }
  if(FunII.Nv>=max(FunII.N.samplesize)){
    print(FunII.Nv)
    print(max(FunII.N.samplesize))
    stop("Error: Nv is larger or equal than the maximum of N.samplesize!")}
  FunII.Ntr <- FunII.N.samplesize-FunII.Nv
  
  if(is.vector(FunII.sigma)){
    FunII.sigma.diag <- FunII.sigma
  }else{
    FunII.sigma.diag <- diag(FunII.sigma)
  }
  
  FunII.XtY <- FunII.beta*FunII.N.samplesize*FunII.sigma.diag
  FunII.Var.Y=quantile(FunII.SE^2*FunII.N.samplesize*FunII.sigma.diag,probs=seq(0,1,0.1))[10]
  FunII.PErr.Avg <- rep(0, length(FunII.beta))
  
  if (FunII.Corr==F){
    # If Corr=F, FunII.VarXY is a vector
    FunII.VarXY <- FunII.N.samplesize*(FunII.SE^2)*(FunII.sigma.diag^2)
    for (FunII.r in 1:FunII.rep){
      FunII.XtY.train.temp <- rnorm(n=length(FunII.beta), mean = (FunII.Ntr/FunII.N.samplesize)*FunII.XtY, sd = sqrt((FunII.Ntr*FunII.Nv/FunII.N.samplesize)*FunII.VarXY))
      FunII.XtY.test.temp <- FunII.XtY - FunII.XtY.train.temp
      FunII.Zscore.temp <- FunII.XtY.train.temp/(sqrt(FunII.N.samplesize*FunII.Ntr)*FunII.SE*FunII.sigma.diag)
      FunII.order.temp <- order(abs(FunII.Zscore.temp),decreasing=T)
      if(is.vector(FunII.sigma)){
        FunII.PErr.Avg <- FunII.PErr.Avg + FunI.PEtilde(FunI.XtY.tr=FunII.XtY.train.temp[FunII.order.temp],
                                                        FunI.XtY.v=FunII.XtY.test.temp[FunII.order.temp],
                                                        FunI.ntr=FunII.Ntr[FunII.order.temp], FunI.nv=FunII.Nv,
                                                        FunI.Sigma=FunII.sigma.diag[FunII.order.temp],FunI.Var.Y=FunII.Var.Y)/FunII.rep
      }else{
        FunII.PErr.Avg <- FunII.PErr.Avg + FunI.PEtilde(FunI.XtY.tr=FunII.XtY.train.temp[FunII.order.temp],
                                                        FunI.XtY.v=FunII.XtY.test.temp[FunII.order.temp],
                                                        FunI.ntr=FunII.Ntr, FunI.nv=FunII.Nv[FunII.order.temp],
                                                        FunI.Sigma=FunII.sigma[FunII.order.temp,FunII.order.temp])/FunII.rep 
      }
    }
  }else{
    # If Corr=T, FunII.VarXY is a matrix
    FunII.VarXY <- sapply(FunII.beta*FunII.sigma.diag, function(x_j){x_j*FunII.beta*FunII.sigma.diag}) - diag((FunII.beta*FunII.sigma.diag)^2) + diag(FunII.N.samplesize*(FunII.SE^2)*(FunII.sigma.diag^2))
    FunII.XtY.train.temp <- FunII.Ntr*mvrnorm(FunII.rep, mu=(1/FunII.N.samplesize)*FunII.XtY, Sigma= (FunII.Nv/(FunII.N.samplesize*FunII.Ntr))*FunII.VarXY)
    FunII.XtY.test.temp <- -sweep(FunII.XtY.train.temp, 2, FunII.XtY)
    FunII.Zscore.temp <- sweep(FunII.XtY.train.temp, 2, sqrt(FunII.N.samplesize*FunII.Ntr)*FunII.SE*FunII.sigma.diag , FUN = "/") # original version
    FunII.order.temp <- apply(abs(FunII.Zscore.temp), 1, function(s){order(s,decreasing=T)})
    if(is.vector(FunII.sigma)){
      #for (FunII.r in 1:FunII.rep){
      for (FunII.r in 1){
        FunII.PErr.Avg <- FunI.PEtilde(FunI.XtY.tr=FunII.XtY.train.temp[FunII.r, FunII.order.temp[,FunII.r]],
                                       FunI.XtY.v=FunII.XtY.test.temp[FunII.r, FunII.order.temp[,FunII.r]],
                                       FunI.ntr=FunII.Ntr, FunI.nv=FunII.Nv[FunII.order.temp[,FunII.r]],
                                       FunI.Sigma=FunII.sigma.diag[FunII.order.temp[,FunII.r]])
      }
    }else{
      for (FunII.r in 1:FunII.rep){
        FunII.PErr.Avg <- FunII.PErr.Avg + FunI.PEtilde(FunI.XtY.tr=FunII.XtY.train.temp[FunII.r, FunII.order.temp[,FunII.r]],
                                                        FunI.XtY.v=FunII.XtY.test.temp[FunII.r, FunII.order.temp[,FunII.r]],
                                                        FunI.ntr=FunII.Ntr, FunI.nv=FunII.Nv[FunII.order.temp[,FunII.r]],
                                                        FunI.Sigma=FunII.sigma[FunII.order.temp[,FunII.r],FunII.order.temp[,FunII.r]])/FunII.rep
      }
    }
  }
  return(FunII.PErr.Avg)
}







