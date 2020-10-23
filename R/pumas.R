#--- libraries: --->
#library(data.table)
#library(R.utils)
#library(ggplot2)

#--- main function: --->
#' @export
pumas.main<-function(input_path,output_path=NULL,beta_header,af_header,se_header,pvalue_header,samplesize_header,n_fold=NULL,odds_ratio=NULL,make_plot=NULL){
  
  
  input_GWAS = read.table(paste0(input_path),header=T,na.strings=c("na","NA",""))
  input_GWAS = na.omit(input_GWAS)
  
  
  if (is.null(n_fold)){n_fold=4}
  if (is.null(make_plot)){make_plot=False}
  if (is.null(output_path)&isTRUE(make_plot)){stop(print("Output path for plot is missing!"))}
  if (!is.character(beta_header)){stop(print("Beta header should be a character!"))}
  if (!is.character(se_header)){stop(print("Se header should be a character!"))}
  if (!is.character(pvalue_header)){stop(print("P-value header should be a character!"))}
  if (!is.character(af_header)){stop(print("Allele frequency header should be a character!"))}
  #TODO to create a MAF
  if (is.null(af_header)){print("Allele frequency header is missing. Please provide the correct header for allele frequency.")}
  if (is.null(beta_header)){print("Beta header is missing. Please provide the correct header for effect sizes.")}
  if (is.null(se_header)){print("Se header is missing. Please provide the correct header for standard error.")}
  if (is.null(pvalue_header)){print("P-value header is missing. Please provide the correct header for p-value.")}
  
  
  if(is.numeric(samplesize_header)){
    N.sample = rep(samplesize_header,nrow(input_GWAS))
  }else{
    N.sample = input_GWAS[,samplesize_header]
  }
  if(isTRUE(odds_ratio)){
    beta=log(input_GWAS[,beta_header])
  }else{
    beta=input_GWAS[,beta_header]
  }
  af=input_GWAS[,af_header]
  se=input_GWAS[,se_header]
  pvalue=input_GWAS[,pvalue_header]
  pvalue=pvalue[order(pvalue,decreasing=F)]
  x_var=2*af*(1-af)
  remove(input_GWAS)
  
  
  Res_TildeRL <- FunII.TildeRL(FunII.beta=beta,FunII.SE=se,FunII.sigma=x_var,FunII.N.samplesize=N.sample,FunII.Nv=(1/n_fold)*min(N.sample),FunII.rep=n_fold)
  
  if(make_plot){FunI.plot(Res_TildeRL=Res_TildeRL,p.value=pvalue,input_path=input_path,output_path = output_path)}
  
  print(list(maximal.R2=max(Res_TildeRL),Optimal.P.value=pvalue[which(Res_TildeRL==max(Res_TildeRL))]))
}


#---plot function: --->
FunI.plot<-function(Res_TildeRL,p.value,input_path,output_path){
  png(paste0(output_path),units="in",width = 9,height = 5,res = 300)
  plot(y=Res_TildeRL,x=log(p.value,base=10),main=paste0(input_path),ylab="R^2",xlab="log(P-value)")
  points(y=max(Res_TildeRL),x=log(p.value,base=10)[which(Res_TildeRL==max(Res_TildeRL))],col="red", pch=19)
  abline(v=log(p.value,base=10)[which(Res_TildeRL==max(Res_TildeRL))],col="red",lty=2)
  text(labels=paste0(p.value[which(Res_TildeRL==max(Res_TildeRL))]),y=0.5*max(Res_TildeRL),x=log(p.value,base=10)[which(Res_TildeRL==max(Res_TildeRL))],pos=3,col="red")
  # legend("p-value is : ", col="red", cex=0.8)
  invisible(dev.off())
}


### ---  Funtion File ---###


#--- Prediction Error: --->
FunI.PEtilde<-function(FunI.XtY.tr, FunI.XtY.v, FunI.ntr, FunI.nv, FunI.Sigma,FunI.Var.Y){
  if (! is.vector(FunI.XtY.tr)){stop(print("X_t_Y.tr should be a vector!"))}
  if (! is.vector(FunI.XtY.v)){stop(print("X_t_Y.v should be a vector!"))}
  if (length(FunI.XtY.tr)!=length(FunI.XtY.v)){stop(print("Dimensions of training and testing summary statistics differ!"))}
  
  FunI.beta_tilde = (FunI.XtY.tr/FunI.ntr)/(FunI.Sigma)
  FunI.Var.Y.Hat = cumsum(FunI.Sigma*(FunI.beta_tilde^2))
  FunI.Cov.Y_Y.Hat = cumsum(FunI.beta_tilde*(FunI.XtY.v/FunI.nv))
  FunI.Cor.squared = FunI.Cov.Y_Y.Hat^2/(FunI.Var.Y.Hat*FunI.Var.Y)
  return(FunI.Cor.squared)
}


#--- Resample Summary Statistics: --->
FunII.TildeRL<-function(FunII.beta, FunII.SE, FunII.sigma, FunII.N.samplesize, FunII.Nv, FunII.rep){
  
  if(length(FunII.N.samplesize)!=length(FunII.beta)){
    stop("Error: Inputs lengths do not match!")
  }
  FunII.PErr.Avg <- rep(0, length(FunII.beta))
  
  
  FunII.Ntr = FunII.N.samplesize-FunII.Nv
  FunII.XtY = FunII.beta*FunII.N.samplesize*FunII.sigma
  FunII.Var.Y = quantile(FunII.SE^2*FunII.N.samplesize*FunII.sigma,probs=seq(0,1,0.1))[10]
  #The default for estimation of variance of Y is 90% quantile
  FunII.VarXY = FunII.N.samplesize*(FunII.SE^2)*(FunII.sigma^2)
  
  
  for (FunII.r in 1:FunII.rep){
    FunII.XtY.train.temp = rnorm(n=length(FunII.beta), mean = (FunII.Ntr/FunII.N.samplesize)*FunII.XtY, sd = sqrt((FunII.Ntr*FunII.Nv/FunII.N.samplesize)*FunII.VarXY))
    FunII.XtY.test.temp = FunII.XtY - FunII.XtY.train.temp
    FunII.Zscore.temp = FunII.XtY.train.temp/(sqrt(FunII.N.samplesize*FunII.Ntr)*FunII.SE*FunII.sigma)
    FunII.order.temp = order(abs(FunII.Zscore.temp),decreasing=T)
    
    FunII.PErr.Avg = FunII.PErr.Avg + FunI.PEtilde(FunI.XtY.tr=FunII.XtY.train.temp[FunII.order.temp],
                                                   FunI.XtY.v = FunII.XtY.test.temp[FunII.order.temp],
                                                   FunI.ntr = FunII.Ntr[FunII.order.temp], FunI.nv=FunII.Nv,
                                                   FunI.Sigma = FunII.sigma[FunII.order.temp],FunI.Var.Y=FunII.Var.Y)/FunII.rep
  }
  return(FunII.PErr.Avg)
}

