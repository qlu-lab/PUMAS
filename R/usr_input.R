source("PUMAS-master/R/pumas.R")
set.seed(0)


pumas.main<-function(input_path, beta_header,maf_header,se_header, sample_size,flag_plot){
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












