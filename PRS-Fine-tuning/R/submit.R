source("/Users/yupeilin/PRS-Fine-tuning/R/prs_tuning.R")
set.seed(0)


##################################################
data.real <- read.table(paste0("/Users/yupeilin/PRS-Fine-tuning/input/T0030_pruned.txt"), header=T,na.strings=c("na","NA",""))
data.real=na.omit(data.real)

beta=data.real$Beta
MAF=data.real$EAF
se=data.real$SE
N.sample=766345
X_VAR=2*MAF*(1-MAF)


#p-value selection
p.value=data.real$Pval
p.value=p.value[order(p.value,decreasing=F)]
index.cons = which(p.value<0.0001)

Res_TildeRL <- FunII.TildeRL(FunII.beta=beta,FunII.SE=se,FunII.sigma=X_VAR,FunII.N.samplesize=N.sample,FunII.Nv=0.25*N.sample,FunII.rep=4, FunII.Corr=F)
#p-value selection
Res_TildeRL=Res_TildeRL[-index.cons]
FunI.plot(Res_TildeRL=Res_TildeRL,p.value=p.value)
