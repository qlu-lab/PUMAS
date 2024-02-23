# GWAS preparation for PUMAS/PUMACUBS
options(stringsAsFactors=F)

if(!require(data.table)){
    install.packages("data.table")
    library(data.table)
}
if(!require(optparse)){
    install.packages("optparse")
    library(optparse)
}

# Read the argument into R
options(stringsAsFactors=F)
option_list = list(
  # required
  make_option("--file_path", action = "store", default = NULL, type = "character"),
  make_option("--frq_path", action = "store", default = NULL, type = "character"),
  make_option("--output_path", action = "store", default = NULL, type = "character"),
  make_option("--snp", action = "store", default = NULL, type = "character"),
  make_option("--a1", action = "store", default = NULL, type = "character"),
  make_option("--a2", action = "store", default = NULL, type = "character"),
  make_option("--stat", action = "store", default = NULL, type = "character"),
  make_option("--p", action = "store", default = NULL, type = "character"),
  # optional
  make_option("--chr", action = "store", default = NULL, type = "character"),
  make_option("--bp", action = "store", default = NULL, type = "character"),
  make_option("--OR", action = "store_true", default = FALSE, type = "logical"),
  make_option("--logit", action = "store_true", default = FALSE, type = "logical"),
  make_option("--se", action = "store", default = NULL, type = "character"),
  make_option("--maf", action = "store", default = NULL, type = "character"),
  make_option("--n.total", action = "store", default = NULL, type = "numeric"),
  make_option("--n.case", action = "store", default = NULL, type = "numeric"),
  make_option("--n.con", action = "store", default = NULL, type = "numeric"),
  make_option("--n.col", action = "store", default = NULL, type = "character"),
  make_option("--n.case.col", action = "store", default = NULL, type = "character"),
  make_option("--n.con.col", action = "store", default = NULL, type = "character")
)

opt = parse_args(OptionParser(option_list=option_list))
# required
file_path <- opt$file_path
frq_path <- opt$frq_path
output_path <- opt$output_path
snp.raw <- opt$snp
a1.raw <- opt$a1
a2.raw <- opt$a2
stat.raw <- opt$stat
p.raw <- opt$p

# optional
n.total.raw <- opt$n.total
chr.raw <- opt$chr
bp.raw <- opt$bp
OR.raw <- opt$OR
binary.raw <- opt$logit
se.raw <- opt$se
maf.raw <- opt$maf
n.case.raw <- opt$n.case
n.con.raw <- opt$n.con
n.col.raw <- opt$n.col
n.case.col.raw <- opt$n.case.col
n.con.col.raw <- opt$n.con.col

gwas_qc <- function(gwas,chr,bp,snp,a1,a2,stat,p,binary,se,maf,n.total,n_case,n_con,n_col,n_case_col,n_con_col,freq,OR){
    
    # 0. make an empty dataframe and initiate values
    gwas.tmp <- data.frame(CHR=rep(NA,nrow(gwas)),BP=rep(NA,nrow(gwas)),SNP=rep(NA,nrow(gwas)),A1=rep(NA,nrow(gwas)),A2=rep(NA,nrow(gwas)),MAF=rep(NA,nrow(gwas)),BETA=rep(NA,nrow(gwas)),SE=rep(NA,nrow(gwas)),P=rep(NA,nrow(gwas)),N=rep(NA,nrow(gwas)),N_imp=rep(NA,nrow(gwas)))
    
    gwas.tmp$SNP <- gwas[,snp]
    gwas.tmp$A1 <- toupper(gwas[,a1])
    gwas.tmp$A2 <- toupper(gwas[,a2])
    gwas.tmp$P <- as.numeric(gwas[,p])
    
    
    if (!is.null(chr)){ # chr
        if (length(grep(x=gwas[,chr],pattern="[Cc][Hh][Rr](.*)"))!=0){
            gwas.tmp$CHR <- as.numeric(gsub(x=gwas[,chr],pattern="[Cc][Hh][Rr](.*)",replacement="\\1"))
        } else {
            gwas.tmp$CHR <- as.numeric(gwas[,chr])
        }
    } else {
        gwas.tmp$CHR <- NULL
    }
    
    if (!is.null(bp)){ # bp
        gwas.tmp$BP <- gwas[,bp]
    } else {
        gwas.tmp$BP <- NULL
    }
    
    if (!binary & OR) {
        stop("When logit = F, cannot set OR = T.")
    }

    # if OR=T, then change OR to beta
    if (OR){ # beta
        gwas.tmp$BETA <- log(as.numeric(x=gwas[,stat]))
    } else {
        gwas.tmp$BETA <- as.numeric(gwas[,stat])
    }
    
    if (!is.null(maf)){ # maf
        gwas.tmp$MAF <- as.numeric(gwas[,maf])
    } else {
        gwas.tmp$MAF <- NULL
    }
    
    if (binary){ # sample size
        if ((!is.null(n_case_col)) & (!is.null(n_con_col))){
            gwas.tmp$N <- as.numeric(gwas[,n_case_col]) + as.numeric(gwas[,n_con_col])
            gwas.tmp$ratio <- as.numeric(gwas[,n_case_col])/gwas.tmp$N
        } else if ((!is.null(n_case)) & (!is.null(n_con))){
            gwas.tmp$N <- ifelse(is.null(n_col), n_case + n_con, as.numeric(gwas[,n_col]))
            gwas.tmp$ratio <- n_case/gwas.tmp$N
        } else {
            stop("User inputs don't comply with the sample size requirement: (1) when logit = F, please provide either n.total or n.col; (2) when logit = T, please provide either (n.case and n.con) or (n.case.col and n.con.col; n.col is recommended to be provided too).")
        }
    } else {
        if (!is.null(n_col)){
           gwas.tmp$N <- as.numeric(gwas[,n_col])
       } else if (!is.null(n.total)) {
           gwas.tmp$N <- n.total
       } else {
            stop("User inputs don't comply with the sample size requirement: (1) when logit = F, please provide either n.total or n.col; (2) when logit = T, please provide either (n.case and n.con) or (n.case.col and n.con.col; n.col is recommended to be provided too).")
        }
    }

    
    # 1. adjust SNPs with 0 p-value
    gwas.tmp$P <- ifelse(as.numeric(gwas.tmp$P==0),min(as.numeric(gwas.tmp$P)[as.numeric(gwas.tmp$P)!=0]),as.numeric(gwas.tmp$P))
    
    
    # 2. initiate SE; impute SE if needed for linesr regression summary statistics
    if (!binary){
        if (is.null(se)){
            gwas.tmp$SE <- abs(gwas.tmp$BETA/qnorm(gwas.tmp$P/2,0,1))
        } else {
            gwas.tmp$SE <- as.numeric(gwas[,se])
        }
        gwas.tmp <- gwas.tmp[!is.na(gwas.tmp$SE),]
    }

    
    # 3. take overlap between hm3 and 1kg CEU SNPs
    rs_overlap <- intersect(gwas.tmp$SNP,freq$SNP)
    gwas.tmp <- gwas.tmp[match(rs_overlap,gwas.tmp$SNP),]
    freq <- freq[match(rs_overlap,freq$SNP),]
    

    # 4. coordinate A1/A2 with 1kg, remove/flip SNPs if needed (beta/MAF, A1 and A2)
    snp.flip <- (freq$A1==gwas.tmp$A2) & (freq$A2==gwas.tmp$A1)
    A1.tmp <- gwas.tmp$A1
    A2.tmp <- gwas.tmp$A2
    gwas.tmp$A1[snp.flip] <- A2.tmp[snp.flip]
    gwas.tmp$A2[snp.flip] <- A1.tmp[snp.flip]
    gwas.tmp$BETA[snp.flip] <- -gwas.tmp$BETA[snp.flip]
    if (!is.null(maf)){
        gwas.tmp$MAF[snp.flip] <- 1 - gwas.tmp$MAF[snp.flip]
    }
    gwas.tmp <- gwas.tmp[(freq$A1==gwas.tmp$A1) & (freq$A2==gwas.tmp$A2),]

    
    # 5. impute MAF if needed and QC based on MAF difference
    freq <- freq[match(gwas.tmp$SNP,freq$SNP),]
    if (is.null(maf)){
        gwas.tmp$MAF <- as.numeric(freq$MAF)
    } else {
        gwas.tmp <- gwas.tmp[abs(gwas.tmp$MAF-as.numeric(freq$MAF))<=0.05,]
    }
    
    
    # 6. transform logistic beta and SE to lpm ##
    if (binary){
        Var.G <- 2*gwas.tmp$MAF*(1-gwas.tmp$MAF)
        gwas.tmp$SE <- sqrt(gwas.tmp$ratio*(1-gwas.tmp$ratio)/(gwas.tmp$N*Var.G))
        gwas.tmp$BETA <- qnorm(gwas.tmp$P/2,0,1,lower.tail=F) * sign(gwas.tmp$BETA) * gwas.tmp$SE
        gwas.tmp$ratio <- NULL
    }
    
    
    # 7. QC on extrem numerical values
    gwas.tmp$BETA <- as.numeric(gwas.tmp$BETA)
    gwas.tmp$SE <- as.numeric(gwas.tmp$SE)
    gwas.tmp$P <- as.numeric(gwas.tmp$P)
    gwas.tmp$MAF <- as.numeric(gwas.tmp$MAF)
    gwas.tmp$N <- as.numeric(gwas.tmp$N)
    gwas.tmp <- gwas.tmp[gwas.tmp$SE!=0 & gwas.tmp$SE!=Inf & gwas.tmp$SE!=-Inf & !is.na(gwas.tmp$SE),]
    gwas.tmp <- gwas.tmp[gwas.tmp$MAF!=0 & gwas.tmp$MAF!=1 & !is.na(gwas.tmp$MAF),]
    gwas.tmp <- gwas.tmp[gwas.tmp$N!=0 & !is.na(gwas.tmp$N),]
    gwas.tmp <- gwas.tmp[gwas.tmp$BETA!=Inf & gwas.tmp$BETA!=-Inf & !is.na(gwas.tmp$BETA),]
    
    
    # 8. impute N_eff based on https://www.biorxiv.org/content/10.1101/2021.03.29.437510v4 and QC on imputed sample size
    Var.G <- 2*gwas.tmp$MAF*(1-gwas.tmp$MAF)
    Var.Y <- quantile((gwas.tmp$SE^2)*gwas.tmp$N*Var.G,probs=seq(0,1,0.1))[10]
    gwas.tmp$N_imp <- (Var.Y/Var.G - gwas.tmp$BETA^2)/(gwas.tmp$SE^2)
    
    # real N
    N_thresh_small <- (quantile(gwas.tmp$N,probs=seq(0,1,0.1))[10])*0.7
    N_thresh_big <- (quantile(gwas.tmp$N,probs=seq(0,1,0.1))[10])*0.9
    N_thresh_lower10 <- (quantile(gwas.tmp$N,probs=seq(0,1,0.1))[2])
    
    # imputed N
    N_imp_thresh_lower10 <- (quantile(gwas.tmp$N_imp,probs=seq(0,1,0.1))[2])

    if ((!binary) & (!is.null(n_col))) { # continuous traits with N column
        if (n.total < 1e5 | N_thresh_lower10 >= N_thresh_big) {
            gwas.tmp <- gwas.tmp[gwas.tmp$N>=N_thresh_big,]
        } else {
            gwas.tmp <- gwas.tmp[gwas.tmp$N>=N_thresh_small,]
        }
        gwas.tmp$N_imp <- NULL

    } else if (binary & (!is.null(n_col)|!is.null(n_case_col))){ # binary traits with N column
        if (n.total < 1e5 | N_thresh_lower10 >= N_thresh_big) {
            gwas.tmp <- gwas.tmp[gwas.tmp$N>=N_thresh_big,]
        } else {
            gwas.tmp <- gwas.tmp[gwas.tmp$N>=N_thresh_small,]
        }
        gwas.tmp$N_imp <- NULL
        
    } else {
        if (n.total < 1e5 | N_imp_thresh_lower10 >= n.total*0.9) {
            gwas.tmp <- gwas.tmp[gwas.tmp$N_imp>=n.total*0.9 & gwas.tmp$N_imp<=n.total*1.1,]
        } else {
            gwas.tmp <- gwas.tmp[gwas.tmp$N_imp>=n.total*0.7 & gwas.tmp$N_imp<=n.total*1.1,]
        }
        gwas.tmp$N <- round(gwas.tmp$N_imp,0)
        gwas.tmp$N_imp <- NULL
    }
    
    
    ## sort gwas by chr and bp
    if (!is.null(chr) & !is.null(bp)){
        gwas.tmp$BP <- as.numeric(gwas.tmp$BP)
        gwas.tmp <- gwas.tmp[with(gwas.tmp,order(CHR,BP,decreasing=F)),]
    }
    
    
    ## return QCed sumstats
    if (nrow(gwas.tmp) < 2e5){
        warning("Number of SNPs remained is smaller than 200,000 after QC, suggesting poor GWAS quality.")
    }
    
    return(gwas.tmp)
}


# executable
frq <- fread(frq_path, header=T) # change to user input folder path

gwas.raw <- as.data.frame(fread(file_path, header=T)) # need user to input $file_path

gwas.QCed <- gwas_qc(gwas=gwas.raw,chr=chr.raw,bp=bp.raw,snp=snp.raw,a1=a1.raw,a2=a2.raw,maf=maf.raw,stat=stat.raw,p=p.raw,binary=binary.raw,se=se.raw,n.total=n.total.raw,n_case=n.case.raw,n_con=n.con.raw,n_col=n.col.raw,n_case_col=n.case.col.raw,n_con_col=n.con.col.raw,freq=frq,OR=OR.raw)

# retrieve trait_name from $file_path
trait_name <- unlist(strsplit(tail(unlist(strsplit(file_path, "/")), 1), "[.]"))[1]

fwrite(gwas.QCed,paste0(output_path,"/",trait_name,".txt.gz"),col.names=T,row.names=F,sep="\t",quote=F,na=NA)
