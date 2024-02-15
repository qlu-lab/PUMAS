#!/s/bin/R35

if(!require(data.table)){
    install.packages("data.table")
    suppressMessages(library(data.table))
}
if(!require(BEDMatrix)){
    install.packages("BEDMatrix")
    suppressMessages(library(BEDMatrix))
}
if(!require(optparse)){
    install.packages("optparse")
    suppressMessages(library(optparse))
}
if(!require(parallel)){
    install.packages("parallel")
    suppressMessages(library(parallel))
}
if(!require(Rcpp)){
    install.packages("Rcpp")
    suppressMessages(library(Rcpp))
}

### import ensemble learning files
source("./evaluation/EN_sumstats.R")
source("./evaluation/SL_sumstats.R")
source("./evaluation/helpers.R")
sourceCpp("./evaluation/CoordDescent.cpp")

### Read input arguments into R
options(stringsAsFactors=F)
option_list = list(
  make_option("--k", action = "store", default = NA, type = "numeric"),
  make_option("--ensemble", action = "store", default = "EN", type = "character"),
  make_option("--ref_path", action = "store", default = NA, type = "character"),
  make_option("--trait_name", action = "store", default = NA, type = "character"),
  make_option("--prs_method", action = "store", default = NA, type = "character"),
  make_option("--xty_path", action = "store", default = NA, type = "character"),
  make_option("--stats_path", action = "store", default = NA, type = "character"),
  make_option("--weight_path", action = "store", default = NA, type = "character"),
  make_option("--full_weight_path", action = "store", default = NA, type = "character"),
  make_option("--output_path", action = "store", default = NA, type = "character"),
  make_option("--thr", action = "store", default = 1e-5, type = "numeric"),
  make_option("--parallel", action = "store_true", default = FALSE),
  make_option("--threads", action = "store", default = NULL, type = "numeric")
)
### assign variables for user provided arguments
opt = parse_args(OptionParser(option_list=option_list))
k <- opt$k
ensemble <- opt$ensemble
ref_path <- opt$ref_path
trait_name <- opt$trait_name
prs_method <- as.character(unlist(lapply((strsplit(opt$prs_method, ',')),trimws)))
xty_path <- opt$xty_path
stats_path <- opt$stats_path
weight_path <- opt$weight_path
full_weight_path <- opt$full_weight_path
output_path <- opt$output_path
threshold <- opt$thr
parallel <- opt$parallel
if (!parallel) {
    cores <- 1
} else {
    if (is.null(opt$threads)) {
        cores <- min(k, detectCores() - 1)
    } else {
        cores <- opt$threads
    }
}

### start time for PUMA-CUBS evaluation
cat("\n\nRunning PUMA-CUBS evaluation ...\n")
begin.time.total = start_time("PUMA-CUBS evaluation")

### Read reference genotype
cat("\n\nReading reference genotypes ...\n")
begin.time <- start_time("Reading reference genotypes")
ref.geno <- BEDMatrix(ref_path)
rs.geno <- fread(paste0(ref_path,".bim"),header=F)$V2
xty.tmp <- fread(paste0(xty_path,trait_name,".xty.omnibus.ite1.txt"),header=T)
ref.geno <- as.matrix(ref.geno[,match(xty.tmp$SNP,rs.geno)])
end_time("Reading reference genotypes", begin.time)

### Sort PRS weights
cat("\n\nSorting PRS weights ...\n")
begin.time <- start_time("Sorting PRS weights")
snp.weights.mat <- sort_sumstats(prs_methods=prs_method,xty.snp=xty.tmp,iterations=k)
end_time("Sorting PRS weights", begin.time)

### Calculate single PRS r2 and exclude some PRS methods
cat("\n\nCalculate single PRS and remove PRS methods ...\n")
begin.time = start_time("Calculate single PRS and remove PRS methods")
prs_exclude_list <- mclapply(1:k, function(ite) {
    tmp <- single_ss_main(X.ref=ref.geno, snp.weight=snp.weights.mat[[ite]], ite=ite)
    tmp.ratio <- (tmp[,2])/(tmp[,1])
    return(which(tmp.ratio >= (median(tmp.ratio) + 1.5*IQR(tmp.ratio))))
}, mc.cores=cores)
prs_exclude <- unique(unlist(prs_exclude_list))
write.table(prs_exclude,paste0(output_path,"/prs_toExclude.txt"),col.names = F, row.names = F,sep = "\n",quote = F)
end_time("Calculate single PRS and remove PRS methods", begin.time)

### Sort PRS Weights for Full Sumstats
if (!is.null(full_weight_path)) {
    cat("\n\nSorting PRS weights for full sumstats...\n")
    begin.time <- start_time("Sorting PRS weights for full sumstats")
    full.snp.weight <- sort_sumstats(prs_methods=prs_method,xty.snp=xty.tmp,iterations=1,full=TRUE)
    full.snp.weight <- apply(full.snp.weight[[1]][,-prs_exclude],2,scale,center=F)
    end_time("Sorting PRS weights for full sumstats", begin.time)
}

### run ensemble learning for chosen ensemble type (EN, SL, all)
if (ensemble == "EN") {
    cat("\n\nEN ensemble learning to calculate weights ...\n")
    begin.time = start_time("EN ensemble learning")
    EN_ensemble_learning()
    end_time("EN ensemble learning", begin.time)
} else if (ensemble == "SL") {
    cat("\n\nSuper learning to calculate weights ...\n")
    begin.time = start_time("Super learning")
    EN_super_learning()
    end_time("Super learning", begin.time)
} else if (ensemble == "all") {
    # run EN first
    cat("\n\nEN ensemble learning to calculate weights ...\n")
    begin.time = start_time("EN ensemble learning")
    EN_ensemble_learning()
    end_time("EN ensemble learning", begin.time)
    # run SL after EN has finished
    cat("\n\nSuper learning to calculate weights ...\n")
    begin.time = start_time("Super learning")
    EN_super_learning()
    end_time("Super learning", begin.time)
} else {
    cat("\nIncorrect ensemble type (", ensemble, ") given. Please input either EN, SL, or all.")
    stop()
}

# report computing time for PUMA-CUBS evaluation
end_time("PUMA-CUBS evaluation", begin.time.total)
