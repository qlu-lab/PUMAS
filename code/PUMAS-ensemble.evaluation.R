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
source("./evaluation/PUMAS-ensemble.R")
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

### start time for PUMAS-ensemble evaluation
cat("\n\nRunning PUMAS-ensemble evaluation...\n")
begin.time.total = start_time("PUMAS-ensemble evaluation")

### Read reference genotype
cat("\n\nReading reference genotypes...\n")
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
# check NA
tmp.to.exclude = c()
for (ite in 1:k) {
    tmp.to.exclude = c(tmp.to.exclude,as.numeric(which(is.na(colSums(snp.weights.mat[[ite]])==T))))
}
end_time("Sorting PRS weights", begin.time)

if (!is.na(full_weight_path)) {
    cat("\n\nSorting PRS weights for full sumstats...\n")
    begin.time <- start_time("Sorting PRS weights for full sumstats")
    full.snp.weight <- sort_sumstats(prs_methods=prs_method,xty.snp=xty.tmp,iterations=1,full=TRUE)[[1]]

    tmp.to.exclude <- cbind(tmp.to.exclude, as.numeric(which(is.na(colSums(full.snp.weight)==T))))
    # remove NA
    if (length(tmp.to.exclude != 0)) {
        full.snp.weight <- full.snp.weight[,-tmp.to.exclude]
    }
    end_time("Sorting PRS weights for full sumstats", begin.time)
} else {
    full.snp.weight <- NULL
}

# remove NA from snp.weights.mat
if (length(tmp.to.exclude != 0)) {
    for (ite in 1:k) {
        snp.weights.mat[[ite]] = snp.weights.mat[[ite]][,-tmp.to.exclude]
    }
}

### Calculate single PRS sd and exclude some PRS methods
cat("\n\nCalculate single PRS and remove PRS methods...\n")
begin.time = start_time("Calculate single PRS and remove PRS methods")
prs_exclude_list <- mclapply(1:k, function(ite) {
    tmp <- apply(snp.weights.mat[[ite]],2,sd)
    threshold <- median(tmp, na.rm = TRUE) + 3 * IQR(tmp, na.rm = TRUE)
    return(which(tmp >= threshold | tmp == 0 | is.na(tmp)))
}, mc.cores=cores)
prs_exclude <- unique(unlist(prs_exclude_list))
if (length(prs_exclude) != 0){
    for (ite in 1:k) {
        snp.weights.mat[[ite]] <- snp.weights.mat[[ite]][,-prs_exclude]
    }
}
end_time("Calculate single PRS and remove PRS methods", begin.time)

### exclude PRS methods for Full Sumstats
if (!is.na(full_weight_path)) {
    full.snp.weight.all <- full.snp.weight
    if (length(prs_exclude) != 0){
        full.snp.weight <- full.snp.weight[,-prs_exclude]
    }
    full.snp.weight <- apply(full.snp.weight,2,scale,center=F)
}

### run ensemble learning for chosen ensemble type (EN, SL, all)
if (ensemble == "EN") {
    cat("\n\nEN ensemble learning to calculate weights...\n")
    begin.time = start_time("EN ensemble learning")
    en_final_weights = EN_en_learning()
    if (!is.na(full_weight_path)) {
        full.snp.out = cbind(full.snp.weight.all, as.numeric(en_final_weights[,"Ensemble_Weight"]))
        full.snp.out = apply(full.snp.out,2,as.numeric)
        colnames(full.snp.out)[ncol(full.snp.out)] = "EN"
    }
    end_time("EN ensemble learning", begin.time)

} else if (ensemble == "SL") {
    cat("\n\nSL ensemble learning to calculate weights...\n")
    begin.time = start_time("SL ensemble learning")
    sl_final_weights = EN_super_learning()
    if (!is.na(full_weight_path)) {
        full.snp.out = cbind(full.snp.weight.all, as.numeric(sl_final_weights[,"Ensemble_Weight"]))
        full.snp.out = apply(full.snp.out,2,as.numeric)
        colnames(full.snp.out)[ncol(full.snp.out)] = "SL"
    }
    end_time("SL ensemble learning", begin.time)

} else if (ensemble == "all") {
    # run EN first
    cat("\n\nEN ensemble learning to calculate weights...\n")
    begin.time = start_time("EN ensemble learning")
    en_final_weights = EN_en_learning()
    if (!is.na(full_weight_path)) {
        full.snp.out = cbind(full.snp.weight.all, as.numeric(en_final_weights[,"Ensemble_Weight"]))
        colnames(full.snp.out)[ncol(full.snp.out)] = "EN"
        end_time("EN ensemble learning", begin.time)
    }

    # run SL after EN has finished
    cat("\n\nSL ensemble learning to calculate weights...\n")
    begin.time = start_time("SL ensemble learning")
    sl_final_weights = EN_super_learning()
    if (!is.na(full_weight_path)) {
        full.snp.out = cbind(full.snp.out, as.numeric(sl_final_weights[,"Ensemble_Weight"]))
        colnames(full.snp.out)[ncol(full.snp.out)] = "SL"
        full.snp.out = apply(full.snp.out,2,as.numeric)
    }
    end_time("SL ensemble learning", begin.time)

} else {
    cat("\nIncorrect ensemble type (", ensemble, ") given. Please input either EN, SL, or all.")
    stop()
}

# write results to file
if (!is.na(full_weight_path)) {
    fwrite(as.data.frame(cbind(xty.tmp[,c("SNP", "A1")],full.snp.out)), paste0(output_path,trait_name,".ensemble.weights.txt"),col.names = T, row.names = F,sep = "\t",quote = F)
}

# report computing time for PUMAS-ensemble evaluation
end_time("PUMAS-ensemble evaluation", begin.time.total)
