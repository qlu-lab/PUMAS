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
if(!require(parallel)){
    install.packages("parallel")
    library(parallel)
}
if(!require(Rcpp)){
    install.packages("Rcpp")
    library(Rcpp)
}

# import ensemble learning files
source("./code/evaluation/Linear_sumstats.R")
source("./code/evaluation/EN_sumstats.R")

# Read the argument into R
options(stringsAsFactors=F)
option_list = list(
  make_option("--k", action = "store", default = NA, type = "numeric"),
  make_option("--regression", action = "store", default = "EN", type = "character"),
  make_option("--ref_path", action = "store", default = NA, type = "character"),
  make_option("--trait_name", action = "store", default = NA, type = "character"),
  make_option("--prs_method", action = "store", default = NA, type = "character"),
  make_option("--xty_path", action = "store", default = NA, type = "character"),
  make_option("--stats_path", action = "store", default = NA, type = "character"),
  make_option("--weight_path", action = "store", default = NA, type = "character"),
  make_option("--output_path", action = "store", default = NA, type = "character")
)
# assign variables for user provided arguments
opt = parse_args(OptionParser(option_list=option_list))
k <- opt$k
regression <- opt$regression
ref_path <- opt$ref_path
trait_name <- opt$trait_name
prs_method <- as.character(unlist(lapply((strsplit(opt$prs_method, ',')),trimws)))
xty_path <- opt$xty_path
stats_path <- opt$stats_path
weight_path <- opt$weight_path
output_path <- opt$output_path

# start time for ensemble learning
begin.time.total = Sys.time()
cat("\nEnsemble learning ... begins at ", format(begin.time.total, format = "%F %R %Z"), ":", sep = "", "\n\n")

# run ensemble learning for chosen regression type
if (regression == "EN") {
    EN_ensemble_learning()
} else if (regression == "linear") {
    linear_ensemble_learning()
} else {
    cat("\nIncorrect regression type (", regression, ") given. Please input either EN or linear.")
    stop()
}

# report computing time for ensemble learnings
end.time.total <- Sys.time()
total.time <- difftime(time1=end.time.total,time2=begin.time.total,units="sec")
hrs <- floor(floor(total.time)/3600)
mins <- floor(floor(floor(total.time) - hrs * 3600)/60)
secs <- total.time - hrs*3600 - mins*60
cat("\nEnsemble learning ended at ", format(end.time.total, format = "%F %R %Z"), ". \n  Total time: ", hrs, " hours, ", mins, " minutes and ", secs, " seconds.\n\n", sep = "")