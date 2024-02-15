#!/s/bin/R35

# PUMAS R2
source("./subsampling/helpers.R")

# Read the argument into R
options(stringsAsFactors=F)
option_list = list(
  make_option("--k", action = "store", default = NA, type = "numeric"),
  make_option("--partitions", action = "store", default = NA, type = "character"),
  make_option("--trait_name", action = "store", default = NA, type = "character"),
  make_option("--gwas_path", action = "store", default = NA, type = "character"),
  make_option("--ld_path", action = "store", default = NA, type = "character"),
  make_option("--output_path", action = "store", default = NA, type = "character"),
  make_option("--chr", action = "store", default = NULL, type = "numeric"),
  make_option("--parallel", action = "store_true", default = FALSE),
  make_option("--threads", action = "store", default = detectCores(), type = "numeric")
)

opt = parse_args(OptionParser(option_list=option_list))
k <- opt$k
partitions <- as.numeric(unlist(lapply((strsplit(opt$partitions, ',')),trimws)))
trait_name <- opt$trait_name
gwas_path <- opt$gwas_path
ld_path <- opt$ld_path
output_path <- opt$output_path
chr <- opt$chr
parallel <- opt$parallel

if (!parallel) {
  threads <- 1
} else {
  if (!opt$threads) {
    threads <- k
  } else {
    threads <- opt$threads
  }
}

# main function
main <- function(){
    if (file.exists(paste0(gwas_path,trait_name,".txt"))) {
      gwas <- fread(paste0(gwas_path,trait_name,".txt"),header=T)
    } else if (file.exists(paste0(gwas_path,trait_name,".gz"))) {
      gwas <- fread(paste0(gwas_path,trait_name,".gz"),header=T)
    } else {
      gwas <- fread(paste0(gwas_path,trait_name,".txt.gz"),header=T)
    }
    if (!is.null(chr)) {
      gwas <- gwas[gwas$`CHR` == chr]
      chr <- paste0(".", chr)
    }

    load(paste0(ld_path, "/ld_ukb.RData"))
    load(paste0(ld_path, "/rs_ukb.RData"))

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
