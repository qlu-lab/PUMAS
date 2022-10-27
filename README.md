# PUMAS/PUMA-CUBS
 * A tutorial for PUMA-CUBS is coming soon
 * Check out our new preprint "Optimizing and benchmarking polygenic risk scores with GWAS summary statistics" [here](https://www.biorxiv.org/content/10.1101/2022.10.26.513833v1).
 * Previous version of PUMAS for fine-tuning P+T/C+T PRSs is available [here](https://github.com/qlu-lab/PUMAS/tree/original).

## Input
### Subsampling step
The input is GWAS summary statistics. This is a flat file containing at least the following fields:

| Name | Description  |
|----------------|------------------------------------------------------------------------------|
| SNP | SNP identifier (rsID) |
| A1  | Effective allele  |
| A2  | Other allele  |
| MAF | Minor allele frequency  |
| BETA |  Effect size |
| SE  | Standard error  |    
| P | P-value |   
| N | Sample size |

```
  CHR     BP      SNP             A1      A2      MAF     BETA    SE      P       N       N_case  N_con   N_imp                   N_imp_bin
  1       779322  rs4040617       G       A       0.1353  -0.001  0.0023  0.65    592497  NA      NA      633045.158654276        3231551.21093348
  1       785989  rs2980300       T       C       0.1473  -0.0011 0.0023  0.61    591333  NA      NA      589656.23993338 3010060.88896371
  1       1003629 rs4075116       C       T       0.283   0.0034  0.0016  0.037   696881  NA      NA      754233.108093159        3850205.67504973
  ...
```

### Evaluation/Ensemble step
A flat file with a header row containing the following fields (there can be multiple columns of SNP weights):

| Name | Description  |
|----------------|------------------------------------------------------------------------------|
| SNP | SNP identifier (rsID) |
| A1  | Effective allele  |
| A2  | Other allele  |
| *weight* |  SNP weights for PRS|

```
  CHR     SNP             A1      A2      s.0.2_lambda.5e.3       s.0.2_lambda.1e.2       s.0.2_lambda.5e.2       s.0.5_lambda.5e.3       s.0.5_lambda.1e.2       s.0.5_lambda.5e.2       s.0.9_lambda.5e.3       s.0.9_lambda.1e.2      s.0.9_lambda.5e.2
  1       rs4040617       G       A       0                       0                       0                       0                       0                 0                       0                       0                      0
  1       rs2980300       T       C       0                       0                       0                       0                       0                 0                       0                       0                      0
  1       rs4075116       C       T       0                       0                       0                       0                       0                 0                       0                       0                      0
  ...
```

## PUMAS
* Subsampling
```
Rscript ./code/PUMAS.subsampling.R \
--k 4 \
--partitions 0.75,0.25 \
--trait_name <trait name> \
--gwas_path <data folder> \
--ld_path <ld folder> \
--output_path <output folder>
```
* Evaluation
```
Rscript ./code/PUMAS.evaluation.R \
--k 4 \
--ref_path <LD ref> \
--trait_name <trait name> \
--prs_method lassosum \
--xty_path <subsampling output folder> \
--stats_path <subsampling output folder> \
--weight_path <self-calculated snp weights> \
--output_path <output folder>
```

## PUMA-CUBS
* Subsampling
```
Rscript ./code/PUMA-CUBS.subsampling.R \
--k 4 \
--partitions 0.6,0.2,0.1,0.1 \
--trait_name <trait name> \
--gwas_path <data folder> \
--ld_path <ld folder> \
--output_path <output folder>
```
* Evaluation
```
Rscript ./code/PUMA-CUBS.evaluation.R \
--k 4 \
--ref_path <LD ref> \
--trait_name <trait name> \
--prs_method lassosum \
--xty_path <subsampling output folder> \
--stats_path <subsampling output folder> \
--weight_path <self-calculated snp weights> \
--output_path <output folder>
```
