# PUMAS/PUMA-CUBS
**PUMAS** and **PUMA-CUBS** are summary-statistcis-based method to fine-tune, combine, and benchmark PRS methods using only GWAS summary statistics and a LD reference panel. If the PRS fine-tuning is the only task, please use **PUMAS** functions. Otherwise to achieve all three objectives, please use **PUMA-CUBS**. A workflow of PUMAS/PUMA-CUBS is shown below ![here](https://github.com/qlu-lab/PUMAS/blob/master/Workflow.png)

## Announcements
* We are currently preparing additional LD reference datasets.
* Previous version of PUMAS for fine-tuning P+T/C+T PRSs is available [here](https://github.com/qlu-lab/PUMAS/tree/original).

## Version History
* 11/04/2022: Uploaded a tutorial for PUMAS and PUMA-CUBS.

## Getting Started
* Clone this repository by `git clone https://github.com/qlu-lab/PUMAS.git`
* Downlaod the LD reference data constructed using 1000 Genomes Project Phase III European ancestry data
  * **Approximately independent LD blocks**:
    * `wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/PUMAS/LD/ld_1kg.RData`
    * `wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/PUMAS/LD/rs_1kg.RData`
  * **Genotype data from the LD panel**:
    * `wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/PUMAS/LD/1kg_hm3_QCed_noM*`

## Using PUMAS
### Subsample training and tuning summary statistics
For PUMAS/PUMA-CUBS to subsample GWAS summary statistics from a full GWAS summary-level data, two datasets are requried:
  * **Approximately independent LD blocks**
  * **GWAS summary statistics**. A tab delimited file in `.txt` or `.gz` format containing at least the following fields:
    * SNP: SNP identifier (rsID)
    * A1: effective allele
    * A2: other allele
    * MAF: minor allele frequency
    * BETA: SNP effect size estimation
    * SE: standard error for BETA
    * P: P-value
    * N: sample size

    * An example of formatted GWAS summary statistics will look like this:
  ```
    CHR     BP      SNP             A1      A2      MAF     BETA    SE      P       N
    1       779322  rs4040617       G       A       0.1353  -0.001  0.0023  0.65    592497
    1       785989  rs2980300       T       C       0.1473  -0.0011 0.0023  0.61    591333
    1       1003629 rs4075116       C       T       0.283   0.0034  0.0016  0.037   696881
    ...
  ```
After formatting GWAS summary statistics and downloading LD blocks, run:
```
Rscript ./code/PUMAS.subsampling.R \
--k <number of folds> \
--partitions <training>,<tuning> \
--trait_name <trait name> \
--gwas_path <GWAS sumstats folder> \
--ld_path <ld folder> \
--output_path <output folder>
```
  * `k`: number of folds in PUMAS's implementation of Monte Carlo cross-validation (e.g., `--k 4`)
  * `partitions`: subsets' sample size proportion compared to total samples (e.g., `--partitions 0.75,0.25`)
  * `trait_name`: file name of GWAS summary statistics
  * `gwas_path`: folder containing GWAS summary statistics
  * `ld_path`: path to approximately independent LD blocks
  * `output_path`: folder to write partitioned GWAS summary statistics

### Evaluate PRS performance
After partitioning summary statistics, users can train any PRS method using the subsampled training summary statistics. Then, to use PUMAS for evaluating and fine-tuning PRS methods, three datasets are required:
  * **Tuning summary statistics** obtained from the subsampling step
  * **Genotype data from the LD panel**
  * **SNP weights**. A flat `.<pre_method>.txt` file containing SNP weights for each tuning parameter in a PRS framework (there can be multiple columns of SNP weights). An example of SNP weights file look like this:
```
  CHR     SNP             A1      A2      s.0.2_lambda.5e.3       s.0.2_lambda.1e.2       s.0.2_lambda.5e.2       s.0.5_lambda.5e.3       s.0.5_lambda.1e.2       s.0.5_lambda.5e.2       s.0.9_lambda.5e.3       s.0.9_lambda.1e.2      s.0.9_lambda.5e.2
  1       rs4040617       G       A       0                       0                       0                       0                       0                 0                       0                       0                      0
  1       rs2980300       T       C       0                       0                       0                       0                       0                 0                       0                       0                      0
  1       rs4075116       C       T       0                       0                       0                       0                       0                 0                       0                       0                      0
  ...
```

After gathering all necessary datasets, run:
```
Rscript ./code/PUMAS.evaluation.R \
--k <number of folds> \
--ref_path <LD ref> \
--trait_name <trait name> \
--prs_method <prs_method> \
--xty_path <subsampled sumstats folder> \
--stats_path <statistics folder> \
--weight_path <SNP weights> \
--output_path <output folder>
```
  * `k`: number of folds in PUMAS's implementation of Monte Carlo cross-validation
  * `ref_path`: path to the LD genotype data
  * `trait_name`: file name for both subsampled summary statistics and SNP weight txt file (e.g., simply put `Height` if the subsampled summary statistics file is named `Height.gwas.ite1.txt` and the SNP weights file is named `Height.lassosum.txt`)
  * `prs_method`: PRS method name (e.g., put `lassosum` if the SNP weights file is named `Height.lassosum.txt`)
  * `xty_path`: folder containing partitioned summary statistics from subsampling step
  * `stats_path`: folder containing statistics from subsampling step (e.g, variance of phenotype and sample size for each sumstats)
  * `weight_path`: folder containing SNP weights file
  * `output_path`: folder to write PRS evalation and model-tuning results by PUMAS
  
## Using PUMA-CUBS

### Subsample training, tuning, ensemble training, and testing summary statistics
PUMA-CUBS uses exactly the same inputs as PUMAS. The only difference between implementation between PUMAS and PUMA-CUBS is scripting. To partition full GWAS summary statistics to four different subsets, run:
```
Rscript ./code/PUMA-CUBS.subsampling.R \
--k <number of folds> \
--partitions <training>,<tuning>,<ensemble training>,<testing> \
--trait_name <trait name> \
--gwas_path <GWAS sumstats folder> \
--ld_path <ld folder> \
--output_path <output folder>
```
  * `k`: number of folds in PUMA-CUBS's implementation of Monte Carlo cross-validation (e.g., `--k 4`)
  * `partitions`: subsets' sample size proportion compared to total samples (e.g., `--partitions 0.6,0.2,0.1,0.1`)
  * `trait_name`: file name of GWAS summary statistics
  * `gwas_path`: folder containing GWAS summary statistics
  * `ld_path`: path to approximately independent LD blocks
  * `output_path`: folder to write partitioned GWAS summary statistics
  
### Construct ensemble PRS and benchmark PRS models
Different from PUMAS's PRS evaluation function, PUMA-CUBS requires SNP weights from each PRS method to be stored in a separate `.<pre_method>.txt` file so that PUMA-CUBS can construct ensemble PRS based on fine-tuned PRS model from each method and benchmark all PRS models. After PRS model training, run:
```
Rscript ./code/PUMA-CUBS.evaluation.R \
--k <number of folds> \
--ref_path <LD ref> \
--trait_name <trait name> \
--prs_method <prs_methods> \
--xty_path <subsampled sumstats folder> \
--stats_path <statistics folder> \
--weight_path <SNP weights> \
--output_path <output folder>
```
  * `k`: number of folds in PUMA-CUBS's implementation of Monte Carlo cross-validation
  * `ref_path`: path to the LD genotype data
  * `trait_name`: file name for all subsampled summary statistics and SNP weight txt file
  * `prs_method`: PRS methods' names (e.g., put `lassosum,prscs,ldpred2` if the SNP weights files are named `Height.lassosum.txt`,`Height.prscs.txt`, and `Height.ldpred2.txt`)
  * `xty_path`: folder containing partitioned summary statistics from subsampling step
  * `stats_path`: folder containing statistics from subsampling step (e.g, variance of phenotype and sample size for each sumstats)
  * `weight_path`: folder containing SNP weights files
  * `output_path`: folder to write PRS model-tuning and benchmarking results by PUMA-CUBS
  
## Output

## Citation
* If you use PUMAS/PUMA-CUBS, please cite:
  Zhao, Z., Gruenloh, T., Wu, Y., Sun, Z., Miao, J., Wu, Y., Song, J., & Lu, Q. (2022). [Optimizing and benchmarking polygenic risk scores with GWAS summary statistics](https://www.biorxiv.org/content/10.1101/2022.10.26.513833v1). *bioRxiv*.
* If you use PUMAS for fine-tuning P+T/C+T PRS, pleasse cite:
  Zhao, Z., Yi, Y., Song, J., Wu, Y., Zhong, X., Hohman, T.J., Fletcher, J., & Lu, Q. (2021). [PUMAS: fine-tuning polygenic risk scores with GWAS summary statistics](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02479-9). *Genome Biology*, 22,257.
  
## Support
Please send questions and issues related to PUMAS/PUMA-CUBS software to Zijie Zhao (zzhao232@wisc.edu) and Qiongshi Lu (qlu@biostat.wisc.edu).
