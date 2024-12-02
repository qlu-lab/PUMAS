# PUMAS/PUMAS-ensemble
**PUMAS** and **PUMAS-ensemble** are summary-statistcis-based method to fine-tune, combine, and benchmark PRS methods using only GWAS summary statistics and a LD reference panel. If the PRS fine-tuning is the only task, please use **PUMAS** functions. Otherwise to achieve all three objectives, please use **PUMAS-ensemble**. A workflow of PUMAS/PUMAS-ensemble is shown below   
![here](https://github.com/qlu-lab/PUMAS/blob/master/Workflow.png)   

## Announcements
* We are currently preparing additional LD reference datasets.
* Previous version of PUMAS with simple ensemble learning is available [here](https://github.com/qlu-lab/PUMAS/tree/PUMASv2).
* Previous version of PUMAS for fine-tuning P+T/C+T PRSs is available [here](https://github.com/qlu-lab/PUMAS/tree/original).

## Version History
* 11/04/2022: Upload a tutorial for PUMAS and PUMAS-ensemble.
* 01/30/2023: Upload a script and tutorial for cleaning GWAS sumamry statistics.
* 12/02/2024: Update PUMAS with PUMAS-EN and PUMAS-SL.

## Getting Started
* Clone this repository by `git clone https://github.com/qlu-lab/PUMAS.git`
* Downlaod the LD reference data constructed using 1000 Genomes Project Phase III European ancestry data
  * **Approximately independent LD blocks**:
    * `wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/PUMAS/LD/ld_1kg.RData`
    * `wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/PUMAS/LD/rs_1kg.RData`
  * **Genotype data from the LD panel**:
    * `wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/PUMAS/LD/1kg_hm3_QCed_noM*`
  * **Frequency data from the LD panel**:
    * `wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/PUMAS/Freq/1kg_hm3_QCed_noM_freq.frq`
  * If `wget` doesn't work, download the data above via [box folder](https://uwmadison.box.com/s/6yv7u8wxm6zutj7763jekdhed47kl0f1).
* Install the following R (>=3.5.1) dependencies by `install.packages()`:
  * optparse
  * data.table
  * BEDMatrix
  * parallel

## GWAS summary statistics preparation
We highly recommend that users clean their summary statistics prior to applying PUMAS/PUMACUBS. Here we provide a GWAS sumstats QC script. **Please make sure that the input GWAS sumstats has rsID for each SNP.** To use the GWAS QC script, run:
```
Rscript ./code/gwas_qc.R \
--file_path <raw GWAS sumstats path> \ # required
--frq_path <frequency data path> \ # required
--output_path <output folder> \ # required
--snp <SNP column name> \ # required
--a1 <A1 column name> \ # required
--a2 <A2 column name> \ # required
--stat <BETA/OR column name> \ # required
--OR \ # use this flag if sumstats are reported as odds ratios
--logit \ # use this flag if the sumstats come from logistic regression, regardless of whether OR or beta is reported
--p <P column name> \ # required
--n.total <total sample size> \ # a number, required
--n.col <N/sample size column name> \ # see below for detailed instruction for sample size
--n.case <case sample size> \ # see below for detailed instruction for sample size
--n.con <control sample size> \ # see below for detailed instruction for sample size
--n.case.col <case sample size column name> \ # see below for detailed instruction for sample size
--n.control.col <control sample size column name> # see below for detailed instruction for sample size
--chr <CHR column name> \ # optional
--bp <BP column name> \ # optional
--se <SE column number> \ # optional (recommended to provide)
--maf <MAF column number> \ # optional (recommended to provide)
```

### Sample size requirement
Sample size information can be often misspecified in reported GWAS summary statistics. Ideally GWAS sumstats contain per-SNP total sample size for linear regression association statistics and per-SNP case and control sample size for logistic regression association statistics. In practice, for linear and logistic summary statistics, users should provide one of the following sample size information respectively with priority shown below:
* **Linear regression**: ```n.col``` > ```n.total```
* **Logistic regression**: ```n.case.col; n.con.col``` > ```n.case; n.con```. If users provide ```n.case; n.con```, ```n.col``` is also recommened to provide.

If the sumstats don't contain any per-SNP sample size information, this script will impute sample size and conduct QC based on imputed sample size. We follow sample size imputation introduced in [Prive et al. (2022)](https://www.cell.com/hgg-advances/fulltext/S2666-2477(22)00052-5).

## Using PUMAS-ensemble

### Subsample training, tuning, ensemble training, and testing summary statistics
PUMAS-ensemble uses the same inputs as PUMAS. The only difference between implementation between PUMAS and PUMAS-ensemble is scripting. To partition full GWAS summary statistics to four different subsets, run:
```
Rscript ./code/PUMAS-ensemble.subsampling.R \
--k <number of folds> \
--partitions <training>,<tuning>,<ensemble training>,<testing> \
--trait_name <trait name> \
--ensemble <ensemble_method> \
--gwas_path <GWAS sumstats folder> \
--ld_path <ld folder> \
--output_path <output folder> \
--parallel \
--threads <number of parallel threads>
```
  * `k`: number of folds in PUMAS-ensemble implementation of Monte Carlo cross-validation (e.g., `--k 4`)
  * `partitions`: subsets' sample size proportion compared to total samples (e.g., `--partitions 0.6,0.2,0.1,0.1`)
  * `trait_name`: file name of GWAS summary statistics
  * `ensemble`: the method of ensemble learning you want to do. Options are `EN` (ElasticNet), `SL` (Super Learning), or `all` (linear, EN, and SL ensemble methods). If you want to use Superlearning, you must pass `SL` in the subsampling step as well to get an additional partition of the data.
  * `gwas_path`: folder containing GWAS summary statistics
  * `ld_path`: path to ld file containing approximately independent LD blocks
  * `output_path`: folder to write partitioned GWAS summary statistics
  * `parallel`: boolean flag for whether to run the Rscript in parallel or not. Default is FALSE. If you run PUMA-CUBS in parallel it will take significantly less time.
  * `threads`: the number of threads to run the Rscript in parallel with.
  
### Construct ensemble PRS and benchmark PRS models
The required input datasets are mostly the same as PUMAS's PRS evaluation function. Different from PUMAS, PUMAS-ensemble requires SNP weights from each PRS method to be stored in a separate `.<pre_method>.txt` file so that PUMAS-ensemble can construct ensemble PRS based on fine-tuned PRS model from each method and benchmark all PRS models (**again, please make sure that SNP weights files have exactly the same set of SNPs, A1, and A2 in the same order as subsampled summary statistics**). After PRS model training, run:
```
Rscript ./code/PUMAS-ensemble.evaluation.R \
--k <number of folds> \
--ref_path <LD ref> \
--trait_name <trait name> \
--prs_method <prs_methods> \
--ensemble <ensemble_method> \
--xty_path <subsampled sumstats folder> \
--stats_path <statistics folder> \
--weight_path <SNP weights> \
--full_weight_path <SNP weights from full sumstats> \
--output_path <output folder> \
--parallel \
--threads <number of parallel threads> \
--thr <threshold for stopping ensemble learning>
```
  * `k`: number of folds in PUMAS-ensemble's implementation of Monte Carlo cross-validation
  * `ref_path`: path to the LD genotype data
  * `trait_name`: file name for all subsampled summary statistics and SNP weight txt file
  * `prs_method`: PRS methods' names (e.g., put `lassosum,prscs,ldpred2` if the SNP weights files are named `Height.lassosum.txt`,`Height.prscs.txt`, and `Height.ldpred2.txt`)
  * `ensemble`: the method of ensemble learning you want to do. Options are `EN` (ElasticNet), `SL` (Super Learning), or `all` (linear, EN, and SL ensemble methods)
  * `xty_path`: folder containing partitioned summary statistics from the subsampling step
  * `stats_path`: folder containing statistics from the subsampling step (e.g, variance of phenotype and sample size for each sumstats)
  * `weight_path`: folder containing SNP weights files for prs methods
  * `full_weight_path`: folder containing SNP weights files for prs methods with full summary statistics. This will be used to calculate the ensemble weights
  * `output_path`: folder to write PRS model-tuning and benchmarking results by PUMA-CUBS
  * `parallel`: boolean flag for whether to run the Rscript in parallel or not. Default is FALSE. If you run PUMA-CUBS in parallel it will take significantly less time.
  * `threads`: the number of threads to run the Rscript in parallel with. Default is k if `parallel` is TRUE, 1 if `parallel` is FALSE.
  * `thr` (optional): threshold for stopping ensemble learning. Default is 1e-3. Don't change this value unless you have a reason to.


## Using PUMAS
### Subsample training and tuning summary statistics
For PUMAS/PUMAS-ensemble to subsample GWAS summary statistics from a full GWAS summary-level data, two datasets are requried:
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
  * `ld_path`: folder containing approximately independent LD blocks
  * `output_path`: folder to write partitioned GWAS summary statistics

### Evaluate PRS performance
After partitioning summary statistics, users can train any PRS method using the subsampled training summary statistics. Then, to use PUMAS for evaluating and fine-tuning PRS methods, three datasets are required:
  * **Remaining summary statistics** obtained from the subsampling step
  * **Genotype data from the LD panel**
  * **SNP weights**. A flat `.<pre_method>.txt` file containing SNP weights for each tuning parameter in a PRS framework (there can be multiple columns of SNP weights). **Importantly, please make sure that SNP weights files have exactly the same set of SNPs, A1, and A2 in the same order as subsampled summary statistics**. An example of SNP weights file looks like this:
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
  * `output_path`: folder to write PRS evaluation and model-tuning results by PUMAS
  
## Output
### PUMAS-ensemble
#### Subsampling
* `<trait_name>.gwas.omnibus.ite<i>.txt`: subsampled training GWAS summary statistics in the same format of input full GWAS summary statistics
* `<trait_name>.xty.omnibus.ite<i>.txt`: subsampled tuning, ensemble training, and testing summary statistics
* `<trait_name>.omnibus.forEVAL.txt`: information including variance of phenotype and each subset of summary statistics' sample size

#### Ensemble PRS construction and PRS benchmarking
* `<trait_name>.<prs_method>.en.r2.txt`: evaluated on the testing dataset. This file includes predictive R2 for ensemble PRS for each fold of Monte Carlo cross-validation from elastic net-based ensemble learning. Each row is a fold in MCCV.
* `<trait_name>.<prs_method>.sl.r2.txt`: evaluated on the testing dataset. This file includes predictive R2 for ensemble PRS for each fold of Monte Carlo cross-validation from super learning-based ensemble learning. Each row is a fold in MCCV.
* `<trait_name>.ensemble.weights.txt`: this file includes the weight of each SNP from all PRS methods, and from `EN` and/or `SL` depending on the ensemble method chosen.

### PUMAS
#### Subsampling
* `<trait_name>.gwas.ite<i>.txt`: subsampled training GWAS summary statistics in the same format of input full GWAS summary statistics
* `<trait_name>.xty.ite<i>.txt`: subsampled tuning sammary statistics
* `<trait_name>.forEVAL.txt`: information including variance of phenotype and each subet of suammry statistics' sample size

#### PRS Fine-tuning
* `<trait_name>.<prs_method>.txt`: predictive R2 for each tuning parameter within a PRS method for each fold of Monte Carlo cross-validation. Each row is a fold in MCCV and each column is a tuning parameter.

## Citation
* If you use PUMAS/PUMAS-ensemble, please cite:
  
  Zhao, Z., Gruenloh, T., Wu, Y., Sun, Z., Miao, J., Wu, Y., Song, J., & Lu, Q. (2024). [Optimizing and benchmarking polygenic risk scores with GWAS summary statistics](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03400-w). *Genome Biology*, 25, 260.

  Zhao, Z., Yi, Y., Song, J., Wu, Y., Zhong, X., Hohman, T.J., Fletcher, J., & Lu, Q. (2021). [PUMAS: fine-tuning polygenic risk scores with GWAS summary statistics](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02479-9). *Genome Biology*, 22, 257.
  
## Support
Please send questions and issues related to PUMAS/PUMAS-ensemble software to Zijie Zhao (zzhao232@wisc.edu), Stephen Dorn (svdorn@wisc.edu) and Qiongshi Lu (qlu@biostat.wisc.edu).
