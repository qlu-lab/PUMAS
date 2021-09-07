# PUMAS
Fine-tuning polygenic risk score models using GWAS summary statistics

## Updates

Please clone the repo to local in order to build and load.

Last update: 10/25/2019. Please update the package if downloaded before 10/26/2019 4:53pm CT.

## Introduction

Our project gives user-friendly function that presents a direct and explicit result of fine-tuning polygenic risk score models. A quick start panel that walks through PUMAS usage is included below. For full explanation and tutorial, please go to  [wiki](https://github.com/qlu-lab/PRS-Fine-tuning/wiki)



### Step 0
Downloading pumas repo form this page by clicking `Clone or Download`



### Step 1
Load our package by double click on `pumas.Rproj` and go to menu of R > Build > Load All

### Step 2
Giving input information to do an analysis from our package in R console. In our example, see details in [wiki](https://github.com/qlu-lab/PUMAS/wiki) 
```
pumas.main(input_path="/working-directory-to-input/T0030_pruned.txt",output_path="/working-directory-to-output/T0030_pruned.png",beta_header="Beta",af_header="EAF",se_header="SE",pvalue_header="Pval",samplesize_header=766345,make_plot=TRUE)
```

# Quick Start

## Note: before applying PUMAS to GWAS summary statistics, the GWAS data needs to be pruned in advance. Do not use clumped GWAS as input.
LD-pruning can be done by [PLINK](https://www.cog-genomics.org/plink/1.9/ld).


## Input Data
The function requires an input GWAS summary statistics file in the form of .txt/.txt.gz. The GWAS summary data should include:

| Parameter                   | Example | Description                                                                  |
|----------------------------|----------------|------------------------------------------------------------------------------|
| input_path            | "input/T0030_pruned.txt"    |  GWAS path |
| beta_header            |  Beta    |  Effect size |
| af_header         | EAF         |    Either allele/minor allele frequency      |
| se_header              | SE        |       Standard error             |    
| pvalue_header              | Pval        |       P-value             |   
| samplesize_header         |766345 | Sample size |


## Optional Features

`samplesize_header` can be either a number or a character of column name. If the input GWAS does not include a column for per-SNP sample sizes, the user can provide a single number (usually reported along with a published GWAS) as the uniform sample size for all SNPs.

`n_fold` is the number of subsets that PUMAS partitions the complete GWAS into. For example, when `n_fold=4` PUMAS will generate the training summary statistics based on 3 subsets and calculate the testing summary statistics of the remaining 1 subset. `n_fold` can be user-specified. When `n_fold` is not specified, PUMAS will use a default of 4 subsets to implement the model tuning approach.

`odds_ratio` is a boolean value that tells PUMAS whether the weight input from GWAS is effect size (quantitative) or odds ratio (binary). When `odds_ratio=T`, PUMAS asusmes that `beta_header` is the header for odds ratio and applies log-transformation on the column. The default is `odds_ratio=F`.

## Make Plot (optional)

PUMAS can output a scatterplot that illustrates a detailed pattern of the predictive performance under each PRS models. Y-axis is the predictive R2 and X-axis is log-transformed p-value cutoff for every model. To make this plot, the user needs to provide 2 more parameters to the main function:

`output_path` is the path to store the plot in `.png` format. For example, `output_path='result/test.png'` is acceptable.

`make_plot` is a boolean value that tells PUMAS whether to output a scatterplot. The default is `make_plot=F`.

For more detailed interpretations please see details in [wiki](https://github.com/qlu-lab/PRS-Fine-tuning/wiki) page.
![Test Image 4](https://github.com/qlu-lab/PUMAS/blob/master/result/T0030_pruned.png)

## Result

PUMAS results are printed in the interface.

`Maximal.R2`is the estimated R2 at for the best model selected.

`Optimal.P.value` is the p-value cutoff for the best model selected.


## Authors

See also the list of [contributors](##) who participated in this project.

## License

This package is under MIT license.

All rights reserved for [Lu-Laboratory](http://qlu-lab.org/)

## Citation

If you use the package, please cite

[Zhao, Z., Yi, Y. et al. PUMAS: fine-tuning polygenic risk scores with GWAS summary statistics. Genome Biology. 2021;22:257.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02479-9)

Cite the code: [![DOI](https://zenodo.org/badge/209808756.svg)](https://zenodo.org/badge/latestdoi/209808756)

## Acknowledgments


