# PUMAS
Fine-tuning polygenic risk score models using GWAS summary statistics

## Updates

Last update: 10/25/2019. Please update the package if downloaded before 10/25/2019 5:00am.

## Introduction

Our project gives user-friendly function that presents a direct and explicit result of fine-tuning polygenic risk score models. For full explanation and tutorial, please go to  [wiki](https://github.com/qlu-lab/PRS-Fine-tuning/wiki)

### Prerequisites
R
The following R packages are required
`devtools`
`PRS-Fine-tuning` 

If you don't have devtools installed, please install as following
##### Mac and Linux:
```
devtools::install_github("hadley/devtools")
```
##### Windows:
```library(devtools)
build_github_devtools()

#### Restart R before continuing ####
install.packages("devtools.zip", repos = NULL)

# Remove the package after installation
unlink("devtools.zip")
```
Then, load the package
```
library(devtools)
```

### Step 0
Downloading PUMAS on working directory

In R, type the following command in console
```
install_github("qlu-lab/PUMAS")
```




### Step 1
After loading the library devtools and downloading our libraries.

```
devtools::load_all(PUMAS)
```
If you are having trouble with working directory, after downloading the package from github, you can double click on `pumas.Rproj` and go to menu of R > Build > Load All

### Step 2
Giving input infortion to do an analysis from our package in R console. In our example, see details in [wiki](https://github.com/qlu-lab/PUMAS/wiki) 
```
pumas.main("/working-directory-to-input/T0030_pruned.txt","/working-directory-to-output/T0030_pruned.png","Beta","EAF","SE",766345,TRUE)
```

# Quick Start

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
| n_fold          |4 | Number of subsets |

## Optional Features

`samplesize_header` can be either a number or a character of column name. If the input GWAS does not include a column for per-SNP sample sizes, the user can provide a single number (usually reported along with a published GWAS) as the uniform sample size for all SNPs.

`n_fold` can be user-specified. When `n_fold` is not specified, PUMAS will use a default of 4 subsets to implement the model tuning approach.

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

All rights reserved for [Lu-Laboratory](http://qlu-lab.org/)

## Citation

If you use the package, please cite

[Zhao, Z. et al. Fine-tuning Polygenic Risk Scores with GWAS Summary Statistics. bioRxiv, 810713 (2019).](https://www.biorxiv.org/content/10.1101/810713v1)

## Acknowledgments


