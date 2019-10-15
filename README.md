# PRS-Fine-tuning
Fine-tuning polygenic risk score models using GWAS summary statistics

## Introduction

Our project gives...

### Prerequisites

The following R packages are required
##### devtools
##### PRS-Fine-tuning (Our package)

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
Downloading the following libraries on working directory

```
install_github("qlu-lab/PRS-Fine-tuning/PRS-Fine-tuning")
```


## Input Data
We requires an input file in a form of csv including those values. The file name is T0030_pruned.txt in our case.
| Parameter                   | Parameter Usage | Description                                                                  |
|----------------------------|----------------|------------------------------------------------------------------------------|
| Beta            |  data.real$Beta     |  |
| MAF         | data.real$EAF          |         |
| Se              | data.real$SE        |                         |                                                    
| Sample          |N.sample=766345 | Sample size that we are interested in |
| Res_TildeRL         |results by FunII.TildeRL |   |
| P-value        |--p_val | The cutoff value|



### Step 1

Ask users to create folder results
```
mkdir results
```

### Step 2
Run p-value fine-tuning

```

```


## Authors

See also the list of [contributors](##) who participated in this project.

## License

All rights reserved for Lu-Laboratory

## Acknowledgments


