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
Downloading PRS-Fine-tuning on working directory

In R, type the following command in console
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
| P-value        |--p_val | The cutoff value |

The input line should be in /Your-Wokring-directory/PRS-Fine-tuning/input

## Output Data
The output will be a png file, for full interpretations please see details in [wiki](https://github.com/qlu-lab/PRS-Fine-tuning/wiki) page.
![Test Image 4](https://github.com/qlu-lab/PRS-Fine-tuning/blob/master/PRS-Fine-tuning/result/T0030_pruned.png)

### Step 1
After loading the library devtools and downloading our libraries.

```
devtools::load_all(PRSFinetuning)
```
If you are having trouble with working directory, go to menu of R > Build > Load All

### Step 2
Choose whether you needs a output graph saved to /Your-Wokring-directory/PRS-Fine-tuning/output
```

```


## Authors

See also the list of [contributors](##) who participated in this project.

## License

All rights reserved for Lu-Laboratory

## Acknowledgments


