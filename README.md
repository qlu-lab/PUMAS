# PRS-Fine-tuning
Fine-tuning polygenic risk score models using GWAS summary statistics

## Introduction

Our project gives...

### Prerequisites

The following R packages are required:


```
Give examples
```

### Step 0
Downloading the following libraries and create a directory called trio-twas and change to this destination directory.

```
wget ..
mkdir prs-tuning
cd prs-tuning
```
Download xx.gz and decompress xx.gz under folder trio-twas
```
unzip xx.gz

```

Change to Scripts folder 
```
cd scripts
```


End with an example of getting some data out of the system or using it for a little demo

## Input Data
| Parameter                   | Directory Name | Description                                                                  |
|----------------------------|----------------|------------------------------------------------------------------------------|
| Beta            |  data.real$Beta     |  |
| MAF         | data.real$EAF          | The chromosome of the gene        |
| Se              | data.real$SE        | Gene name                        |                                                    
| Sample          |N.sample=766345 | The path of the fam file, the ids in the fam file should fit the ids from the vcf file  |
| Res_TildeRL         |results by FunII.TildeRL | The path of the original prediciton matrix  |
|P-value        |--p_val | The cutoff value|

Explain how to run the automated tests for this system

### Step 1

Ask users to create folder results
```
mkdir results
```

### Step 2
Run p-value fine-tuning

```
R 
submit.R
```


## Authors

See also the list of [contributors](##) who participated in this project.

## License

All rights reserved for Lu-Laboratory

## Acknowledgments


