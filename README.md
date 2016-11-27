# BhGLM - NBMM (Negative Binomial Mixed Model)

## Introduction
Recent advances in next-generation sequencing (NGS) technology enable researchers to collect a large volume of metagenomic sequencing data. These data provide valuable resources for investigating interactions between the microbiome and host environmental/clinical factors. In addition to the well-known properties of microbiome count measurements, for example, varied total sequence reads across samples, over-dispersion and zero-inflation, microbiome studies usually collect samples with hierarchical structures, which introduce correlation among the samples and thus further complicate the analysis and interpretation of microbiome count data. We propose negative binomial mixed models (NBMMs) for detecting the association between the microbiome and host environmental/clinical factors for correlated microbiome count data. The proposed mixed-effects models incorporate random effects into the commonly used fixed-effects negative binomial model to account for correlation among the samples. 

The details of the statistical model are as follows:

<img src="image/NBMM.PNG" width="600" align="center">

## Installation
You can install our BhGLM package by downloading BhGLM_1.1.0.zip
```r
install.packages("./BhGLM")
library(BhGLM)
```

## Contact
Feel free to contact me by nyi AT uab.edu
