# creNet
creNET is an R package developed for prediction of response or risk factors using base-line gene expession data. The method is a  biological network-based regularization model for regression and classification. The model uses and ''ovelap group-norm penalty'' to penalize gene sets rather than individual genes. 

Gene ssets are defined based on the topology of the gene interaction network. The penalty term is weighted to inforce differential shrinkage of gene sets based on their biological relevance. The biological relevance is measured using the R package QuaternaryProd: https://github.com/carltonyfakhry/QuaternaryProd

## Installation
To install the package directly from github, you need the devtools libray.
```{R}
library(devtools)
install_github("kouroshz/creNet", local = FALSE)
```
## Usage
Please see creNetScripts for example usage: https://github.com/kouroshz/creNetScripts
