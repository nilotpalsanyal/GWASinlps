
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" align="right" style="float:right; height:129px; padding-top:3px;" />

# GWASinlps: Nonlocal Prior Based Iterative Variable Selection Tool for Genome-Wide Association Studies

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/GWASinlps)](https://CRAN.R-project.org/package=GWASinlps)
[![CodeFactor](https://www.codefactor.io/repository/github/nilotpalsanyal/gwasinlps/badge)](https://www.codefactor.io/repository/github/nilotpalsanyal/gwasinlps)
[![R-CMD-check](https://github.com/nilotpalsanyal/GWASinlps/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nilotpalsanyal/GWASinlps/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

GWASinlps performs Bayesian nonlocal prior based iterative variable
selection for data from genome-Wide association studies (GWAS), or other
high-dimensional data, as described in Sanyal & Ferreira
([2012](#ref-paper)).

## Installation

### Install from CRAN

``` r
install.packages("GWASinlps")
```

### Install from GitHub

``` r
# install.packages("devtools")
devtools::install_github("nilotpalsanyal/GWASinlps")
```

## The main function:

`GWASinlps()` is the main function which accepts continuous or binary
data (such as phenotype data) and a matrix with the independent variable
values (SNP genotypes). The function also needs as input values for
scaling parameter of the selected nonlocal prior and the tuning
paramters. These should be fixed based on exploratory study and/or
subject-specific heuristics. For example, in GWAS analysis, as the GWAS
effect sizes are generally very small (typical effect size of a SNP is
around 0.05% of the total phenotypic variance for quantitative traits),
the scaling parameter can be chosen such that the non-local prior allows
at least 1% chance of a standardized effect size being 0.05 or less in
absolute value. Such estimates of the scaling parameter for the MOM and
iMOM priors are 0.022 and 0.008, respectively.

Here is a simple illistration of the use the `GWASinlps()` function for
both continous and binary phenotypes.

<!-- For a detailed manual, see the package [vignette](https://cran.r-project.org/web/packages/BHMSMAfMRI/vignettes/BHMSMAfMRIvignette.pdf){target="_blank"}.  -->

### GWASinlps analysis with continuous data/phenotypes

``` r
library(GWASinlps)
#> Loading required package: mombf
#> Loading required package: mvtnorm
#> Loading required package: ncvreg
#> Loading required package: mgcv
#> Loading required package: nlme
#> This is mgcv 1.8-40. For overview type 'help("mgcv-package")'.
#> Loading required package: fastglm
#> Loading required package: bigmemory
#> 
#>  Welcome to GWASinlps...Happy selection!
#>  
#>  Website: https://nilotpalsanyal.github.io/GWASinlps/
#>  Bug report: https://github.com/nilotpalsanyal/GWASinlps/issues

# Generate design matrix (genotype matrix)
n = 200   #number of subjects
p = 10000 #number of variables/SNPs
m = 10    #number of true variables/SNPs

set.seed(1) 
f = runif( p, .1, .2 )  #simulate minor allele frequency
x = matrix(nrow = n, ncol = p)
for(j in 1:p)
  x[,j] = rbinom(n, 2, f[j])  #simulate genotypes
colnames(x) = 1:p

# Generate true effect sizes
causal_snps = sample(1:p, m)
beta = rep(0, p)
beta[causal_snps] = rnorm(m, mean = 0, sd = 2 )

# Generate continuous (phenotype) data
y = x %*% beta + rnorm(n, 0, 1) 

# GWASinlps analysis
inlps <- GWASinlps(y, x, family="normal", prior="mom", tau=0.2, 
          k0=1, m=50, rxx=0.2)
#> =================================
#> Number of selected variables: 9
#> Time taken: 0.04 min
#> =================================

# Results
cat( "GWASinlps with continuous phenotypes selected", length(inlps$selected), 
"SNPs with", length(intersect(inlps$selected, causal_snps)), 
"true positive(s) and", length(setdiff(causal_snps, inlps$selected)), "false 
negative(s) out of a pool of", p, "SNPs with data from", n, "persons." )
#> GWASinlps with continuous phenotypes selected 9 SNPs with 8 true positive(s) and 2 false 
#> negative(s) out of a pool of 10000 SNPs with data from 200 persons.

# Compare with LASSO
library(glmnet)
#> Loading required package: Matrix
#> Loaded glmnet 4.1-4
fit.cvlasso = cv.glmnet( x, y, alpha = 1 )
l.min = fit.cvlasso $lambda.min # lambda that gives minimum cvm
l.1se = fit.cvlasso $lambda.1se  # largest lambda such that error is 
                                 # within 1 se of the minimum

lasso_min = which( as.vector( coef( fit.cvlasso, s = l.min ) )[-1] != 0 )
cat( "LASSO with lambda.min selected", length(lasso_min), "SNPs with", 
     length(intersect(lasso_min, causal_snps)), "true positives and",
    length(setdiff(causal_snps, inlps$selected)), "false negative(s)." )
#> LASSO with lambda.min selected 190 SNPs with 8 true positives and 2 false negative(s).

lasso_1se = which( as.vector( coef( fit.cvlasso, s = l.1se ) )[-1] != 0 )
cat( "LASSO with lambda.1se selected", length(lasso_1se), "SNPs with", 
     length(intersect(lasso_1se, causal_snps)), "true positives and",
    length(setdiff(causal_snps, inlps$selected)), "false negative(s)." )
#> LASSO with lambda.1se selected 44 SNPs with 8 true positives and 2 false negative(s).
```

### GWASinlps analysis with binary data/phenotypes

``` r
library(GWASinlps)

# Generate design matrix (genotype matrix)
n = 200   #number of subjects
p = 1000 #number of variables/SNPs
m = 10    #number of true variables/SNPs

set.seed(1) 
f = runif( p, .1, .2 )  #simulate minor allele frequency
x = matrix(nrow = n, ncol = p)
for(j in 1:p)
  x[,j] = rbinom(n, 2, f[j])  #simulate genotypes
colnames(x) = 1:p

# Generate true effect sizes
causal_snps = sample(1:p, m)
beta = rep(0, p)
beta[causal_snps] = rnorm(m, mean = 0, sd = 2 )

# Generate binary (phenotype) data
prob = exp(x %*% beta)/(1 + exp(x %*% beta))
y = sapply(1:n, function(i)rbinom(1,1,prob[i]) )

# GWASinlps analysis
mode(x) = "double"  #needed for fastglm() function below
mmle_xy = apply( x, 2, function(z) coef( fastglm(y=y, 
  x=cbind(1,matrix(z)), family = binomial(link = "logit")) )[2] ) 
  #pre-compute MMLEs of betas as it takes time

inlps <- GWASinlps(y, x, mmle_xy=mmle_xy, family="binomial", 
prior="mom", tau=0.2, k0=1, m=50, rxx=0.2)
#> =================================
#> Number of selected variables: 3
#> Time taken: 0.1 min
#> =================================

# Results
cat( "GWASinlps with binary phenotypes selected", 
length(inlps$selected), "SNPs with", 
length(intersect(inlps$selected, causal_snps)), "true positive(s) and", 
length(setdiff(causal_snps, inlps$selected)), 
"false negative(s) out of a pool of", p, "SNPs with data from", n, "persons." )
#> GWASinlps with binary phenotypes selected 3 SNPs with 3 true positive(s) and 7 false negative(s) out of a pool of 1000 SNPs with data from 200 persons.
```

## References:

<div id="refs" class="references">

<div id="ref-paper">

Sanyal et al. (2018), “GWASinlps: Nonlocal prior based iterative SNP
selection tool for genome-wide association studies”. Bioinformatics,
35(1), 1-11. <span
target="_blank"><https://doi.org/1093/bioinformatics/bty472></span>
