\name{GWASinlps-package}

\alias{GWASinlps-package}

\docType{package}

\title{Nonlocal prior based iterative variable selection tool for genome-wide association study data, or other high-dimensional data}

\description{The \pkg{GWASinlps} package performs variable selection for data from genome-wide association studies (GWAS), or other high-dimensional data, combining in an iterative framework, the computational efficiency of the screen-and-select variable selection approach based on some association learning and the parsimonious uncertainty quantification provided by the use of nonlocal priors as described in Sanyal et al. (2018).
}

\details{
\tabular{ll}{
Package: \tab GWASinlps\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2022-09-12\cr
License: \tab GPL (>= 2)\cr
}

The main function: \cr
\code{\link{GWASinlps}}

The main function calls the following functions: \cr
\code{\link{nlpsLM}} \cr
\code{\link{nlpsGLM}}
}

\author{Nilotpal Sanyal <nilotpal.sanyal@gmail.com>

Maintainer: Nilotpal Sanyal <nilotpal.sanyal@gmail.com>
}

\references{
Sanyal et al. (2018), "GWASinlps: Nonlocal prior based iterative SNP selection tool for genome-wide association studies". Bioinformatics, 35(1), 1-11.
}











