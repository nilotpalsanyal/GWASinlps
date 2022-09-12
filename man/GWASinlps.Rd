\name{GWASinlps}

\alias{GWASinlps}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Nonlocal prior based iterative variable selection for GWAS data, or other high-dimensional data}

\description{\code{GWASinlps} performs variable selection with data from Genome-wide association studies (GWAS), or other high-dimensional data, combining in an iterative framework, the computational efficiency of the screen-and-select variable selection approach based on some association learning and the parsimonious uncertainty quantification provided by the use of nonlocal priors as described in Sanyal et al. (2018).}

\usage{
GWASinlps(
  y, 
  x, 
  family = c("normal","binomial"),
  cor_xy = NULL,
  mmle_xy = NULL,
  prior = c("mom", "imom", "zellner", "horeshoe"), 
  tau, 
  priorDelta = modelbbprior(1,1), 
  k0, 
  m, 
  rxx, 
  nskip = 3, 
  niter = 2000, 
  verbose = FALSE, 
  seed = NULL, 
  tau.hs.method = "halfCauchy", 
  sigma.hs.method = "Jeffreys"
)
}

\arguments{
  \item{y}{The vector of continuous or binary response (phenotype).}
  \item{x}{The design matrix with subjects in rows and independent variables (SNPs) in columns. Missing values are not accepted currently.}
  \item{family}{\code{"normal"} for continuous data, \code{"binomial"} for binary data (logit link is used).}
  \item{cor_xy}{Used only when \code{family="normal"}. Vector of (Pearson) correlation coefficients of \code{y} with the individual columns of \code{x}.}
  \item{mmle_xy}{Used only when \code{family="binomial"}. Vector of maximum marginal likelihood estimates of the regression parameters corresponding to the \code{x} variables (SNPs). These are obtained from GLM fits of \code{y} with the individual columns of \code{x} including an intercept.}
  \item{prior}{\code{"mom"} for pMOM prior, \code{"imom"} for piMOM prior, \code{"zellner"} for Zellner's g-prior, \code{"horseshoe"} for horseshoe prior. For \code{family="binomial"}, \code{"zellner"} considers group Zellner prior and \code{"horseshoe"} is not available.}
  \item{tau}{The value of the scale parameter tau of the nonlocal prior. }
  \item{priorDelta}{Prior for model space. Defaults to \code{modelbbprior(1,1)}, which is beta-binomial(1,1) prior.}
  \item{k0}{GWASinlps tuning parameter denoting the number of leading SNPs (see Details).}
  \item{m}{GWASinlps tuning parameter, denoting the maximum number of SNPs to be selected.}
  \item{rxx}{GWASinlps tuning parameter denoting the correlation threshold to determine leading sets (see References)}
  \item{nskip}{GWASinlps tuning parameter denoting the maximum allowed count of skipping an iteration that does not select any variable (SNP) (see References).}
  \item{niter}{Number of MCMC iterations for nonlocal prior based Bayesian variable selection. Defaults to 2000.}
  \item{verbose}{If TRUE, prints result from the iterations progressively. FALSE by default.}
  \item{seed}{For reproducibility. If provided, the random seed is set to this value at the beginning of the function.}
  \item{tau.hs.method}{Necessary only when \code{prior="horseshoe"}. See \code{horseshoe} function reference.}
  \item{sigma.hs.method}{Necessary only when \code{prior="horseshoe"}. See \code{horseshoe} function reference.}
}

\details{The GWASinlps method selects variables (SNPs) iteratively. 

For continuous response (phenotype), the procedure starts with an initial set of independent variables (SNPs), a design matrix (SNP genotype matrix) \code{x} and a response (phenotype) vector \code{y}. An iteration proceeds by determining the \code{k0} \emph{leading SNPs} having the highest association with \code{y}. The measure of association is the absolute value of the Pearson's correlation coefficient \code{cor_xy}. These \code{k0} leading SNPs, in turn, determine \code{k0} \emph{leading sets}, where each leading set consists of all SNPs with absolute correlation coefficient more than or equal to \code{rxx} with the correspondng leading SNP. Then within each leading set, non-local prior based Bayesian variable selection is run (using package \pkg{mombf}) and the variables (SNPs) appearing in the highest posterior probability model (HPPM) are considered selected in the current iteration. Note that a single variable (SNP) can be selected from multiple leading sets. The selected variables (SNPs) are regressed out from \code{y} using \code{lm} fit. The variables (SNPs) that are included in one or more \emph{leading sets} but do not appear in any HPPM are dropped from further analysis. With updated \code{y} and variable (SNP) set, next iteration proceeds similarly. The procedure continues until the stopping point, determined by the GWASinlps tuning parameters \code{m}, \code{rxx}, and \code{nskip}, is reached. For more details, see References.

For binary response (phenotype), the procedure starts with an initial set of variables (SNPs), a design matrix (SNP genotype matrix) \code{x} and a response (phenotype) vector \code{y}. The first iteration proceeds by determining the \code{k0} \emph{leading SNPs} having the highest association with \code{y}, where the measure of association is the absolute value of the maximum marginal likelihood estimate \code{mmle_xy}. These \code{k0} leading SNPs, in turn, determine \code{k0} \emph{leading sets}, where each leading set consists of all SNPs with absolute correlation coefficient more than or equal to \code{rxx} with the correspondng leading SNP. Then within each leading set, non-local prior based Bayesian variable selection is run (using package \pkg{mombf}) and the variables (SNPs) appearing in the HPPM are considered selected in the first iteration. Note that a single variable (SNP) can be selected from multiple leading sets. The variables (SNPs) which are included in one or more \emph{leading sets} but do not appear in any HPPM are dropped from further analysis. Unlike the case of continuous data, for binary data the selected variables (SNPs) cannot be regressed out from \code{y} as the definition of residual is not as straightforward. Here, we account for the selected variables (SNPs) by including them in \code{glm} fits of \code{y} with each of the remaining variables (SNPs). The \code{glm} coefficients of the remaining variables, thus obtained, reflect their contribution in presence of the selected variables (SNPs) of the first iteration. Considering the absolute values of these \code{glm} coefficients as the measure of association, we proceed with the second iteration with updated variable (SNP) set. And so on in this manner. The procedure continues until the stopping point, determined by the GWASinlps tuning parameters \code{m}, \code{rxx}, and \code{nskip}, is reached.

For horseshoe prior, package \pkg{horseshoe} is used.
}

\value{
  A list containing
    \item{selected }{Vector with names of the GWASinlps selected variables (SNPs) in the order they were selected. }
    \item{selected_iterwise }{List with selected variables (SNPs) from each iteration. }
}

\references{
Sanyal et al. (2018), "GWASinlps: Nonlocal prior based iterative SNP selection tool for genome-wide association studies". Bioinformatics, 35(1), 1-11.
}


\seealso{
\code{\link{nlpsLM}}, \code{\link{nlpsGLM}}, \code{\link[mombf]{modelSelection}}
} 


\author{Nilotpal Sanyal <nilotpal.sanyal@gmail.com>
}


\examples{
n = 200
p = 1000
m = 10

# Generate design matrix (genotype matrix)
set.seed(1) 
f = runif( p, .1, .2 ) # simulate minor allele frequency
x = matrix( nrow = n, ncol = p )
colnames(x) = 1:p
for(j in 1:p)
  x[,j] = rbinom( n, 2, f[j] )

# Generate true effect sizes
causal_snps = sample( 1:p, m )
beta = rep( 0, p )
set.seed(1)
beta[causal_snps] = rnorm(m, mean = 0, sd = 2 )

# Generate continuous (phenotype) data
y = x \%*\% beta + rnorm(n, 0, 1) 

# Fix scale parameter tau 
tau = 0.2

# GWASinlps analysis
inlps = GWASinlps(y, x, family="normal", prior="mom", tau=tau, k0=1, m=50, 
        rxx=0.2)
cat( "GWASinlps selected", length(inlps$selected), "SNPs with", 
    length(intersect(inlps$selected, causal_snps)), "true positive(s) and",
    length(setdiff(causal_snps, inlps$selected)), "false negative(s) out  
    of a pool of", p, "SNPs with data from", n, "persons." )

# Compare with LASSO
library(glmnet)
fit.cvlasso = cv.glmnet( x, y, alpha = 1 )
l.min = fit.cvlasso $lambda.min # lambda that gives minimum cvm
l.1se = fit.cvlasso $lambda.1se  # largest lambda such that error is 
                                 # within 1 se of the minimum

lasso_min = which( as.vector( coef( fit.cvlasso, s = l.min ) )[-1] != 0 )  
cat( "LASSO with lambda.min selected", length(lasso_min), "SNPs with", 
     length(intersect(lasso_min, causal_snps)), "true positives and",
    length(setdiff(causal_snps, inlps$selected)), "false negative(s)." )

lasso_1se = which( as.vector( coef( fit.cvlasso, s = l.1se ) )[-1] != 0 )
cat( "LASSO with lambda.1se selected", length(lasso_1se), "SNPs with", 
     length(intersect(lasso_1se, causal_snps)), "true positives and",
    length(setdiff(causal_snps, inlps$selected)), "false negative(s)." )
}

