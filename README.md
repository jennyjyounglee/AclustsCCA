
<!-- README.md is generated from README.Rmd. Please edit that file -->

## 0. What is AclustsCCA?

The framework of `AclustsCCA` consists of two parts The goal of
AclustsCCA is to …

## 1. Installation

You can install the released version of AclustsCCA from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("AclustsCCA")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jennyjyounglee/AclustsCCA")
library(AclustsCCA)
```

## 2. Input Parameters for AclustsCCA

The below is a list of input parameters:

### 2.1 Input data

**`X`** by exposure data matrix, where is sample size and is number of
exposures.

**`Y`** by outcome data matrix, where is sample size and is number of
outcomes.

**`Z`** by confounder data matrix, where is sample size and is number of
confounders. If `NULL`, partial residuals are used for SparseCCA
analysis.

**`clusters.list`** A list of clusters with CpG sites obtained using
A-clustering, each item is a cluster that contains a set of probes.
A-clustering is implemented if `NULL` or can be provided by users.

**`annot`** A preloaded annotation file that includes columns “IlmnID”,
“Coordinate\_37”, “Islands\_Name”, “Relation\_to\_Island”,
“UCSC\_RefGene\_Name”. Only needed if **`clusters.list`** is `NULL`.

### 2.2 Input parameters associated with A-clustering:

**`dist.type`** A type of similarity distance function. Options are
“spearman” (default), “pearson” (correlation measures) or “euclid”.

**`Aclust.method`** A type of clustering function. Options are “single”,
“complete” or “average” (default).

**`thresh.dist`** A similarity distance threshold. Two neighboring
clusters are merged to a single cluster if the similarity distance
between them is above dist.thresh. Corresponds to *D̄* in the paper and
the default is 0.2

**`max.dist`** Optional maximum length between neighboring variables
permitting to cluster them together. Corresponds to *d̄*<sub>*b**p*</sub>
in the paper and the default is 1000.

**`bp.thresh.dist`** A distance in chromosomal location. Any set of
methylation sites within an interval smaller or equal to bp.dist will be
potentially merged, depending on the similarity between sites at the
ends of the interval. Corresponds to $\\underline{d}\_{bp}$ in the paper
and the default is 999.

### 2.3 Input parameters associated with SparseCCA:

**`Xmethod`** A penalty function for the exposure, i.e. penalty function
when regressing Y onto X. Options are “lasso”, “alasso”,“gglasso”, and
“SGL” (default).

**`Ymethod`** A penalty function for the outcome, i.e. penalty function
when regressing X onto Y. Options are “lasso”, “alasso”,“gglasso”,
“SGL”, and “OLS” (default).

**`init.method`** Initialization method. Options are “lasso”, “OLS”, and
“SVD” (default).

**`X.groupidx`** A vector of length that indicates grouping structure of
exposure .

**`Y.groupidx`** A vector of length that indicates grouping structure of
outcome .

**`standardize`** A logical flag for exposure and outcome
standardization, prior to fitting the model.

**`max.iter`** A maximum number of iterations of SparseCCA. The default
is 100.

**`conv`** A tolerance value for convergence of SparseCCA. The default
is 0.01.

## 1.4 Input parameters associated permutation test for AclustsCCA:

**`maxnum`** A maximal total number of permutations across all the
clusters.

**`maxB`** A maximal number of permutations for a single cluster.

**`permute.tmp.filepath`** Filepath to save intermittent permutation
results.

**`permute`** A logical flag for whether to run permutation test or not.

**`nthread`** A number of threads to parallelize permutation test and
implementation of SparseCCA across all the clusters.

**`FDR.thresh`** FDR threshold. The default is 0.05.

## 3. Implementation of AclustsCCA

The framework of `AclustsCCA` consists of two parts:

1.  Implement A-clustering on DNA methylation data

2.  Implement SparseCCA on each cluster identified by A-clustering

-   When implementing SparseCCA, partial residuals are used to adjust
    for potential confounders

-   For statistical inference, permutation test is performed

This entire process can be done using the code below, but I personally
DO NOT suggest this for computational time and memory. Please read the
next section for the suggested way of running `AclustsCCA`.

``` r
# Load annotation file
data(annot) # row: CpG sites
# Load sample data
data(sample.data)
DATA.X <- sample.data$DATA.X # row: subjects (n), column: exposures (p)
DATA.Y <- sample.data$DATA.Y # row: subjects (n), column: CpG sites (q)
DATA.Z <- sample.data$DATA.Z # row: subjects (n), column: confounders (r)

# Settings for Aclust
dist.type <- "spearman"
Aclust.method <- "average"
thresh.dist <- 0.2
max.dist <- 1000
bp.thresh.dist <- 999

# Settings for SparseCCA
Xmethod <- "SGL"
Ymethod <- "lasso"
X.groupidx <- c(rep(1,5),rep(2,5),rep(3,5),rep(4,5))
maxB <- 500
nthread <- 1

AclustsCCA.result <- AclustsCCA(X=DATA.X,
                                Y=DATA.Y,
                                Z=DATA.Z,
                                clusters.list=NULL,
                                annot=annot,
                                # parameters for A-clustering
                                dist.type = dist.type,
                                Aclust.method = Aclust.method,
                                thresh.dist = thresh.dist,
                                max.dist = max.dist,
                                bp.thresh.dist = bp.thresh.dist,
                                # parameters for SparseCCA
                                Xmethod=Xmethod,
                                Ymethod=Ymethod,
                                X.groupidx=X.groupidx,
                                # parameters for permutation test for AclustsCCA
                                permute=TRUE,
                                maxB=maxB,
                                nthread=nthread)

summary.AclustsCCA(AclustsCCA.result)
# AclustsCCA.result$permutation.result@B
# AclustsCCA.result$permutation.result
# test.idx <- which(hBH(AclustsCCA.result$permutation.result@g / AclustsCCA.result$permutation.result@num,0.05))
# true.idx <- which(sapply(AclustsCCA.result$clusters.list,function(x) sum(x%in%TRUE.CpGs))!=0)
# lapply(AclustsCCA.result$ALPHA.observed[test.idx],function(x) names(x)[x!=0])
# lapply(AclustsCCA.result$BETA.observed[test.idx],function(x) names(x)[x!=0])
```

## 4. Suggestions on how to run AclustsCCA

The framework of `AclustsCCA` consists of two parts:

1.  Implement A-clustering on DNA methylation data

2.  Implement SparseCCA on each cluster identified by A-clustering

-   When implementing SparseCCA, partial residuals are used to adjust
    for potential confounders

-   For statistical inference, permutation test is performed

All of the above procedures can be implemented at once using
`AclustsCCA` function. Among these steps, permutation test is definitely
a part that takes the most computational time to run and requires
attention when running.

Therefore, I personally suggest running each part separately to save
computational time and memory.

### 4.1 Implement A-clustering and save it

A-clustering is implemented to identify cluster and this part have to be
ran only once. I personally suggest to run this part separately and save
it for future use.

``` r
# Implement A-clustering
all.clusters.list <- Aclust::assign.to.clusters(betas = t(DATA.Y),
                                                annot = annot,
                                                dist.type = dist.type,
                                                method = Aclust.method,
                                                thresh.dist = thresh.dist,
                                                max.dist = max.dist,
                                                bp.thresh.dist = bp.thresh.dist)
# We only need clusters with at least two probes
clusters.list <- all.clusters.list[sapply(all.clusters.list,length)!=1]
```

### 4.2 Obtain partial residuals to adjust for potential confounders

When implementing SparseCCA, partial residuals are used to adjust for
potential confounders. Again, computing partial residuals have to be ran
only once. I personally suggest to run this part separately and save it
for future use.

``` r
X.resid <- partial.residual(data=DATA.X,Z=DATA.Z,nthread=1)
Y.resid <- partial.residual(data=DATA.Y,Z=DATA.Z,nthread=1)
```

### 4.3 Implement SparseCCA and permutation test on each cluster identified by A-clustering

As list of clusters are provided and partial residuals are computed,
`AclustsCCA` will only run SparseCCA on each cluster identified by
A-clustering.

``` r
AclustsCCA.result <- AclustsCCA(X=X.resid,
                                Y=Y.resid,
                                clusters.list=clusters.list,
                                # parameters for SparseCCA
                                Xmethod=Xmethod,
                                Ymethod=Ymethod,
                                X.groupidx=X.groupidx,
                                # parameters for permutation test for AclustsCCA
                                maxB=maxB,
                                permute=TRUE,
                                nthread=nthread) 
```
