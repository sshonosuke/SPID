# SPID
This package implements Bayesian spatial estimaton for area-wise income distributions from grouped data, as proposed by the following papers.

Sugasawa, S., Kobayashi, G. and Kawakubo, Y. (2019). Estimation and inference for area-wise spatial income distributions from grouped data. https://arxiv.org/abs/1904.11109

Functions are implemented in SPID-function.R available in the repository.

```{r}
source("SPID-function.R")   # require "MCMCpack", "sparseM" and "statmod" packages
```

Load demo dataset
```{r}
load("demo-data.RData")
```

Fit the proposed models (PWD and PWL)

Input of `PWD` and `PWL`

- `Data`: (m,N)-matrix of observed counts (m: number of areas; N: number of groups)
- `Z`: vector of boudary values for grouping 
- `W`: adjacent matrix
- `mcmc`: length of MCMC 
- `burn`: burn-in period
- `print`: Number of iterations of MCMC is shown if `T`

Output of `PWD`: List object of MCMC results

- `Beta`: regression coeffieicnts in the outcome model
- `Sigma`: squared value of error variance in the outcome model
- `Phi`: coefficients of the polynomial term of `Y` in the response model
- `Gamma`: coefficients of the spline term of `Y` in the response model
- `Delta`: coefficients for the covaraites in the response model
- `Lam`: precision parameter in the prior for `Gamma`
- `a`: scale parameter controlling locations of knots

```{r}
set.seed(1)
qq=2
K=10
mc=7000
bn=2000
fit=BSS.LM(Y,X,Z,S,q=qq,K=K,mc=mc,burn=bn)
```



