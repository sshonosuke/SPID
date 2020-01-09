# SPID
This package implements Bayesian spatial estimaton for area-wise income distributions from grouped data, as proposed by the following papers.

Sugasawa, S., Kobayashi, G. and Kawakubo, Y. (2020). Estimation and inference for area-wise spatial income distributions from grouped data. Computational Statistics & Data Analysis, to appear  (https://arxiv.org/abs/1904.11109)

Functions are implemented in SPID-function.R available in the repository.

```{r}
source("SPID-function.R")   # require "MCMCpack", "sparseM" and "statmod" packages
```

Load demo dataset
```{r}
load("demo-data.RData")
```

Fit the proposed models (PWD and PWL)

Input of `PWD` 

- `Data`: (m,N)-matrix of observed counts (m: number of areas; N: number of groups)
- `Z`: vector of boudary values for grouping 
- `W`: adjacent matrix
- `mcmc`: length of MCMC 
- `burn`: burn-in period
- `print`: Number of iterations of MCMC is shown if `T`

Output of `PWD`: List object of MCMC results

- `U`: (mc,m,p)-array of posterior samples of area-wise parameters (mc: number of posterior samples; p: dimensinon of area-wise parameter)
- `Mu`: (mc,p)-matrix of posterior samples of grand means
- `Mu`: (mc,p)-matrix of posterior samples of grand means
- `Tau`: (mc,p)-matrix of posterior samples of precision parameters in independent random effect part
- `Lam`: (mc,p)-matrix of posterior samples of precision parameters in spatial difference shrinkage part
- `ML`: (m,p)-matrix of area-wise maximum likelihood estimates 
- `Hessian`: (m,p,p)-array of Hessian matrices for maximum likelihood estimates

```{r}
mc=1500   # length of MCMC
bn=500   # burn-in period

set.seed(1)
fit1=PWD(Data,Z,W,mc,bn)
```

Input of `PWL` 

- `Data`: (m,N)-matrix of observed counts (m: number of areas; N: number of groups)
- `Z`: vector of boudary values for grouping 
- `W`: adjacent matrix
- `mcmc`: length of MCMC 
- `burn`: burn-in period
- `rn`: number of Monte Carlo samples to approximate normalizing constants (default=10)
- `print`: Number of iterations of MCMC is shown if `T`

Output of `PWL`: List object of MCMC results

- `U`: (mc,m,p)-array of posterior samples of area-wise parameters (mc: number of posterior samples; p: dimensinon of area-wise parameter)
- `Mu`: (mc,p)-matrix of posterior samples of grand means
- `Mu`: (mc,p)-matrix of posterior samples of grand means
- `Tau`: (mc,p)-matrix of posterior samples of precision parameters in independent random effect part
- `Lam`: Vector of posterior samples of precision parameter in spatial difference shrinkage part
- `ML`: (m,p)-matrix of area-wise maximum likelihood estimates 
- `Hessian`: (m,p,p)-array of Hessian matrices for maximum likelihood estimates

```{r}
fit2=PWL(Data,Z,W,mc,bn)
```

We can also apply independent random effects model as a submodel of the proposed models.  
Input and Output are almost the same as `PWD`.
```{r}
fit3=IRE(Data,Z,mc,bn)
```

Posterior means of area-wise parameters 
```{r}
apply(fit1$U,c(2,3),mean)  # PWD
apply(fit2$U,c(2,3),mean)  # PWL   
apply(fit3$U,c(2,3),mean)  # IRE
fit1$ML       # area-wise maximum likelihood  
```

Posterior means of precision parameters (Tau) 
```{r}
apply(fit1$Tau,2,mean)   
apply(fit2$Tau,2,mean)
apply(fit3$Tau,2,mean)
```

Posterior means of precision parameters (Lambda) 
```{r}
apply(fit1$Lam,2,mean)
mean(fit2$Lam)
```






