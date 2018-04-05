---
title: NBZIMM NBMMs
nav: true
---

# NBZIMM - NBMM (Negative Binomial Mixed Model)

## Installation
You can install our NBZIMM package by downloading NBZIMM_1.0.zip.
```r
install.packages("NBZIMM")
library(NBZIMM)
```

## Usage
```r
glmm.nb(fixed, random, data, subset, correlation, weights, control, niter = 30, epsilon = 1e-05, verbose = TRUE)
```
## Arguments

- **fixed, random, data, subset, correlation, weights, control**: These arguments are the same as in the function lme in the package nlme.
- **niter**: maximum number of iterations.
- **epsilon**: positive convergence tolerance.
- **verbose**: logical. If TRUE, print out number of iterations and computational time.

## Examples

We use the following simple example to show the use of NBMM.
```r
library(NBZIMM)

# parameter settings
  n = 200    
  n.dam = 20   
  b0 = c(0.4, 0.55)
  cor = c(0.5, 0.8)

# data simulation
  corr = runif(1, cor[1], cor[2])
  x = sim.x(n = n, m = 2, corr = corr)
  q = rep(1/n.dam, n.dam-1)
  q = cumsum(q)
  quantiles = quantile(x[,1], q)
  dam = as.numeric( factor(cut(x[,1], breaks = c(-Inf, quantiles, Inf))) )   
  quantiles = quantile(x[,2], 0.45)
  diet = as.numeric( factor(cut(x[,2], breaks = c(-Inf, quantiles, Inf))) )
  diet = diet - 1   
  
  da = rep(NA, n.dam)
  sigma = runif(1, 0.5, 1)
  for (j in 1:n.dam) da[j] = rnorm(1, 0, sigma)
  mu0 = runif(n, 0.1, 3.5)
  theta = runif(1, 0.1, 5) 
  b = runif(1, b0[1], b0[2]) 
  ys = sim.y(x = diet, mu = mu0 + da[dam], sigma = 1, coefs = b, p.neg = 0, nb.theta = theta) 
  y0 = ys$y.nb
  N = exp(mu0)

# model fitting and summary
  f = glmm.nb(y0 ~ offset(log(N)) + diet, random = ~ 1 | dam, verbose = F) 
  out = summary(f)
  out
```

## Simulation Studies
The following R code is used for simulation studies in a manuscript in revision and the citation will added later.

```r
rm(list=ls())

library(NBZIMM)
library(BhGLM)
library(nlme)

## parameter settings ###
n = 200       # 400, 200
n.dam = 20    # 40, 20

#b0 = c(0, 0)                    
#b0 = c(0.2, 0.35)
b0 = c(0.4, 0.55)

#cor = c(-0.1, 0.1)      
#cor = c(0.5, 0.8) 
cor = c(-0.8, -0.5)

f1 = f2 = f3 = list()
out1 = out2 = out3 = NULL
v0 = v = vresi = NULL
th0 = th = NULL
bb = NULL

#### simulation ###
n.sims = 5000 
start.time <- Sys.time()
for (i in 1:n.sims){
  
  corr = runif(1, cor[1], cor[2])
  x = sim.x(n = n, m = 2, corr = corr)
  q = rep(1/n.dam, n.dam-1)
  q = cumsum(q)
  quantiles = quantile(x[,1], q)
  dam = as.numeric( factor(cut(x[,1], breaks = c(-Inf, quantiles, Inf))) )   
  quantiles = quantile(x[,2], 0.45)
  diet = as.numeric( factor(cut(x[,2], breaks = c(-Inf, quantiles, Inf))) )
  diet = diet - 1   
  
  da = rep(NA, n.dam)
  sigma = runif(1, 0.5, 1); v0 = c(v0, sigma^2); sigma^2
  for (j in 1:n.dam) da[j] = rnorm(1, 0, sigma)
  logT = runif(n, 7.1, 10.5)
  u = -7
  mu0 = logT + u
  theta = runif(1, 0.1, 5); th0 = c(th0, theta); theta
  b = runif(1, b0[1], b0[2]); bb = c(bb, b); b
  ys = sim.y(x = diet, mu = mu0 + da[dam], sigma = 1, coefs = b, p.neg = 0, nb.theta = theta) 
  y0 = ys$y.nb
  N = exp(logT)
  
  tryAgain = TRUE
  infiniteloopcounter = 1
  while (tryAgain & infiniteloopcounter < 5) {
    y1 = asin(sqrt((y0)/(N)))
    if (anyNA(y1) == TRUE) {
      tryAgain = TRUE
      infiniteloopcounter = infiniteloopcounter + 1
    } else {
      tryAgain = FALSE
    }
  }
  if (infiniteloopcounter >= 5) {
    stop("Consistent error found during simulation. Need to investigate cause.")
  }
  
  y1 = y1/sd(y1)
  f1 = lme(y1 ~ diet, random = ~ 1 | dam, verbose = F) 
  out1 = rbind(out1, summary(f1)$tTable["diet", ][c(1,2,5)])
  out1
  
  y2 = log((y0 + 1)/N)
  y2 = y2/sd(y2)
  f2 = lme(y2 ~ diet, random = ~ 1 | dam, verbose = F) 
  out2 = rbind(out2, summary(f2)$tTable["diet", ][c(1,2,5)])
  out2
  
  f3 = glmm.nb(y0 ~ offset(log(N)) + diet, random = ~ 1 | dam, family = "nb", verbose = F) 
  out3 = rbind(out3, summary(f3)$tTable["diet", ][c(1,2,5)])
  out3
  th = c(th, f3$nb.theta); th 
  v = c(v, getVarCov(f3)); v
  vresi = c(vresi, f3$sigma); vresi
    
}
stop.time <- Sys.time()
round(difftime(stop.time, start.time, units = "min"), 3)

pow = cbind(out2[,3], out1[,3], out3[,3])
est = out3[, 1]
est = est - bb
res = sim.out(coefs.p = t(pow), coefs.est = t(est), alpha = c(0.05, 0.01, 0.005, 0.001))
rownames(res[[1]]) = c("LMM_arcsine", "LMM_log", "NBMM")
res
```

## Real Data Analysis
The following R code is used for real data analysis in a manuscript in revision and the citation will added later. The dataset we analyzed was published in Leamy, L.J., et al. Host genetics and diet, but not immunoglobulin A expression, converge to shape compositional features of the gut microbiome in an advanced intercross population of mice. Genome Biol 2014;15(12):552. 

```r
library(NBZIMM)
setwd()

pheno = read.csv("F10_Species_All.csv", header = T) ##or genus family order class phylum 

dim(pheno)
colnames(pheno)
rownames(pheno) = pheno[, 1]

yy = as.matrix(pheno[, 2:(ncol(pheno)-6)])
yy = ifelse(is.na(yy), 0, yy)
zero.p = apply(yy, 2, function(x) {length(x[x != 0])/length(x)} )
zero.p = sort(zero.p, decreasing = T)
zero.p = data.frame(zero.p)
zero.p$id = rownames(zero.p)
zero.p = data.frame(zero.p[zero.p$zero.p>0.1, ])
yy = yy[, rownames(zero.p)]

diet = pheno[, "Diet"]
table(diet)
dam = pheno[, "DamID"]
table(dam); length(table(dam))
sire = pheno[, "SireID"]
table(sire); length(table(sire))
parity = pheno[, "Parity"]
table(parity)

N = pheno[, "TotalReads"]  # total reads
mean(N); sd(N)
mean(log(N)); sd(log(N))

##########################
# NBMM model

di = theta = p.adjust = NULL
for (j in 1:ncol(yy)){
  y = yy[, j]
  f = glmm.nb(y ~ diet + offset(log(N)), random = ~ 1 | dam)
  res = summary(f)
  di = c(di, res$tTable[nrow(res$tTable), 5])
  theta = c(theta, res$nb.theta)
}
p.adjust <- as.vector(as.numeric(p.adjust(di, method = "fdr", n = length(di))))
res1 = cbind(di, p.adjust)
rownames(res1) = colnames(yy)
result1 <- data.frame(res1)
result1$Method <- "NBMM"
result1$names <- rownames(result1)
result1 <- result1[result1$p.adjust < 0.05, ]

# LMM arcsine model
di = p.adjust = NULL
for (j in 1:ncol(yy)){
  y = yy[, j]
  y0 = asin(sqrt(y/N))
  f = lme(y0 ~ diet, random = ~ 1 | dam)
  res = summary(f)
  di = c(di, res$tTable[nrow(res$tTable), 5])
}
p.adjust <- as.vector(as.numeric(p.adjust(di, method = "fdr", n = length(di))))
res2 = cbind(di, p.adjust)
rownames(res2) = colnames(yy)
result2 <- data.frame(res2)
result2$Method <- "LMM_arcsine"
result2$names <- rownames(result2)
result2 <- result2[result2$p.adjust < 0.05, ]

# LMM log model
di = p.adjust = NULL
for (j in 1:ncol(yy)){
  y = yy[, j]
  y0 = log((y + 1)/N)
  f = lme(y0 ~ diet, random = ~ 1 | dam)
  res = summary(f)
  di = c(di, res$tTable[nrow(res$tTable), 5])
}
p.adjust <- as.vector(as.numeric(p.adjust(di, method = "fdr", n = length(di))))
res3 = cbind(di, p.adjust)
rownames(res3) = colnames(yy)
result3 <- data.frame(res3)
result3$Method <- "LMM_log"
result3$names <- rownames(result3)
result3 <- result3[result3$p.adjust < 0.05, ]

result <- rbind(result1, result2)
result$log10p <- -log10(result$p.adjust)
#colnames(result)[8] <- "-log10(p)"

### generate the plot for manuscript ###
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(result$log10p, result$Method, function(x) max(x))
x = sort(x, TRUE)
result$Method = factor(as.character(result$Method), levels=names(x))
x = tapply(result$log10p, result$names, function(x) max(x))
x = sort(x, TRUE)
result$Species = factor(as.character(result$names), levels=names(x))
result$level = "Species"
p <- ggplot(result, aes(x=log10p, y=Species, color=Method)) + 
     geom_point(size=4) + theme(legend.position = "none") + 
     theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5)) +
     xlab(expression('-log'[10]*'(P)')) + theme(text = element_text(size=20)) +
     facet_grid(. ~ level) + coord_fixed()
p
```

## Contact
Feel free to contact me by nyi AT uab.edu
