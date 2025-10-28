# MGARCH

This project was created out of boredom. One day, I wanted to compare the DCC-GARCH model of Engle (2002) with the DCC-GARCH model of Tse and Tsui (2002) but unfortunately there is no implementation of the latter. The only one I have found was part of the MTS R-package, which was not working on my end and even if it would have worked was quite restrictive. Thus, I decided to implement my own version and make it more flexible. 

## Step 1: Install the devtools package

To install a R package, start by installing the devtools package. The best way to do this is from CRAN, by typing:

```r
install.packages("devtools")
```

## Step 2: Install the package of interest from GitHub

Install the package of interest from GitHub using the following code, where you need to remember to list both the author and the name of the package (in GitHub jargon, the package is the repo, which is short for repository). In this example, we are installing the ConnectednessApproach package created by GabauerDavid.

```r
library(devtools)
install_github("GabauerDavid/MGARCH")
```

## Step 3: Example code
```r

rm(list=ls())

library("zoo")
library("MGARCH")
library("rugarch")
library("rmgarch")
library("quantmod")

getSymbols(c("^GSPC", "^GDAXI"), from="2023-01-01", to="2025-01-01")
X = na.omit(diff(log(cbind(GSPC[,4], GDAXI[,4]))))
colnames(X) = c("SP500", "DAX")
k = ncol(X)

# DCC-GARCH (Tse and Tsui, 2002)
ugarch_spec = ugarchspec(mean.model=list(armaOrder=c(0,0)), 
                         variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                         distribution.model="norm")
garch_fits = lapply(1:k, function(i) ugarchfit(spec=ugarch_spec, data=X[,i]))
residuals = sapply(garch_fits, function(fit) residuals(fit))
sigma = sapply(garch_fits, function(fit) sigma(fit))
colnames(residuals) = colnames(sigma) = colnames(X)
sigma = zoo(sigma, index(X))
residuals = zoo(residuals, index(X))

dcc_tt = DCCtt(residuals, sigma, window.size=1)
H_tt = dcc_tt$H
R_tt = dcc_tt$R
dates = as.Date(dimnames(R_tt)[[3]])

# DCC-GARCH (Engle, 2002)
dcc_spec = dccspec(multispec(replicate(k, ugarch_spec)), dccOrder=c(1,1), 
                   distribution="mvnorm")
dcc_fit = dccfit(dcc_spec, data=X)
H_dcc = rcov(dcc_fit)
R_dcc = rcor(dcc_fit)

# Dynamic conditional volatility
par(mfrow = c(1, 1), oma = c(0, 2, 0, 0) + 0.5, mar = c(1, 1, 1, 1) + 0.5, mgp = c(1, 0.4, 0))
plot(dates, H_tt[1,1,], type='l', las=1, xaxs='i', xlab='', ylab='', yaxs='i', lwd=2, tck=-0.005)
title("Dynamic conditional volatility", adj=0)
grid(NA, NULL)
lines(dates, H_dcc[1,1,], col=2)
legend('topleft', c('Tse and Tsui (2002)', 'Engle (2002)'), fill=1:2, bty='n')
box()

# Dynamic conditional covariance
plot(dates, H_tt[1,2,], type='l', las=1, xaxs='i', xlab='', ylab='', yaxs='i', lwd=2, tck=-0.005)
title("Dynamic conditional covariance", adj=0)
grid(NA, NULL)
lines(dates, H_dcc[1,2,], col=2)
legend('topleft', c('Tse and Tsui (2002)', 'Engle (2002)'), fill=1:2, bty='n')
box()

# Dynamic conditional correlation
plot(dates, R_tt[1,2,], type='l', las=1, ylim=c(0,1), xaxs='i', xlab='', ylab='', yaxs='i', lwd=2, tck=-0.005)
title("Dynamic conditional correlation", adj=0)
grid(NA, NULL)
lines(dates, R_dcc[1,2,], col=2)
legend('topleft', c('Tse and Tsui (2002)', 'Engle (2002)'), fill=1:2, bty='n')
box()
```


## BibTeX Citation

If you use this package in a scientific publication, I would appreciate if you use the following citation:

```
@article{tse2002multivariate,
  title={A multivariate generalized autoregressive conditional heteroscedasticity model with time-varying correlations},
  author={Tse, Yiu Kuen and Tsui, Albert K C},
  journal={Journal of Business and Economic Statistics},
  volume={20},
  number={3},
  pages={351--362},
  year={2002},
  publisher={Taylor \& Francis}
}

@article{gabauer2025,
  title={Package ‘MGARCH’},
  author={Gabauer, David},
  year={2025}
}
```
