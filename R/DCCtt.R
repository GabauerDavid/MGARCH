#' @title DCC-GARCH of Tse and Tsui (2002)
#' @description This function estimates the DCC-GARCH model of Tse and Tsui (2002)
#' @param residuals Residuals of the univariate GARCH models
#' @param sigma Conditional volatility of univariate GARCH models
#' @param window.size Window size. Default is 1.
#' @return returns DCC-GARCH results
#' @importFrom zoo index
#' @importFrom stats cor cov optim pt
#' @references
#' Tse, Y. K., & Tsui, A. K. C. (2002). A multivariate generalized autoregressive conditional heteroscedasticity model with time-varying correlations. Journal of Business and Economic Statistics, 20(3), 351-362.
#' @author David Gabauer
#' @export
DCCtt = function(residuals, sigma, window.size=1) {
  
  std_residuals = residuals / sigma
  
  dcc_tt_loglik = function(params, std_residuals, window.size=1) {
    std_residuals = as.matrix(std_residuals)
    alpha = params[1]
    beta = params[2]
    
    if (alpha < 0 || beta < 0 || (alpha + beta) > 1) return(1e6)
    
    n = nrow(std_residuals)
    p = ncol(std_residuals)
    
    R = array(0, dim = c(p, p, n))
    R_bar = cor(std_residuals)
    R[,,1:window.size] = R_bar
    
    log_likelihood = 0
    for (t in 2:n) {
      if (t<=window.size) {
        Rs = R_bar
      } else {
        E = std_residuals[(t-window.size):(t-1),,drop=F]
        Q = t(E) %*% E
        B = diag(1/sqrt(diag(Q)))
        Rs = B %*% Q %*% B
      }
      R[,,t] = (1-alpha-beta)*R_bar + alpha*Rs + beta*R[,,t-1]
      det_R = det(R[,,t])
      if (det_R <= 0) return(1e6)
      if (t>window.size) {
        log_likelihood = log_likelihood + log(det_R) + c(t(std_residuals[t,]) %*% solve(R[,,t]) %*% t(t(std_residuals[t,])))
      }
    }
    return(log_likelihood)
  }
  
  params = c(0.005, 0.99)
  opt_result = optim(par=params, 
                     fn=dcc_tt_loglik, 
                     std_residuals=std_residuals,
                     window.size=window.size,
                     method="L-BFGS-B", 
                     lower=c(0,0), 
                     upper=c(1,1),
                     hessian=TRUE)
  
  hessian_matrix = opt_result$hessian
  standard_errors = sqrt(diag(solve(hessian_matrix)))
  param_estimates = opt_result$par
  log.likelihoods = opt_result$value
  tstat = param_estimates / standard_errors
  pval = 2 * (1 - pt(tstat, df=nrow(residuals)))
  matcoef = cbind(opt_result$par, standard_errors, tstat, pval)
  rownames(matcoef) = c("dcca1", "dccb1")
  colnames(matcoef) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  alpha = opt_result$par[1]
  beta  = opt_result$par[2]
  NAMES = colnames(std_residuals)
  dates = as.character(index(std_residuals))
  
  # Get correlation and covariance matrix
  T = nrow(std_residuals)
  k = ncol(std_residuals)
  H = R = array(NA, dim = c(k, k, T), dimnames=list(NAMES, NAMES, dates))
  R_bar = cor(std_residuals)
  R[,,1:window.size] = R_bar
  H[,,1:window.size] = cov(residuals)
  for (t in 2:T) {
    if (t<=window.size) {
      R[,,t] = R_bar
    } else {
      E = std_residuals[(t-window.size):(t-1),,drop=F]
      Q = t(E) %*% E
      B = diag(1/sqrt(diag(Q)))
      Rs = B %*% Q %*% B
      R[,,t] = (1-alpha-beta)*R_bar + alpha*Rs + beta*R[,,t-1]
    }
    D_inv = diag(as.numeric(sigma[t,]))
    H[,,t] = D_inv%*%R[,,t]%*%D_inv
  }
  
  return(list(H=H, R=R, matcoef=matcoef, std_residuals=std_residuals, log.likelihoods=log.likelihoods))
}
