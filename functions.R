#___________Packages_________
library(fda)
library(wle)
library(robustX)  
library(matrixStats) 
library(pROC)
#____________________________

options(warn = -1)

#____________________L1_median_function_________________

med_L1 = function(X){
  nc = dim(X)[2]
  nr = dim(X)[1]
  
  medin = apply(X, 2, median, na.rm = TRUE)
  loc = function(X, nc) X - rep(nc, each = nr)
  scl = apply(abs(loc(X, medin)), 2, mean, trim=0.2, na.rm = TRUE)
  ps = scl>0
  
  if(sum(ps)){
    medin[ps] = L1median(X[, ps, drop=FALSE], m.init = medin[ps],
                           pscale=scl[ps])$estimate
  }
  return(medin)
}
#________________________________________________________


#______________________Obtain_H_matrix___________________
getHmat = function(data, nbasis, rangeval){
  n = dim(data)[1]
  p = dim(data)[2]
  dimnames(data)=list(as.character(1:n), as.character(1:p))
  grid_points = seq(rangeval[1], rangeval[2], length.out = p)
  bs_basis = create.bspline.basis(rangeval, nbasis = nbasis)
  evalbase = eval.basis(grid_points, bs_basis)
  innp_mat = inprod(bs_basis, bs_basis)
  fdobj = fdPar(bs_basis, int2Lfd(2), lambda=0)
  pcaobj = smooth.basisPar(grid_points, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  H = t(pcaobj$coefs) %*% innp_mat
  return(list(Hmat = H, Bsf = evalbase))
}
#________________________________________________________


#____________________________Logistic_Functional_PLS_proposed_based_on_Dodge_and_Whittaker (2009)_____________________________
log_fpls2 = function(Y, Hmat, Hmat_test, alpha, hmax, Bsf, model){
  
  Htest = Hmat_test
  
  if(model == "ls"){
    Hmat = scale(Hmat)
    Hmat_test = scale(Hmat_test)
  }else  if(model == "wle"){
    minw = 0.1
    medH = med_L1(Hmat)
    medH_test = med_L1(Hmat_test)
    cH = sapply(1:ncol(Hmat), function(i) (Hmat[,i]-medH[i]))
    cH_test = sapply(1:ncol(Hmat_test), function(i) (Hmat_test[,i]-medH_test[i]))
    Hmat = sapply(1:ncol(cH), function(i) cH[,i]/median(abs(cH[,i]))/0.6745)
    Hmat_test = sapply(1:ncol(cH_test), function(i) cH_test[,i]/median(abs(cH_test[,i]))/0.6745)
    wH = sapply(1:ncol(Hmat), function(i) wle.weights(Hmat[,i], y=NULL, smooth=0.0031, 1, raf=1, 
                                                      location=FALSE, max.iter=1000, tol=10^(-6))$weights)
    wH_test = sapply(1:ncol(Hmat_test), function(i) wle.weights(Hmat_test[,i], y=NULL, smooth=0.0031, 1, raf=1, 
                                                      location=FALSE, max.iter=1000, tol=10^(-6))$weights)
    wH = wH/max(wH)
    wH_test = wH_test/max(wH_test)
    wH[wH >= 1-minw] = 1
    wH_test[wH_test >= 1-minw] = 1
    wH[wH <= minw] = 0
    wH_test[wH_test <= minw] = 0
    wH = rowMedians(wH)
    wH_test = rowMedians(wH_test)
    Hmat = sapply(1:ncol(Hmat), function(i) sqrt(wH)*Hmat[,i])
    Hmat_test = sapply(1:ncol(Hmat_test), function(i) sqrt(wH_test)*Hmat_test[,i])
  }

  V = NULL
  T = NULL
  T_test = NULL
  
  nbf = dim(Hmat)[2]
  zval = abs(qnorm(alpha/2))
  for(j in 1:hmax){
    try({
      delta = numeric()
      
      for(i in 1:nbf){
        if(model == "ls"){
          log_mod = glm(Y~Hmat[,i], family = binomial)
          logCf = log_mod$coefficients[2]
          sdCf = summary(log_mod)$coefficients[2,2]
          Wval = logCf / sdCf
          if(abs(Wval) <= abs(zval)) logCf = 0
        }else if(model == "wle"){
          log_mod = wle.glm(Y~Hmat[,i], family = binomial)
          logCf = log_mod$root1$coefficients[2]
          sdCf = summary(log_mod)$coefficients[2,2]
          Wval = logCf / sdCf
          if(abs(Wval) <= abs(zval)) logCf = 0
        }
        
        delta[i] = logCf
      }
      delta = as.matrix(delta)
      delta = delta / norm(delta, type = "F")
      
      V = cbind(V, delta)
      T = cbind(T, Hmat %*% delta)
      T_test = cbind(T_test, Hmat_test %*% delta)
      
      for(jj in 1:nbf){
        Hmat[,jj] = Hmat[,jj] - mean(Hmat[,jj]) - as.numeric(cov(Hmat[,jj], T[,j]) / var(T[,j])) * (T[,j] - mean(T[,j]))
        Hmat_test[,jj] = Hmat_test[,jj] - mean(Hmat_test[,jj]) - as.numeric(cov(Hmat_test[,jj], T_test[,j]) / var(T_test[,j])) * (T_test[,j] - mean(T_test[,j]))
      }
      
    }, silent = TRUE)
  }
  
  Vmat = V
  Vmat = Vmat[,!apply(Vmat,2,mean) %in% NaN]
  
  Tmat = T
  Tmat = Tmat[,!apply(Tmat,2,mean) %in% NaN]
  
  Tmat_test = T_test
  Tmat_test = Tmat_test[,!apply(Tmat_test,2,mean) %in% NaN]
  
  if(model == "ls"){
    Fin_log_mod = glm(Y~Tmat, family = binomial)
    Gam_par = Fin_log_mod$coefficients[-1]
    fits = round(Fin_log_mod$fitted.values)
    Gam_constant = Fin_log_mod$coefficients[1]
    
  }else if(model == "wle"){
    Fin_log_mod = wle.glm(Y~Tmat, family = binomial)
    Gam_par = Fin_log_mod$root1$coefficients[-1]
    fits = round(Fin_log_mod$root1$fitted.values)
    Gam_constant = Fin_log_mod$root1$coefficients[1]
  }
  
  preds1 = cbind(1,Tmat_test) %*% as.matrix(c(Gam_constant, Gam_par))
  preds = round(exp(preds1)/(1+exp(preds1)))
  
  Gam_par <- as.matrix(Gam_par)
  
  beta_hat = Vmat%*%Gam_par
  beta_hat_s = Bsf%*%beta_hat
  
  return(list(fits = fits, beta_hat_s = beta_hat_s, beta_hat=beta_hat,
              Gam_constant=Gam_constant, nbf=nbf, preds=preds))
  
}
##____________________________________________________________________________________________________________________________


#_______________Correct_classification_rate_____________
ccr = function(Y, Ypred){
  n = length(Y)
  p = 0
  for(i in 1:n){
    if(Y[i] == Ypred[i])
      p = p+1
  }
  ccr = p / n
  ce = 1 - ccr
  return(list(ccr = ccr, ce = ce))
}
#_______________________________________________________


getPCA = function(X_tr, X_te, nbasis, ncomp, rangeval){
  data = rbind(X_tr, X_te)
  n_train = dim(X_tr)[1]
  n_test = dim(X_te)[1]
  n = dim(data)[1]
  p = dim(data)[2]
  dimnames(data)=list(as.character(1:n), as.character(1:p))
  grid_points = seq(rangeval[1], rangeval[2], length.out = p)
  bs_basis = create.bspline.basis(rangeval, nbasis = nbasis)
  evalbase = eval.basis(grid_points, bs_basis)
  fdobj = fdPar(bs_basis, int2Lfd(2), lambda=0)
  pcaobj = smooth.basisPar(grid_points, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  dpca = pca.fd(pcaobj[1:n_train,], nharm = ncomp, fdobj)
  PCAscore = dpca$scores
  PCAcoef = dpca$harmonics$coefs
  fdXnew = pcaobj[-(1:n_train),]
  fdXnew$coefs = t(scale(t(fdXnew$coefs), scale = F))
  PCAscore_test = inprod(fdXnew, dpca$harmonics)
  return(list(PCAcoef = PCAcoef, PCAscore = PCAscore, PCAscore_test = PCAscore_test, evalbase = evalbase))
}

log_pca = function(Y, X, X_test, nbasis, ncomp, rangeval){
  pcaX = getPCA(X_tr = X, X_te = X_test, nbasis = nbasis, ncomp = ncomp, rangeval = rangeval)
  scoX = pcaX$PCAscore
  scoXtest = pcaX$PCAscore_test
  
  log_model = glm(Y~scoX, family = binomial)
  Gam_par = log_model$coefficients[-1]
  fits = round(log_model$fitted.values)
  Gam_constant = log_model$coefficients[1]
  
  beta_hat = pcaX$PCAcoef%*%Gam_par
  beta_hat_s = pcaX$evalbase%*%beta_hat
  
  preds1 = cbind(1, scoXtest) %*% c(Gam_constant, Gam_par) 
  preds = round(exp(preds1) / (1 + exp(preds1)))
  
  return(list(fits = fits, beta_hat_s = beta_hat_s, preds = preds))
  
}
