# load the auxilary functions
source("functions.R")
source("dgp.R")
source("dgp2.R")

rangeval = c(0,10)
out_p = 0.1 # outlier proportion 10%

# data generation
set.seed(123)
dat = dgp(500)
X = dat$X
Y = dat$Y
dat_test = dgp(500)
X_test=dat_test$X
Y_test=dat_test$Y

# outlier generation
dat2=dgp2(500)
X_out = dat2$X

index = which(Y == 1)

out_index=sample(index, out_p * 500, replace = F)

X[out_index,] = X_out[out_index,]
Y[out_index] = 0

# get the desing matrix (H)
H = getHmat(data = X, nbasis = 15, rangeval = rangeval)
Hmat=H$Hmat
Bsf=H$Bsf

H2 = getHmat(data = X_test, nbasis = 15, rangeval = rangeval)
Hmat_test = H2$Hmat

# run the pca, pls, and rpls methods
pca_model = log_pca(Y = Y, X = X, X_test = X_test, nbasis = 15, ncomp = 4, rangeval = rangeval)
pls2_model = log_fpls2(Y = Y, Hmat = Hmat, Hmat_test = Hmat_test, alpha = 0.05, hmax = 10, Bsf=Bsf,  model = "ls")
rpls_model = log_fpls2(Y = Y, Hmat = Hmat, Hmat_test = Hmat_test, alpha = 0.05, hmax = 10, Bsf=Bsf,  model = "wle")

beta_t=dat$Bt # true parameter function
# mise values
mise_pca = mean((c(beta_t$data)-c(pca_model$beta_hat_s))^2)   # 0.5033169
mise_pls2 = mean((c(beta_t$data)-c(pls2_model$beta_hat_s))^2) # 0.5780343 
mise_rpls = mean((c(beta_t$data)-c(rpls_model$beta_hat_s))^2) # 0.1654021

# ccr values
pca_ccr = ccr(Y_test, pca_model$preds)
pls2_ccr = ccr(Y_test, pls2_model$preds)
rpls_ccr = ccr(Y_test, rpls_model$preds)

ccr_pca = pca_ccr$ccr   # 0.656
ccr_pls2 = pls2_ccr$ccr # 0.596
ccr_rpls = rpls_ccr$ccr # 0.788

# auc values
auc_pca = auc(Y_test ~ c(pca_model$preds))[1]  # 0.6458333
auc_pls2 = auc(Y_test ~ pls2_model$preds)[1]   # 0.5905449
auc_rpls = auc(Y_test ~ rpls_model$preds)[1]   # 0.7844551
