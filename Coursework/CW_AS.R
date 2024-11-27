
S_1 = solve(opt_res_q2_1$hessian)

nll_gr_q2_1 = nll_pack_q2_1$gr(opt_res_q2_1$par, dataQ2, apply.sum = F)

Khat_1 = t(nll_gr_q2_1) %*% nll_gr_q2_1

cov_1 = S_1 %*% Khat_1 %*% S_1
typeof(cov_1)
cov_11 <- matrix(nrow=2, ncol=2)
## To obtain the variance of theta square we use the fact that the square of
## a normal distribution is non central chi squared with mean=(mu^2 + variance) 
## and variance = 2*(variance^2). covariance between theta and theta^2 
## is 2*mu*var
cov_11 <- matrix(data=c(cov_1, 2*opt_res_q2_1$par*cov_1, 2*opt_res_q2_1$par*cov_1, 2*(cov_1^2)), nrow=2, ncol=2)

for (i in 1:n_grid){
  
  vec.x <- c(1, x[i])
  
  est_1 <- crossprod(vec.x, c(opt_res_q2_1$par, -(opt_res_q2_1$par^2)-(cov_1^2)))
  
  se_1 <- sqrt(crossprod(vec.x, cov_11) %*% vec.x)
  
  P23_bands_1[i,] <- 1/c(est_1 + 1.96 * se_1,
                         
                         est_1,
                         
                         est_1 - 1.96 * se_1)
  
}
plot(y~x, dataQ2)
lines(x, P23_bands_1[,1], col="red")
lines(x, P23_bands_1[,2], col="green")
lines(x, P23_bands_1[,3], col="blue")

S_2 = solve(opt_res_q2_2$hessian)

nll_gr_q2_2 = nll_pack_q2_2$gr(opt_res_q2_2$par, dataQ2, apply.sum = F)

Khat_2 = t(nll_gr_q2_2) %*% nll_gr_q2_2

cov_2 = S_2 %*% Khat_2 %*% S_2

for (i in 1:n_grid){
  
  vec.x <- c(1, x[i])
  
  est_2 <- crossprod(vec.x, opt_res_q2_2$par)
  
  se_2 <- sqrt(crossprod(vec.x, cov_2) %*% vec.x)
  
  P23_bands_2[i,] <- 1/c(est_2 + 1.96 * se_2,
                         
                         est_2,
                         
                         est_2 - 1.96 * se_2)
  
}

plot(y~x, dataQ2)
lines(x, P23_bands_2[,1], col="red")
lines(x, P23_bands_2[,2], col="green")
lines(x, P23_bands_2[,3], col="blue")

S_3 = solve(opt_res_q2_3$hessian)

nll_gr_q2_3 = nll_pack_q2_3$gr(opt_res_q2_3$par, dataQ2, apply.sum = F)

Khat_3 = t(nll_gr_q2_1) %*% nll_gr_q2_1

cov_3 = S_1 %*% Khat_1 %*% S_1

cov_33 <- matrix(nrow=2, ncol=2)
cov_33 <- matrix(data=c(cov_1[1:1], 2*opt_res_q2_1$par[1]*cov_1[1:1], 2*opt_res_q2_1$par[1]*cov_1[1:1], 2*(cov_1[1:1]^2)), nrow=2, ncol=2)

for (i in 1:n_grid){
  
  vec.x <- c(1, x[i])
  
  est_1 <- crossprod(vec.x, c(opt_res_q2_3$par[1], -(opt_res_q2_3$par[1]^2)))
  
  se_1 <- sqrt(crossprod(vec.x, cov_33) %*% vec.x)
  
  P23_bands_3[i,] <- 1/c(est_1 + 1.96 * se_1,
                         
                         est_1,
                         
                         est_1 - 1.96 * se_1)
  
}
plot(y~x, dataQ2)
lines(x, P23_bands_3[,1], col="red")
lines(x, P23_bands_3[,2], col="green")
lines(x, P23_bands_3[,3], col="blue")


n_grid <- 100
x <- seq(0,1,length=n_grid)
mean_KL1   <-rep(NA,n_grid)
mean_KLfull <- rep(NA, n_grid)
## Plot of full model for comparison
for (i in 1:n_grid){
  vec.x <- c(1, x[i])
  est_full <- crossprod(vec.x,c(2/3, -1/3))
  mean_KLfull[i] <- 1/est_full
}
ci_full <- matrix(NA,
                  nrow = n_grid,
                  ncol = 2)
for (i in 1:n_grid){
  vec.x <- c(1, x[i])
  est_full <- crossprod(vec.x,c(opt_res_q2_kl1$par, -opt_res_q2_kl1$par^2))
  mean_KL1[i] <- 1/est_full
}
plot(y~x, dataQ2, ylim=c(0,5))
lines(x, mean_KL1, col="red")
lines(x, mean_KLfull, col="green")

KL_pack_2 = KL_pack(nll_pack_q2_2)
opt_res_q2_kl2 = fit_optim(2, KL_pack_2, dataQ2, mean = 0.5, sd = 0.25, silent = F)
n_grid <- 100
x <- seq(0,1,length=n_grid)
mean_KL2   <-rep(NA,n_grid)
for (i in 1:n_grid){
  vec.x <- c(1, x[i])
  est_full <- crossprod(vec.x,opt_res_q2_kl2$par)
  mean_KL2[i] <- 1/est_full
}
plot(y~x, dataQ2, ylim=c(0,5))
plot(y~x, dataQ2)
lines(x, mean_KL2, col="red")
lines(x, mean_KLfull, col="green")

KL_pack_3 = KL_pack(nll_pack_q2_3)
opt_res_q2_kl3 = fit_optim(2, KL_pack_3, dataQ2, mean = 0.5, sd = 0.25, silent = F)
n_grid <- 100
x <- seq(0,1,length=n_grid)
mean_KL3   <-rep(NA,n_grid)
for (i in 1:n_grid){
  vec.x <- c(1, x[i])
  est_full <- crossprod(vec.x,c(opt_res_q2_kl3$par[1], -(opt_res_q2_kl3$par[1]^2)))
  mean_KL3[i] <- 1/est_full
}
plot(y~x, dataQ2, ylim=c(0,5))
plot(y~x, dataQ2)
lines(x, mean_KL3, col="red")
lines(x, mean_KLfull, col="green")

## Based on plots, the best worst mean function are from the first and third 
## model depending on use. The first model works equally well across all values
## x while the thrid model performs better at lower values of x but not as well
## for larger values.