nll_expr_q2_full = expr(lgamma(exp(-theta)) - exp(-theta) * (log(alpha + beta * x) - theta) + (1 - exp(-theta)) * log(y) + exp(-theta) * (alpha + beta * x) * y)
nll_pack_q2_full = deriv_pack(nll_expr_q2_full)
opt_res_q2_full = fit_optim(3, nll_pack_q2_full, dataQ2)
se <- sqrt(solve(opt_res_q2_full$hessian)*((nll_pack_q2_full$gr(opt_res_q2_full$par, dataQ2)^2)*solve(opt_res_q2_full$hessian)))
n_grid <- 100
x <- seq(0,1,length=n_grid)
mean_full   <-rep(NA,n_grid)
ci_full <- matrix(NA,
                  nrow = n_grid,
                  ncol = 2)
for (i in 1:n_grid){
  vec.x <- c(1, x[i])
  est_full <- crossprod(vec.x,opt_res_q2_full$par[1:2])
  se_full <- sqrt(crossprod(vec.x, se[1:2, 1:2]) %*% vec.x)
  mean_full[i] <- 1/est_full
  ci_full[i,] <- c(1/(est_full - 1.96 * se_full),1/(est_full + 1.96 * se_full)) 
}
plot(y~x, dataQ2, ylim=c(0,5))
lines(x, mean_full, col="red")
lines(x, ci_full[,1], col="green")
lines(x, ci_full[,2], col="blue")