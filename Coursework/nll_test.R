## Question 1 - Test
alpha = rnorm(1)
beta = rnorm(1)
theta = rnorm(1)
x = rnorm(1)
y = rnorm(1)

dnorm(y, exp(alpha + beta * x), exp(theta + alpha + beta * x), log = T)
nll_pack_q1_1$fn(c(alpha, beta, theta), list("x" = x, "y" = y))

dnorm(log(y), alpha + beta * x, exp(theta), log = T)
nll_pack_q1_2$fn(c(alpha, beta, theta), list("x" = x, "y" = y))

dlogis(log(y), alpha + beta * x, exp(theta), log = T)
nll_pack_q1_3$fn(c(alpha, beta, theta), list("x" = x, "y" = y))

# Question 2.1
# nll_expr_q2_full = expr(lgamma(k) - k * log(lambda) + (1 - k) * log(y) + lambda * y)
# Reparameterise phi as exp(theta)

nll_expr_q2_full = expr(lgamma(exp(-theta)) - exp(-theta) * (log(alpha + beta * x) - theta) + (1 - exp(-theta)) * log(y) + exp(-theta) * (alpha + beta * x) * y)
nll_pack_q2_full = deriv_pack(nll_expr_q2_full)

nll_expr_q2_1 = expr(lgamma(exp(-theta1)) - exp(-theta1) * (log(theta1 - x * theta1 ^ 2) - theta1) + (1 - exp(-theta1)) * log(y) + exp(-theta1) * (theta1 - x * theta1 ^ 2) * y)
nll_pack_q2_1 = deriv_pack(nll_expr_q2_1,
                           c("theta1"),
                           theta_init = 2)

nll_expr_q2_2 = expr(lgamma(exp(-theta2)) - exp(-theta2) * (log(theta2 + tau2 * x) - theta2) + (1 - exp(-theta2)) * log(y) + exp(-theta2) * (theta2 + tau2 * x) * y)
nll_pack_q2_2 = deriv_pack(nll_expr_q2_2,
                           c("theta2", "tau2"),
                           theta_init = rep(1, 2))

nll_expr_q2_3 = expr(lgamma(exp(-theta4)) - exp(-theta4) * (log(theta3 - x * theta3 ^ 2) - theta4) + (1 - exp(-theta4)) * log(y) + exp(-theta4) * (theta3 - x * theta3 ^ 2) * y)
nll_pack_q2_3 = deriv_pack(nll_expr_q2_3,
                           c("theta3", "theta4"),
                           theta_init = rep(2, 2))

# Test
-dgamma(y, exp(-theta), (alpha + beta * x) * exp(-theta), log=T)
nll_pack_q2_full$fn(c(alpha, beta, theta), list("x" = x, "y" = y))

-dgamma(y, exp(-theta), (theta - x * theta ^ 2) * exp(-theta), log = T)
nll_pack_q2_1$fn(theta, list("x" = x, "y" = y))

-dgamma(y, exp(-theta), (theta + beta * x) * exp(-theta), log = T)
nll_pack_q2_2$fn(c(theta, beta), list("x" = x, "y" = y))

-dgamma(y, exp(-theta), (alpha - 1 * alpha ^ 2) * exp(-theta), log = T)
nll_pack_q2_3$fn(c(alpha, theta), list("x" = 1, "y" = y))

# Q2.4 (Test using exercise from Lab 5)
mu_expr = expr(theta1)
sigma_expr = expr(exp(theta2))
delta_expr = expr(theta3)
tau_expr = expr(exp(theta4))
expr1 = expr((y - !!mu_expr)/!!sigma_expr)
asinh_expr = expr(log(!!expr1 + (1 + (!!expr1) ^ 2) ^ 0.5))
S_expr = expr(sinh(!!tau_expr * !!asinh_expr - !!delta_expr))
C_expr = expr((1 + (!!S_expr) ^ 2) ^ 0.5)
nll_test_expr = expr(-log(!!tau_expr) - log(!!C_expr) + 0.5 * (!!S_expr) ^ 2 + log(!!sigma_expr) + 0.5 * log(2 * pi) + 0.5 * log(1 + (!!expr1) ^ 2))
nll_test_pack = deriv_pack(nll_test_expr,
                           c("theta1", "theta2", "theta3", "theta4"),
                           c("y"),
                           rep(2, 4))
opt_res_test = fit_optim(4, nll_test_pack, list(y = y_sample_q1), silent = F)

n <- 500
y_sample_q6 <- rlogis(n,location = 0, scale= 1)

KL_test_pack = function(nll_pack) {
  fn_integrand = function(y, theta) {
    (nll_pack$fn(theta, list("y" = y), apply.sum = F) +
      dlogis(y, log = T)) * dlogis(y)
  }
  gr_integrand = function(y, theta, col_int = 1) {
    aux = nll_pack$gr(theta, list("y" = y), apply.sum = F) *
      dlogis(y)
    aux[, col_int]
  }
  KL_fn = function(theta) {
    integrate(fn_integrand, -50, 50, theta = theta)$value
  }
  
  KL_gr = function(theta) {
    res = 1:length(theta)
    for(j in res) {
      res[j] = integrate(gr_integrand, -50, 50,
                              theta = theta, col_int = j)$value
    }
    res
  }
  
  list(fn = KL_fn, gr = KL_gr)
}

KL_pack_5 = KL_test_pack(nll_test_pack)
optim(c(0, 0, 0, 0),
      KL_pack_5$fn,
      hessian = T,
      control=list(trace=10,REPORT=1))

test = function(y, theta) {
  nll_test_pack$fn(theta, list("y" = y), apply.sum = F)
}
integrate(test, -10, 10, theta = c(1, 1, 1, 1))

optim_q6_0<-
  optim(par=c(0,0,0,0),
      fn=KL_pack_5$fn,
      hessian = T,
      control=list(trace=10,REPORT=1))
optim_q6<-fit_optim(4, KL_pack_5, NULL, silent = F)

xx <- seq(-10,10,length=100)
yy_dens <- rep(NA,100)
yy_logis <- dens_logis(xx) 

for (i in 1:100){
  yy_dens[i] <- exp(-nll_test_pack$fn(optim_q6$par, list(y = xx[i])))
}


plot(xx,yy_dens,type="l")
lines(xx,yy_logis,lty=2)  

# Q2.5
unique(dataQ2$x)
mean_invQ2 = aggregate(list(y = dataQ2$y), list(x = dataQ2$x), function(x) 1 / mean(x))
fit = lm(y ~ x, mean_invQ2)
summary(fit)

# Q2.3

S_1 = solve(opt_res_q2_1$hessian)
nll_gr_q2_1 = nll_pack_q2_1$gr(opt_res_q2_1$par, dataQ2, apply.sum = F)
Khat_1 = crossprod(nll_gr_q2_1)
cov_1 = S_1 %*% Khat_1 %*% S_1

theta1 = opt_res_q2_1$par
theta.hat = c(theta1 = theta1)
se_theta1 = sqrt(solve(opt_res_q2_1$hessian))
var_theta1 = solve(opt_res_q2_1$hessian)
se_S = sqrt(deltavar(1/(theta1 - x * theta1 ^ 2),
                     meanval = theta.hat,
                     Sigma = cov_1))
test = P23_bands_1
est_1 = 1/(theta1 - x * theta1 ^ 2)
test[, 1] = est_1 - 1.96 * se_S
test[, 2] = est_1
test[, 3] = est_1 + 1.96 * se_S
plot(dataQ2, pch = 19)
lines(x, test[, 1])
lines(x, test[, 2])
lines(x, test[, 3])
plot_bands(test, dataQ2, main = "")
