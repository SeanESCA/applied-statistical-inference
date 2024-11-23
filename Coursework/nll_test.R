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