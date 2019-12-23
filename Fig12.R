###################################################################################
#                                                                                 #
#   Filename    :       Fig12.R                          						  #						                                                                                     
#   Project     :       BiomJ article "A Utility Approach to Individualized       #
#                       Optimal Dose Selection Using Biomarkers    				  #
#   Authors     :       Pin Li                                                    #    
#   Date        :       1.13.2019                                                 #
#   Required R packages :  tmvtnorm, stats, ggplot2                               #
###################################################################################

rm(list = ls())
library(tmvtnorm)
library(stats)
library(ggplot2)

#define the invlogit function
invlogit <- function(x){
	return(exp(x) / (exp(x) + 1))
}

#sampling with sample size n and number of covariates p
my.sample <- function(n, p){
	x <- rmvnorm(n * 2, mean = rep(0, p))
	x <- x[abs(0.4 * rowSums(x[, 1:3]) - 0.8 * x[, 4]) < 1, ]
	x <- x[1:n, ]

	a <- runif(n, -1, 1)

	h_y <- (1 + 0.4 * rowSums(x[, 1:3]) - 0.8 * x[, 4]) * a
	h_r <- (1 - 0.4 * rowSums(x[, 1:3]) + 0.8 * x[, 4]) * a

	y <- rbinom(n, 1, invlogit(x[, 1] + h_y))
	r <- rbinom(n, 1, invlogit(-1.386 - x[, 1] + h_r))

	data.frame(y, r, a, x)
}


#get the matrix with interactions
design <- function(mat){
	int_mat <- model.matrix( ~ . ^ 2, data = data.frame(mat))
	return(int_mat[, 2:(2 * p + 2)])
}


set.seed(1)
n <- 200
p <- 5
case1 <- my.sample(n, p)
case1_val <- case1


#true coefficients
y_coef <- c(0, 1, 1, rep(0, p - 1), rep(0.4, 3), -0.8, rep(0, p - 4))
r_coef <- c(-1.386, 1, -1, rep(0, p - 1), rep(-0.4, 3), 0.8, rep(0, p - 4))

p_y_0 <- invlogit(cbind(1, design(cbind(a = -1, case1_val[, -c(1:3)]))) %*% y_coef)
p_r_0 <- invlogit(cbind(1, design(cbind(a = -1, case1_val[, -c(1:3)]))) %*% r_coef)

#define function to find optimal dose for the population at given lambda
my.opt <- function(lambda) {
    D <- rep(NA, n)
    for (j in 1:n) {
        f <- function(a) invlogit(c(1, design(cbind(a, case1_val[, -c(1:3)]))[j, ]) %*% y_coef)- 
            lambda * invlogit(c(1, design(cbind(a, case1_val[, -c(1:3)]))[j, ]) %*% r_coef) - 
            p_y_0[j] + lambda * p_r_0[j]
        D[j] <- optimize(f, c(-1, 1), tol = 0.001, maximum = TRUE)$maximum
    }
    return(D)
}


head(as.matrix(case1_val[, -c(1:3)]) %*% y_coef[8:12] + y_coef[2])
head(as.matrix(case1_val[, -c(1:3)]) %*% r_coef[8:12] + r_coef[2])

# choose sub 4, 11, 15 for Fig 1
ind <- 11

#get the predicted outcome at given dose
a <- seq(-1, 1, 0.01)
y_1 <- rep(NA, 201)
r_1 <- rep(NA, 201)
for (ai in 1:201) {
    y_1[ai] <- invlogit(c(1, design(cbind(a[ai], case1_val[, -c(1, 2, 3)]))[ind, ]) %*% y_coef)
    r_1[ai] <- invlogit(c(1, design(cbind(a[ai], case1_val[, -c(1, 2, 3)]))[ind, ]) %*% r_coef)
}

ggplot(data = data.frame(a, y_1, r_1), aes(x = a)) + 
	geom_line(aes(y = 0.5 + (y_1 - y_1[1]) - (r_1 - r_1[1]), linetype = "a")) + 
    geom_line(aes(y = y_1 - y_1[1], linetype = "b")) + 
    geom_line(aes(y = r_1 - r_1[1], linetype = "c")) + 
    scale_linetype_manual(values = c(1, 2, 3), name = "", labels = expression(Utility, delta[E], delta[T])) +
    xlab("dose") + 
    scale_y_continuous(limits = c(0, 1), "Predicted outcome", sec.axis = sec_axis(~. - 0.5, name = "Utility")) +
    theme_classic() + 
    theme(legend.position = "top")
    
    
#get the optimal dose at given lambda
lambda <- seq(0.1, 4, length = 100)
d_1 <- rep(NA, 100)
for (li in 1:100) {
    f <- function(a) invlogit(c(1, design(cbind(a, case1_val[, -c(1, 2, 3)]))[ind, ]) %*% y_coef) - lambda[li] * invlogit(c(1, design(cbind(a, case1_val[, -c(1, 2, 3)]))[ind, ]) %*% r_coef) - p_y_0[ind] + lambda[li] * p_r_0[ind]
    d_1[li] <- optimize(f, c(-1, 1), tol = 0.001, maximum = TRUE)$maximum
}


ggplot(data.frame(lambda, d_1), aes(y = d_1, x = lambda)) + 
    ylim(-1, 1) + 
    geom_line() + 
    theme_classic() + 
    xlab(expression(theta)) + 
    ylab("Optimal dose")


#get the population average outcome at given lambda for Fig 2
y_d <- rep(NA, 100)
r_d <- rep(NA, 100)
for (li in 1:100) {
    D_d <- my.opt(lambda[li])
    r_d[li] <- mean(invlogit(cbind(1, design(cbind(a = D_d, case1_val[, -c(1, 2, 3)]))) %*% r_coef))
    y_d[li] <- mean(invlogit(cbind(1, design(cbind(a = D_d, case1_val[, -c(1, 2, 3)]))) %*% y_coef))
}
 

ggplot(data = data.frame(lambda, y_d, r_d), aes(x = lambda)) + 
    geom_line(aes(y = y_d, linetype = "a")) + 
    geom_line(aes(y = r_d, linetype = "b")) + 
    ylim(0, 1) + 
    xlab(expression(theta)) + 
    ylab("Average outcome") + 
    theme_classic() + 
    theme(legend.position = c(0.25, 0.85)) + 
    scale_linetype_manual(values = c(2, 3), name = "", 
                          labels = expression(E[x](Pr(paste(italic(E), "|", italic(d^{opt})),italic(x))),
                                              E[x](Pr(paste(italic(T), "|", italic(d^{opt})),italic(x)))))

ggplot(data.frame(y_d, r_d), aes(y = y_d, x = r_d)) + 
    geom_line() + 
    theme_classic() + 
    xlab("Average toxicity") + 
    ylab("Average efficacy")                              


