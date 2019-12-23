###################################################################################
#                                                                                 #
#   Filename    :       Simulation.R                          						  #						                                                                                     
#   Project     :       BiomJ article "A Utility Approach to Individualized       #
#                       Optimal Dose Selection Using Biomarkers "   				  #
#   Authors     :       Pin Li                                                    #    
#   Date        :       8.13.2019                                                 #
#   Required R packages :  tmvtnorm, stats, glmnet, BB, numDeriv                  #
###################################################################################


rm(list = ls())
library(tmvtnorm)
library(stats)
library(glmnet)
library(BB)
library(numDeriv)

###### define the invlogit function
invlogit <- function(x){
	return(exp(x) / (exp(x) + 1))
}

###### sampling with sample size n and number of covariates p
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


###### get the matrix with interactions
design <- function(mat){
	int_mat <- model.matrix( ~ . ^ 2, data = data.frame(mat))
	return(int_mat[, 2:(2 * p + 2)])
}


###### function to find optimal dose
opt.dose <- function(tau, r_coef, y_coef, N, case1_val){
	#true coefficients
	y_coef_t <- c(0, 1, 1, rep(0, p - 1), rep(0.4, 3), -0.8, rep(0, p - 4))
	r_coef_t <- c(-1.386, 1, -1, rep(0, p - 1), rep(-0.4, 3), 0.8, rep(0, p - 4))

	p_y_0 <- invlogit(cbind(1, design(cbind(a = -1, case1_val[, -c(1:3)]))) %*% y_coef)
	p_r_0 <- invlogit(cbind(1, design(cbind(a = -1, case1_val[, -c(1:3)]))) %*% r_coef)

    if(mean(p_r_0) > tau){
        return(list(lambda_h = NA, opt_d = NA, sd_d = NA, pred_y = NA, pred_r = NA))
    } else if (mean(invlogit(cbind(1, design(cbind(a = 1, case1_val[, -c(1:3)]))) %*% r_coef)) < tau){
        return(list(lambda_h = 0, opt_d = 1, sd_d = 0, 
                    pred_y = mean(invlogit(cbind(1, design(cbind(a = 1, case1_val[, -c(1:3)]))) %*% y_coef)),
                    pred_r = mean(invlogit(cbind(1, design(cbind(a = 1, case1_val[, -c(1:3)]))) %*% r_coef))))
    } else {
    	#define function to find optimal dose for the population at given lambda
        my.opt <- function(lambda){
            D <- rep(NA, N)
            for (j in 1:N){
                f <- function(a) invlogit(c(1, design(cbind(a, case1_val[, -c(1:3)]))[j, ]) %*% y_coef) -
                lambda * invlogit(c(1, design(cbind(a, case1_val[, -c(1:3)]))[j, ]) %*% r_coef) -
                p_y_0[j] + lambda * p_r_0[j]
                D[j] <- optimize(f, c(-1, 1), tol = 0.001, maximum = TRUE)$maximum
            }
            return(D)
        } 

    	#gloden search to find best lambda
        lambda_l <- 0
        lambda_h <- 100
        lambda_d <- lambda_l + (sqrt(5) - 1) / 2 * (lambda_h - lambda_l)
    
    
        D_d <- my.opt(lambda_d)
        r_d <- mean(invlogit(cbind(1, design(cbind(a = D_d, case1_val[, -c(1:3)]))) %*% r_coef_t))
        
        while (abs(r_d - tau) > 0.003){
            if (r_d > tau){
                lambda_l <- lambda_d
            } else {
                lambda_h <- lambda_d
            }
            lambda_d <- lambda_l + (sqrt(5) - 1) / 2 * (lambda_h - lambda_l)
            D_d <- my.opt(lambda_d)
            r_d <- mean(invlogit(cbind(1, design(cbind(a = D_d, case1_val[, -c(1:3)]))) %*% r_coef_t))
            #print(c(lambda_d, r_d))
        }
        
        y_d <- mean(invlogit(cbind(1, design(cbind(a = D_d, case1_val[, -c(1:3)]))) %*% y_coef_t))
        
        return(list(lambda_h = lambda_d, opt_d = mean(D_d), sd_d = sd(D_d), pred_y = y_d, pred_r = r_d))
    }
}

###### for illustration, run 10 replications
it_n <- 10

###### scenario 0  as example
p <- 5

theory <- matrix(NA, it_n, 5)

ignore <- matrix(NA, it_n, 2)

min <- matrix(NA, it_n, 2)

mid <- matrix(NA, it_n, 2)

fix <- matrix(NA, it_n, 5)

fs <- matrix(NA, it_n, 5)
fs_s <- matrix(NA, it_n, 4 * p + 4)

lasso <- matrix(NA, it_n, 5)
lasso_s <- matrix(NA, it_n, 4 * p + 4)

classo <- matrix(NA, it_n, 5)
classo_s <- matrix(NA, it_n, 4 * p + 6)

for (it in 1:it_n){
    n <- 200
    N <- n
    
    set.seed(it)
    case1 <- my.sample(n, p)
    case1_val <- case1
    
    case1_X <- data.frame(design(case1[, -c(1:2)]))
    
    ###theoretic
    y_coef <- c(0, 1, 1, rep(0, p - 1), rep(0.4, 3), -0.8, rep(0, p - 4))
    r_coef <- c(-1.386, 1, -1, rep(0, p - 1), rep(-0.4, 3), 0.8, rep(0, p - 4))
    
    ###if ignore risk
    ignore[it, 1] <- mean(invlogit(cbind(1, design(cbind(a = 1, case1_val[, -c(1:3)]))) %*% y_coef))
    ignore[it, 2] <- mean(invlogit(cbind(1, design(cbind(a = 1, case1_val[, -c(1:3)]))) %*% r_coef))
    
    ###if choose min dose
    min[it, 1] <- mean(invlogit(cbind(1, design(cbind(a = -1, case1_val[, -c(1:3)]))) %*% y_coef))
    min[it, 2] <- mean(invlogit(cbind(1, design(cbind(a = -1, case1_val[, -c(1:3)]))) %*% r_coef))
    
    ###if choose mid dose
    mid[it, 1] <- mean(invlogit(cbind(1, design(cbind(a = 0, case1_val[, -c(1:3)]))) %*% y_coef))
    mid[it, 2] <- mean(invlogit(cbind(1, design(cbind(a = 0, case1_val[, -c(1:3)]))) %*% r_coef))
    
    ###theoretic optimal dose
    my_opt <- opt.dose(0.2, r_coef, y_coef, N, case1_val)
    theory[it, 1] <- my_opt$opt_d
    theory[it, 2] <- my_opt$sd_d
    theory[it, 3] <- my_opt$pred_y
    theory[it, 4] <- my_opt$pred_r
    theory[it, 5] <- my_opt$lambda_h
    
    ###Fit regression model for Y given (A, X, AX), R given (A, X, AX)
    md_y <- glm(case1$y ~ . , data = case1_X, family = 'binomial')
     
    md_r <- glm(case1$r ~ . , data = case1_X, family = 'binomial')
    
    ###forward selection, force DOSE in
    null_y <- glm(case1$y ~ a, data = case1_X, family = 'binomial')
    y_forwards <- step(null_y, scope = list(lower = null_y, upper = md_y), direction = 'forward', trace = 0) 
    y_coef <- rep(0, 2 * p + 2)
    y_coef[match(names(y_forwards$coefficients), names(md_y$coeff))] <- y_forwards$coefficients
     
    null_r <- glm(case1$r ~ a, data = case1_X, family = 'binomial')
    r_forwards <- step(null_r, scope = list(lower = null_r, upper = md_r), direction = 'forward', trace = 0) 
    r_coef <- rep(0, 2 * p + 2)
    r_coef[match(names(r_forwards$coefficients), names(md_r$coeff))] <- r_forwards$coefficients
     
    fs_s[it, ] <- c(y_coef, r_coef)
    my_opt <- opt.dose(0.2, r_coef, y_coef, N, case1_val)
    fs[it, 1] <- my_opt$opt_d
    fs[it, 2] <- my_opt$sd_d
    fs[it, 3] <- my_opt$pred_y
    fs[it, 4] <- my_opt$pred_r
    fs[it, 5] <- my_opt$lambda_h
    
     
    
    ###fixed model with only dose
    y_coef <- c(null_y$coefficients, rep(0, 2 * p))
    
    r_coef <- c(null_r$coefficients, rep(0, 2 * p))
    
    my_opt <- opt.dose(0.2, r_coef, y_coef, N, case1_val)
    fix[it, 1] <- my_opt$opt_d
    fix[it, 2] <- my_opt$sd_d
    fix[it, 3] <- my_opt$pred_y
    fix[it, 4] <- my_opt$pred_r
    fix[it, 5] <- my_opt$lambda_h
    
    ###lasso, with no penalty for dose
    
    cv_fit_y <- cv.glmnet(as.matrix(case1_X), case1$y, family = 'binomial', penalty.factor = c(0, rep(1, 2 * p)))
    y_coef <- as.vector(coef(cv_fit_y, s = 'lambda.min'))
    
    cv_fit_r <- cv.glmnet(as.matrix(case1_X), case1$r, family = 'binomial', penalty.factor = c(0, rep(1, 2 * p)))
    r_coef <- as.vector(coef(cv_fit_r, s = 'lambda.min'))
    
    lasso_s[it, ] <- c(y_coef, r_coef)
    my_opt <- opt.dose(0.2, r_coef, y_coef, N, case1_val)
    lasso[it, 1] <- my_opt$opt_d
    lasso[it, 2] <- my_opt$sd_d
    lasso[it, 3] <- my_opt$pred_y
    lasso[it, 4] <- my_opt$pred_r
    lasso[it, 5] <- my_opt$lambda_h
    
    
    ###classo,use 10-cv to choose best lambda with min deviance
    ##model y
    folds <- cut(seq(1, nrow(case1)), breaks = 10, labels = FALSE)
    dev <- matrix(NA, 10, 60)
    for(f in 1:10){
        testIndexes <- which(folds == f, arr.ind = TRUE)
        x <- scale(case1_X[-testIndexes, ]) 
        y <- case1$y[-testIndexes]
         
        mu_x <- apply(case1_X[-testIndexes, ], 2, mean)
        sd_x <- apply(case1_X[-testIndexes, ], 2, sd) 
         
        X <- cbind(1, x, -x[, -1]) #don't shrink dose
         
        for(j in 1:60){
            lambda <- cv_fit_y$lambda.min * exp((j - 30) * 0.1)
            #likelihood of logistic regression
            fr <- function(b){
                l <- rep(NA, n * 0.9)
                for(i in 1:(n * 0.9)){
                    l[i] <- y[i] * sum(X[i, ] * b) - log(1 + exp(sum(X[i, ] * b)))
                }
                lambda * sum(b[-c(1, 2)]) - mean(l)
            }
             
            grr <- function(b){
                d <- matrix(NA, n * 0.9, (4 * p + 2))
                for(i in 1:(n * 0.9)){
                    d[i, ]<- y[i] * X[i, ] - invlogit(sum(X[i, ] * b)) * X[i, ]
                }
                lambda * c(0, 0, rep(1, 4 * p)) - colMeans(d)
            }
             
            Amat <- rbind(cbind(0, 0, diag(1, (4 * p))), cbind(0, 1 / sd_x[1], matrix(0, n * 0.9, p), 
            	as.matrix(case1[-testIndexes, 4:(p + 3)]) / sd_x[(p + 2):(2 * p + 1)], matrix(0, n * 0.9, p), 
            	-as.matrix(case1[-testIndexes, 4:(p + 3)]) / sd_x[(p + 2):(2 * p + 1)]))
            b <- c(rep(0, n * 0.9 + 4 * p))
            meq <- 0
            coef <- spg(par = c(0.1, 0.1, rep(0.1, (2 * p)), rep(0.01, (2 * p))), fn = fr, gr = grr, 
            	control = list(trace = FALSE), project = 'projectLinear', projectArgs = list(A = Amat, b = b, meq = meq))$par
             
            # solution
            B_m <- cbind(diag(1, (2 * p + 1)), diag(-1, (2 * p + 1))[, -1]) #calculate beta
            y_coef <- c(coef[1] - sum(t(B_m %*% coef[-1]) / sd_x * mu_x), t(B_m %*% coef[-1]) / sd_x)
             
            #prediction of deviance
            dev[f, j] <- -2 * sum(case1$y[testIndexes] * (as.matrix(cbind(1, case1_X[testIndexes, 1:(2 * p + 1)])) %*% y_coef ) 
            	- log(1 + exp(as.matrix(cbind(1, case1_X[testIndexes, 1:(2 * p + 1)])) %*% y_coef)))
        }
    }
    
    lambda_min <- cv_fit_y$lambda.min * exp((which.min(colMeans(dev)) - 30) * 0.1)
    x <- scale(case1_X)
    y <- case1$y
    
    mu_x <- apply(case1_X, 2, mean)
    sd_x <- apply(case1_X, 2, sd) 
    
    X <- cbind(1, x, -x[, -1]) #don't shrink dose
    #likelihood of logistic regression
    fr <- function(b){
        l <- rep(NA, n)
        for(i in 1:n){
            l[i] <- y[i] * sum(X[i, ] * b) - log(1 + exp(sum(X[i, ] * b)))
        }
        lambda_min * sum(b[-c(1, 2)]) - mean(l)
    }
    
    grr <- function(b){
        d <- matrix(NA, n, (4 * p + 2))
        for(i in 1:n){
            d[i, ] <- y[i] * X[i, ] - invlogit(sum(X[i, ] * b)) * X[i, ]
        }
        lambda_min * c(0, 0, rep(1, 4 * p)) - colMeans(d)
    }
    
    Amat <- rbind(cbind(0, 0, diag(1, (4 * p))), cbind(0, 1 / sd_x[1], matrix(0, n, p), 
    	as.matrix(case1[, 4:(p + 3)]) / sd_x[(p + 2):(2 * p + 1)], matrix(0, n, p), 
    	-as.matrix(case1[, 4:(p + 3)]) / sd_x[(p + 2):(2 * p + 1)]))
    b <- c(rep(0, n + 4 * p))
    meq <- 0
    coef <- spg(par = c(0.1, 0.1, rep(0.1, (2 * p)), rep(0.01, (2 * p))), fn = fr, gr = grr, 
    	control = list(trace = FALSE), project = 'projectLinear', projectArgs = list(A = Amat, b = b, meq = meq))$par
    
    #solution
    B_m <- cbind(diag(1, (2 * p + 1)), diag(-1, (2 * p + 1))[, -1]) #calculate beta
    y_coef <- c(coef[1] - sum(t(B_m %*% coef[-1]) / sd_x * mu_x), t(B_m %*% coef[-1]) / sd_x)
    
    classo_s[it, 1:(2 * p + 3)] <- c(which.min(colMeans(dev)), y_coef)
    
    ##model r
    folds <- cut(seq(1, nrow(case1)), breaks = 10, labels = FALSE)
    dev <- matrix(NA, 10, 60)
    for(f in 1:10){
        testIndexes <- which(folds == f, arr.ind = TRUE)
        x <- scale(case1_X[-testIndexes, ]) 
        y <- case1$r[-testIndexes]
        
        mu_x <- apply(case1_X[-testIndexes, ], 2, mean)
        sd_x <- apply(case1_X[-testIndexes, ], 2, sd) 
        
        X <- cbind(1, x, -x[, -1]) #don't shrink dose
        
        for(j in 1:60){
            lambda <- cv_fit_r$lambda.min * exp((j - 30) * 0.1)
            #likelihood of logistic regression
            fr <- function(b){
                l <- rep(NA, n * 0.9)
                for(i in 1:(n * 0.9)){
                    l[i] <- y[i] * sum(X[i, ] * b) - log(1 + exp(sum(X[i, ] * b)))
                }
                lambda * sum(b[-c(1, 2)]) - mean(l)
            }
            
            grr <- function(b){
                d <- matrix(NA, n * 0.9, (4 * p + 2))
                for(i in 1:(n * 0.9)){
                    d[i, ]<- y[i] * X[i, ] - invlogit(sum(X[i, ] * b)) * X[i, ]
                }
                lambda * c(0, 0, rep(1, 4 * p)) - colMeans(d)
            }
            
            Amat <- rbind(cbind(0, 0, diag(1, (4 * p))), cbind(0, 1 / sd_x[1], matrix(0, n * 0.9, p), 
            	as.matrix(case1[-testIndexes, 4:(p + 3)]) / sd_x[(p + 2):(2 * p + 1)], matrix(0, n * 0.9, p), 
            	-as.matrix(case1[-testIndexes, 4:(p + 3)]) / sd_x[(p + 2):(2 * p + 1)]))
            b <- c(rep(0, n * 0.9 + 4 * p))
            meq <- 0
            coef <- spg(par = c(0.1, 0.1, rep(0.1, (2 * p)), rep(0.01, (2 * p))), fn = fr, gr = grr, 
            	control = list(trace = FALSE), project = 'projectLinear', projectArgs = list(A = Amat, b = b, meq = meq))$par
            
            # solution
            B_m <- cbind(diag(1, (2 * p + 1)), diag(-1, (2 * p + 1))[, -1]) #calculate beta
            r_coef <- c(coef[1] - sum(t(B_m %*% coef[-1]) / sd_x * mu_x), t(B_m %*% coef[-1]) / sd_x)
            
            #prediction of deviance
            dev[f, j] <- -2 * sum(case1$r[testIndexes] * (as.matrix(cbind(1, case1_X[testIndexes, 1:(2 * p + 1)])) %*% r_coef ) 
            	- log(1 + exp(as.matrix(cbind(1, case1_X[testIndexes, 1:(2 * p + 1)])) %*% r_coef)))
        }
    }
    
    lambda_min <- cv_fit_r$lambda.min * exp((which.min(colMeans(dev)) - 30) * 0.1)
    x <- scale(case1_X)
    y <- case1$r
    
    mu_x <- apply(case1_X, 2, mean)
    sd_x <- apply(case1_X, 2, sd) 
    
    X <- cbind(1, x, -x[, -1]) #don't shrink dose
    #likelihood of logistic regression
    fr <- function(b){
        l <- rep(NA, n)
        for(i in 1:n){
            l[i] <- y[i] * sum(X[i, ] * b) - log(1 + exp(sum(X[i, ] * b)))
        }
        lambda_min * sum(b[-c(1, 2)]) - mean(l)
    }
    
    grr <- function(b){
        d <- matrix(NA, n, (4 * p + 2))
        for(i in 1:n){
            d[i, ] <- y[i] * X[i, ] - invlogit(sum(X[i, ] * b)) * X[i, ]
        }
        lambda_min * c(0, 0, rep(1, 4 * p)) - colMeans(d)
    }
    
    Amat <- rbind(cbind(0, 0, diag(1, (4 * p))), cbind(0, 1 / sd_x[1], matrix(0, n, p), 
    	as.matrix(case1[, 4:(p + 3)]) / sd_x[(p + 2):(2 * p + 1)], matrix(0, n, p), 
    	-as.matrix(case1[, 4:(p + 3)]) / sd_x[(p + 2):(2 * p + 1)]))
    b <- c(rep(0, n + 4 * p))
    meq <- 0
    coef <- spg(par = c(0.1, 0.1, rep(0.1, (2 * p)), rep(0.01, (2 * p))), fn = fr, gr = grr, 
    	control = list(trace = FALSE), project = 'projectLinear', projectArgs = list(A = Amat, b = b, meq = meq))$par
    
    #solution
    B_m <- cbind(diag(1, (2 * p + 1)), diag(-1, (2 * p + 1))[, -1]) #calculate beta
    r_coef <- c(coef[1] - sum(t(B_m %*% coef[-1]) / sd_x * mu_x), t(B_m %*% coef[-1]) / sd_x)
    
    classo_s[it, (2 * p + 4):(4 * p + 6)] <- c(which.min(colMeans(dev)), r_coef)
    
    my_opt <- opt.dose(0.2, r_coef, y_coef, N, case1_val)
    classo[it, 1] <- my_opt$opt_d
    classo[it, 2] <- my_opt$sd_d
    classo[it, 3] <- my_opt$pred_y
    classo[it, 4] <- my_opt$pred_r
    classo[it, 5] <- my_opt$lambda_h

}

simulation <- data.frame(cbind(theory, ignore, min, mid, fs, lasso, classo, fix))
selection <- data.frame(cbind(fs_s, lasso_s, classo_s))

###efficacy
boxplot(cbind(theory[,3],fs[,3],lasso[,3],classo[,3],fix[,3]),names=c("Theory","FS","LASSO","cLASSO","FD"), ylab="Mean E")
points(colMeans(cbind(theory[,3],fs[,3],lasso[,3],classo[,3],fix[,3]),na.rm=T),pch=18)
round(colMeans(cbind(theory[,3],fs[,3],lasso[,3],classo[,3],fix[,3]),na.rm=T),3)

###percentage of improvement
round(colMeans((cbind(fs[,3],lasso[,3],classo[,3])-fix[,3])/(theory[,3]-fix[,3]),na.rm=T),3)

###how often CLASSO beats LASSO
summary(classo[,3]>lasso[,3])

###toxicity
round(colMeans(cbind(theory[,4],fs[,4],lasso[,4],classo[,4],fix[,4]),na.rm=T),3)

