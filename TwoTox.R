###################################################################################
#                                                                                 #
#   Filename    :       TwoTox.R                          						  #						                                                                                     
#   Project     :       BiomJ article "A Utility Approach to Individualized       #
#                       Optimal Dose Selection Using Biomarkers "   				  #
#   Authors     :       Pin Li                                                    #    
#   Date        :       8.13.2019                                                 #
#   Required R packages :  tmvtnorm, stats, glmnet, BB, numDeriv                  #
###################################################################################

rm(list=ls())

# # to run parallel simulations on the cluster, grab the array id value from the environment variable passed from sbatch
# slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# 
# # coerce the value to an integer
# taskid <- as.numeric(slurm_arrayid)
# set.seed(taskid)

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
	x <- rmvnorm(n * 3, mean = rep(0, p))
    x <- x[abs(0.4 * rowSums(x[, 1:3]) - 0.8 * x[, 4]) < 1 & (0.5 * x[,1]-0.5 * x[,5])<1, ]
    x <- x[1:n, ]
    
    a <- runif(n, -1, 1)
    
    h_y <- (1 + 0.4 * rowSums(x[, 1:3]) - 0.8 * x[, 4]) * a
    h_r1 <- (1 - 0.4 * rowSums(x[, 1:3]) + 0.8 * x[, 4]) * a
    h_r2 <- (1 - 0.5 * x[,1] + 0.5 * x[, 5]) * a
    
    y <- rbinom(n, 1, invlogit(x[, 1] + h_y))
    r1 <- rbinom(n, 1, invlogit(-1.386 - x[, 1] + h_r1))
    r2 <- rbinom(n, 1, invlogit(-1.2 - x[, 1] + h_r2))
    
    data.frame(y, r1, r2, a, x)
}



###### get the matrix with interactions
design <- function(mat){
    int_mat <- model.matrix( ~ . ^ 2, data = data.frame(mat))
    return(int_mat[, 2:(2 * p + 2)])
}

###### random walk to find best lambda
opt.dose<-function(tau1,tau2, r1_coef,r2_coef, y_coef, N, case1_val){
    
    
    y_coef_t <- c(0, 1, 1, rep(0, p - 1), rep(0.4, 3), -0.8, rep(0, p - 4))
	r1_coef_t <- c(-1.386, 1, -1, rep(0, p - 1), rep(-0.4, 3), 0.8, rep(0, p - 4))
    r2_coef_t <- c(-1.2, 1, -1, rep(0, p - 1), -0.5, rep(0, 3), 0.5, rep(0, p - 5))
    
    p_y_0 <- invlogit(cbind(1, design(cbind(a = -1, case1_val[, -c(1:4)]))) %*% y_coef)
    p_r1_0 <- invlogit(cbind(1, design(cbind(a = -1, case1_val[, -c(1:4)]))) %*% r1_coef)
    p_r2_0 <- invlogit(cbind(1, design(cbind(a = -1, case1_val[, -c(1:4)]))) %*% r2_coef)
    
    
    if(mean(p_r1_0) > tau1 | mean(p_r2_0) > tau2){
        return(list(lambda_h = c(NA,NA), opt_d = NA, pred_y = NA, pred_r1 = NA, pred_r2 = NA))
    }else if (mean(invlogit(cbind(1, design(cbind(a = 1, case1_val[, -c(1:4)]))) %*% r1_coef))<tau1
              & mean(invlogit(cbind(1, design(cbind(a = 1, case1_val[, -c(1:4)]))) %*% r2_coef))<tau2){
        return(list(lambda_h =c(0, 0), opt_d = 1, 
                 pred_y = mean(invlogit(cbind(1, design(cbind(a = 1, case1_val[, -c(1:4)]))) %*% y_coef)),
                 pred_r1 = mean(invlogit(cbind(1, design(cbind(a = 1, case1_val[, -c(1:4)]))) %*% r1_coef)),
                 pred_r2 = mean(invlogit(cbind(1, design(cbind(a = 1, case1_val[, -c(1:4)]))) %*% r2_coef))))
    }else{
        my.opt<-function(lambda1,lambda2){
            D<-rep(NA,N)
            for (j in 1:N){
                f<-function(a) invlogit(c(1,design(cbind(a,case1_val[,-c(1:4)]))[j,]) %*% y_coef)-
                    lambda1*invlogit(c(1,design(cbind(a,case1_val[,-c(1:4)]))[j,]) %*% r1_coef)-           
                    lambda2*invlogit(c(1,design(cbind(a,case1_val[,-c(1:4)]))[j,]) %*% r2_coef)-
                    p_y_0[j]+lambda1*p_r1_0[j]+lambda2*p_r2_0[j]
                D[j]<-optimize(f, c(-1, 1), tol = 0.001,maximum = TRUE)$maximum
            }
            return(D)
        }
        
        
        #random walk
        save<-matrix(NA,1001,6)
        lambda<-c(1,1) #initial value
        D_d<-my.opt(lambda[1],lambda[2])
        r1_d<-mean(invlogit(cbind(1,design(cbind(a=D_d,case1_val[,-c(1:4)]))) %*% r1_coef_t))
        r2_d<-mean(invlogit(cbind(1,design(cbind(a=D_d,case1_val[,-c(1:4)]))) %*% r2_coef_t))
        y_d<-mean(invlogit(cbind(1,design(cbind(a=D_d,case1_val[,-c(1:4)]))) %*% y_coef_t))
        save[1,]<-c(lambda,y_d,r1_d,r2_d,mean(D_d))
        
        for(i in 1:1000){
            lambda_s<-rtmvnorm(n=1, as.vector(lambda), diag(c(0.04,0.04)), upper=c(2,2),lower=c(0,0))
            D_s<-my.opt(lambda_s[1],lambda_s[2])
            r1_d_s<-mean(invlogit(cbind(1,design(cbind(a=D_s,case1_val[,-c(1:4)]))) %*% r1_coef_t))
            r2_d_s<-mean(invlogit(cbind(1,design(cbind(a=D_s,case1_val[,-c(1:4)]))) %*% r2_coef_t))
            y_d_s<-mean(invlogit(cbind(1,design(cbind(a=D_s,case1_val[,-c(1:4)]))) %*% y_coef_t))
            u<-runif(1)
            
            if( (r1_d_s-tau1<0.005) & (r2_d_s-tau2<0.005) & (u<= min(y_d_s/y_d,1))){
                lambda<-lambda_s
                y_d<-y_d_s}
            save[(i+1),]<-c(lambda,y_d,r1_d_s,r2_d_s,mean(D_s))
            i<-i+1
        }
        index<-which.max(save[,3])
        return(list(lambda_h=save[index,1:2],
                    opt_d=save[index,6],
                    pred_y=save[index,3],
                    pred_r1=save[index,4],
                    pred_r2=save[index,5]))
    }
}

        
###### for illustration, run 10 replications
it_n <- 1

###### scenario 0  as example
p <- 5

theory <- matrix(NA, it_n, 6)

ignore <- matrix(NA, it_n, 3)

min <- matrix(NA, it_n, 3)

mid <- matrix(NA, it_n, 3)

fix <- matrix(NA, it_n, 6)

fs <- matrix(NA, it_n, 6)

lasso <- matrix(NA, it_n, 6)

classo <- matrix(NA, it_n, 6)

for (it in 1:it_n){

    n <- 200
    N <- n
    
    
    case1 <- my.sample(n, p)
    case1_val <- case1
    
    case1_X <- data.frame(design(case1[, -c(1:3)]))
    
    ###theoretic
    y_coef <- c(0, 1, 1, rep(0, p - 1), rep(0.4, 3), -0.8, rep(0, p - 4))
	r1_coef <- c(-1.386, 1, -1, rep(0, p - 1), rep(-0.4, 3), 0.8, rep(0, p - 4))
    r2_coef <- c(-1.2, 1, -1, rep(0, p - 1), -0.5, rep(0, 3), 0.5, rep(0, p - 5))
    
    ###if ignore risk
    ignore[it, 1] <- mean(invlogit(cbind(1, design(cbind(a = 1, case1_val[, -c(1:4)]))) %*% y_coef))
    ignore[it, 2] <- mean(invlogit(cbind(1, design(cbind(a = 1, case1_val[, -c(1:4)]))) %*% r1_coef))
    ignore[it, 3] <- mean(invlogit(cbind(1, design(cbind(a = 1, case1_val[, -c(1:4)]))) %*% r2_coef))
    
    ###if choose min dose
    min[it, 1] <- mean(invlogit(cbind(1, design(cbind(a = -1, case1_val[, -c(1:4)]))) %*% y_coef))
    min[it, 2] <- mean(invlogit(cbind(1, design(cbind(a = -1, case1_val[, -c(1:4)]))) %*% r1_coef))
    min[it, 3] <- mean(invlogit(cbind(1, design(cbind(a = -1, case1_val[, -c(1:4)]))) %*% r2_coef))
    
    ###if choose mid dose
    mid[it, 1] <- mean(invlogit(cbind(1, design(cbind(a = 0, case1_val[, -c(1:4)]))) %*% y_coef))
    mid[it, 2] <- mean(invlogit(cbind(1, design(cbind(a = 0, case1_val[, -c(1:4)]))) %*% r1_coef))
    mid[it, 3] <- mean(invlogit(cbind(1, design(cbind(a = 0, case1_val[, -c(1:4)]))) %*% r2_coef))
    
    ###theoretic optimal dose
    my_opt<-opt.dose(0.2, invlogit(-1.2), r1_coef,r2_coef,y_coef, N, case1_val)
    theory[it,1]<-my_opt$opt_d
    theory[it,2]<-my_opt$pred_y
    theory[it,3]<-my_opt$pred_r1
    theory[it,4]<-my_opt$pred_r2
    theory[it,5:6]<-my_opt$lambda_h
    
    ###Fit regression model for Y given (A, X, AX), R given (A, X, AX)
    md_y <- glm(case1$y ~ . , data = case1_X, family = 'binomial')
    
    md_r1 <- glm(case1$r1 ~ . , data = case1_X, family = 'binomial')
    
    md_r2 <- glm(case1$r2 ~ . , data = case1_X, family = 'binomial')
    
    ###forward selection, force DOSE in
    null_y <- glm(case1$y ~ a, data = case1_X, family = 'binomial')
    y_forwards <- step(null_y, scope = list(lower = null_y, upper = md_y), direction = 'forward', trace = 0) 
    y_coef <- rep(0, 2 * p + 2)
    y_coef[match(names(y_forwards$coefficients), names(md_y$coeff))] <- y_forwards$coefficients
    
    null_r1 <- glm(case1$r1 ~ a, data = case1_X, family = 'binomial')
    r1_forwards <- step(null_r1, scope = list(lower = null_r1, upper = md_r1), direction = 'forward', trace = 0) 
    r1_coef <- rep(0, 2 * p + 2)
    r1_coef[match(names(r1_forwards$coefficients), names(md_r1$coeff))] <- r1_forwards$coefficients
    
    
    null_r2 <- glm(case1$r2 ~ a, data = case1_X, family = 'binomial')
    r2_forwards <- step(null_r2, scope = list(lower = null_r2, upper = md_r2), direction = 'forward', trace = 0) 
    r2_coef <- rep(0, 2 * p + 2)
    r2_coef[match(names(r2_forwards$coefficients), names(md_r2$coeff))] <- r2_forwards$coefficients
    
    
    my_opt<-opt.dose(0.2, invlogit(-1.2), r1_coef,r2_coef,y_coef, N, case1_val)
    fs[it,1]<-my_opt$opt_d
    fs[it,2]<-my_opt$pred_y
    fs[it,3]<-my_opt$pred_r1
    fs[it,4]<-my_opt$pred_r2
    fs[it,5:6]<-my_opt$lambda_h
    
    
    
    ###fixed model with only dose
    y_coef <- c(null_y$coefficients, rep(0, 2 * p))
    
    r1_coef <- c(null_r1$coefficients, rep(0, 2 * p))
    
    r2_coef <- c(null_r2$coefficients, rep(0, 2 * p))
    
    my_opt<-opt.dose(0.2, invlogit(-1.2), r1_coef,r2_coef,y_coef, N, case1_val)
    fix[it,1]<-my_opt$opt_d
    fix[it,2]<-my_opt$pred_y
    fix[it,3]<-my_opt$pred_r1
    fix[it,4]<-my_opt$pred_r2
    fix[it,5:6]<-my_opt$lambda_h
    
    
    ###lasso, with no penalty for dose
    
    cv_fit_y <- cv.glmnet(as.matrix(case1_X), case1$y, family = 'binomial', penalty.factor = c(0, rep(1, 2 * p)))
    y_coef <- as.vector(coef(cv_fit_y, s = 'lambda.min'))
    
    cv_fit_r1 <- cv.glmnet(as.matrix(case1_X), case1$r1, family = 'binomial', penalty.factor = c(0, rep(1, 2 * p)))
    r1_coef <- as.vector(coef(cv_fit_r1, s = 'lambda.min'))
    
    cv_fit_r2 <- cv.glmnet(as.matrix(case1_X), case1$r2, family = 'binomial', penalty.factor = c(0, rep(1, 2 * p)))
    r2_coef <- as.vector(coef(cv_fit_r2, s = 'lambda.min'))
    

    my_opt<-opt.dose(0.2, invlogit(-1.2), r1_coef,r2_coef,y_coef, N, case1_val)
    lasso[it,1]<-my_opt$opt_d
    lasso[it,2]<-my_opt$pred_y
    lasso[it,3]<-my_opt$pred_r1
    lasso[it,4]<-my_opt$pred_r2
    lasso[it,5:6]<-my_opt$lambda_h
    
    
    
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
                                                               as.matrix(case1_X[-testIndexes, 2:(p + 1)]) / sd_x[(p + 2):(2 * p + 1)], matrix(0, n * 0.9, p), 
                                                               -as.matrix(case1_X[-testIndexes, 2:(p + 1)]) / sd_x[(p + 2):(2 * p + 1)]))
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
                                                       as.matrix(case1_X[, 2:(p + 1)]) / sd_x[(p + 2):(2 * p + 1)], matrix(0, n, p), 
                                                       -as.matrix(case1_X[, 2:(p + 1)]) / sd_x[(p + 2):(2 * p + 1)]))
    b <- c(rep(0, n + 4 * p))
    meq <- 0
    coef <- spg(par = c(0.1, 0.1, rep(0.1, (2 * p)), rep(0.01, (2 * p))), fn = fr, gr = grr, 
                control = list(trace = FALSE), project = 'projectLinear', projectArgs = list(A = Amat, b = b, meq = meq))$par
    
    #solution
    B_m <- cbind(diag(1, (2 * p + 1)), diag(-1, (2 * p + 1))[, -1]) #calculate beta
    y_coef <- c(coef[1] - sum(t(B_m %*% coef[-1]) / sd_x * mu_x), t(B_m %*% coef[-1]) / sd_x)

    ##model r1
    folds <- cut(seq(1, nrow(case1)), breaks = 10, labels = FALSE)
    dev <- matrix(NA, 10, 60)
    for(f in 1:10){
        testIndexes <- which(folds == f, arr.ind = TRUE)
        x <- scale(case1_X[-testIndexes, ]) 
        y <- case1$r1[-testIndexes]
        
        mu_x <- apply(case1_X[-testIndexes, ], 2, mean)
        sd_x <- apply(case1_X[-testIndexes, ], 2, sd) 
        
        X <- cbind(1, x, -x[, -1]) #don't shrink dose
        
        for(j in 1:60){
            lambda <- cv_fit_r1$lambda.min * exp((j - 30) * 0.1)
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
                                                               as.matrix(case1_X[-testIndexes, 2:(p + 1)]) / sd_x[(p + 2):(2 * p + 1)], matrix(0, n * 0.9, p), 
                                                               -as.matrix(case1_X[-testIndexes, 2:(p + 1)]) / sd_x[(p + 2):(2 * p + 1)]))
            b <- c(rep(0, n * 0.9 + 4 * p))
            meq <- 0
            coef <- spg(par = c(0.1, 0.1, rep(0.1, (2 * p)), rep(0.01, (2 * p))), fn = fr, gr = grr, 
                        control = list(trace = FALSE), project = 'projectLinear', projectArgs = list(A = Amat, b = b, meq = meq))$par
            
            # solution
            B_m <- cbind(diag(1, (2 * p + 1)), diag(-1, (2 * p + 1))[, -1]) #calculate beta
            r_coef <- c(coef[1] - sum(t(B_m %*% coef[-1]) / sd_x * mu_x), t(B_m %*% coef[-1]) / sd_x)
            
            #prediction of deviance
            dev[f, j] <- -2 * sum(case1$r1[testIndexes] * (as.matrix(cbind(1, case1_X[testIndexes, 1:(2 * p + 1)])) %*% r_coef ) 
                                  - log(1 + exp(as.matrix(cbind(1, case1_X[testIndexes, 1:(2 * p + 1)])) %*% r_coef)))
        }
    }
    
    lambda_min <- cv_fit_r1$lambda.min * exp((which.min(colMeans(dev)) - 30) * 0.1)
    x <- scale(case1_X)
    y <- case1$r1
    
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
                                                       as.matrix(case1_X[, 2:(p + 1)]) / sd_x[(p + 2):(2 * p + 1)], matrix(0, n, p), 
                                                       -as.matrix(case1_X[, 2:(p + 1)]) / sd_x[(p + 2):(2 * p + 1)]))
    b <- c(rep(0, n + 4 * p))
    meq <- 0
    coef <- spg(par = c(0.1, 0.1, rep(0.1, (2 * p)), rep(0.01, (2 * p))), fn = fr, gr = grr, 
                control = list(trace = FALSE), project = 'projectLinear', projectArgs = list(A = Amat, b = b, meq = meq))$par
    
    #solution
    B_m <- cbind(diag(1, (2 * p + 1)), diag(-1, (2 * p + 1))[, -1]) #calculate beta
    r1_coef <- c(coef[1] - sum(t(B_m %*% coef[-1]) / sd_x * mu_x), t(B_m %*% coef[-1]) / sd_x)
    
    ##model r2
    folds <- cut(seq(1, nrow(case1)), breaks = 10, labels = FALSE)
    dev <- matrix(NA, 10, 60)
    for(f in 1:10){
        testIndexes <- which(folds == f, arr.ind = TRUE)
        x <- scale(case1_X[-testIndexes, ]) 
        y <- case1$r2[-testIndexes]
        
        mu_x <- apply(case1_X[-testIndexes, ], 2, mean)
        sd_x <- apply(case1_X[-testIndexes, ], 2, sd) 
        
        X <- cbind(1, x, -x[, -1]) #don't shrink dose
        
        for(j in 1:60){
            lambda <- cv_fit_r2$lambda.min * exp((j - 30) * 0.1)
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
                                                               as.matrix(case1_X[-testIndexes, 2:(p + 1)]) / sd_x[(p + 2):(2 * p + 1)], matrix(0, n * 0.9, p), 
                                                               -as.matrix(case1_X[-testIndexes, 2:(p + 1)]) / sd_x[(p + 2):(2 * p + 1)]))
            b <- c(rep(0, n * 0.9 + 4 * p))
            meq <- 0
            coef <- spg(par = c(0.1, 0.1, rep(0.1, (2 * p)), rep(0.01, (2 * p))), fn = fr, gr = grr, 
                        control = list(trace = FALSE), project = 'projectLinear', projectArgs = list(A = Amat, b = b, meq = meq))$par
            
            # solution
            B_m <- cbind(diag(1, (2 * p + 1)), diag(-1, (2 * p + 1))[, -1]) #calculate beta
            r_coef <- c(coef[1] - sum(t(B_m %*% coef[-1]) / sd_x * mu_x), t(B_m %*% coef[-1]) / sd_x)
            
            #prediction of deviance
            dev[f, j] <- -2 * sum(case1$r2[testIndexes] * (as.matrix(cbind(1, case1_X[testIndexes, 1:(2 * p + 1)])) %*% r_coef ) 
                                  - log(1 + exp(as.matrix(cbind(1, case1_X[testIndexes, 1:(2 * p + 1)])) %*% r_coef)))
        }
    }
    
    lambda_min <- cv_fit_r2$lambda.min * exp((which.min(colMeans(dev)) - 30) * 0.1)
    x <- scale(case1_X)
    y <- case1$r2
    
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
                                                       as.matrix(case1_X[, 2:(p + 1)]) / sd_x[(p + 2):(2 * p + 1)], matrix(0, n, p), 
                                                       -as.matrix(case1_X[, 2:(p + 1)]) / sd_x[(p + 2):(2 * p + 1)]))
    b <- c(rep(0, n + 4 * p))
    meq <- 0
    coef <- spg(par = c(0.1, 0.1, rep(0.1, (2 * p)), rep(0.01, (2 * p))), fn = fr, gr = grr, 
                control = list(trace = FALSE), project = 'projectLinear', projectArgs = list(A = Amat, b = b, meq = meq))$par
    
    #solution
    B_m <- cbind(diag(1, (2 * p + 1)), diag(-1, (2 * p + 1))[, -1]) #calculate beta
    r2_coef <- c(coef[1] - sum(t(B_m %*% coef[-1]) / sd_x * mu_x), t(B_m %*% coef[-1]) / sd_x)
    
    
    
    my_opt<-opt.dose(0.2, invlogit(-1.2), r1_coef,r2_coef,y_coef, N, case1_val)
    classo[it,1]<-my_opt$opt_d
    classo[it,2]<-my_opt$pred_y
    classo[it,3]<-my_opt$pred_r1
    classo[it,4]<-my_opt$pred_r2
    classo[it,5:6]<-my_opt$lambda_h
    
}

simulation <- data.frame(cbind(theory, ignore, min, mid, fs, lasso, classo, fix))

###efficacy
round(colMeans(cbind(theory[,2],fs[,2],lasso[,2],classo[,2],fix[,2]),na.rm=T),3)

###percentage of improvement
round(colMeans((cbind(fs[,2],lasso[,2],classo[,2])-fix[,2])/(theory[,2]-fix[,2]),na.rm=T),3)

###how often CLASSO beats LASSO
summary(classo[,2]>lasso[,2])

###toxicity
round(colMeans(cbind(theory[,3],fs[,3],lasso[,3],classo[,3],fix[,3]),na.rm=T),3)
round(colMeans(cbind(theory[,4],fs[,4],lasso[,4],classo[,4],fix[,4]),na.rm=T),3)
