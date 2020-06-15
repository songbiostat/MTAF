#  ----------- PART I: Simulations ---------------
library(aSPU)
library(invgamma)
library(mvtnorm)
library(statmod)
library(MSKAT)
library(CompQuadForm)
library(MultiPhen)
source("tates.R")

# ------------------ PART 1: Simulation for type I & II error -------------

##### Test power with correlated Z & CS

SIM_covariate_power_Z <- function(betaz=1,constr="ar",cor,beta,size=1000,traits=10,nonzero_ratio=0.2,cat_ratio=0,MAF=0.3,n.cov=1){
  
  n = size # sample size
  k = traits # traits
  cat <- round(k*cat_ratio)
  
  signif = round(k*nonzero_ratio)  # # of significant beta
  
  #k.con = 1 # continuous traits
  k.cat = round(signif*cat_ratio) # categorical traits
  
  
  # generate genotye score
  
  X = rbinom(n,2,MAF)
  
  # generate covariates
  

  Z <-  matrix(rnorm(n*n.cov),ncol = n.cov)
 # Z <- cbind(rbinom(n,1,0.6),rbinom(n,1,0.5))
  beta_Z = matrix(runif(n.cov*k,0,betaz),ncol = n.cov)
  
  
  # generate Y's 
  
  if(constr == "ar"){
    pho = cor ^ as.matrix(dist(1:k))   # correlation matrix AR(1) with 0.9
  }else if(constr == "cs"){
    pho = rep(1,k) %*% t(rep(1,k)) * cor # Compound Symmetry
    diag(pho) <- 1
  }
  
  var = sqrt(rinvgamma(k,4,4))    # arbitrary variance for 10 traits
  Sigma = matrix((var %*% t(var)) * pho,ncol = k)   # recover covariance matrix
  
  epsilon = rmvnorm(n,rep(0,k),Sigma)   # generate correlated random error 
  
  Y <- X %o% beta + epsilon + Z %*% t(beta_Z)  # generate Y's 
  
  if(cat_ratio > 0 ){
    col.cat <- sample(1:signif,k.cat)
    col.cat.notsig <- sample((signif+1):traits,cat-k.cat)
    Y_cat <- as.numeric(Y[,col.cat])
    Y_bin <- matrix(rbinom(length(Y_cat),1,prob = 1/(1+exp(-Y_cat))),ncol = k.cat)
    Y[,col.cat] <- Y_bin
    
    Y_cat <- as.numeric(Y[,col.cat.notsig])
    Y_bin <- matrix(rbinom(length(Y_cat),1,prob = 1/(1+exp(-Y_cat))),ncol = cat-k.cat)
    Y[,col.cat.notsig] <- Y_bin
    
    return(list(Y=Y,X=X,Z=Z,beta=beta,epsilon=epsilon,col_binary=sort(c(col.cat,col.cat.notsig))))
  }else{
    return(list(Y=Y,X=X,Z=Z,beta=beta,epsilon=epsilon))
  }
  
  
}


SIM_binary_power_noZ <- function(constr="cs",cor,beta,size=1000,traits=10,MAF=0.3,n.cov=2){
  
  n = size # sample size
  k = traits # traits
  
  # generate genotye score
  
  X = rbinom(n,2,MAF)
  
  
  # generate Y's 
  
  if(constr == "ar"){
    pho = cor ^ as.matrix(dist(1:k))   # correlation matrix AR(1) with rho
  }else if(constr == "cs"){
    pho = rep(1,k) %*% t(rep(1,k)) * cor # Compound Symmetry
    diag(pho) <- 1
  }
  
  var = sqrt(rinvgamma(k,4,4))    # arbitrary variance for 10 traits
  Sigma = matrix((var %*% t(var)) * pho,ncol = k)   # recover covariance matrix
  
  epsilon = rmvnorm(n,rep(0,k),Sigma)   # generate correlated random error 
  
  
  Y <-  X %o% beta + epsilon # generate Y's 
  
  Y.bin <- matrix(rbinom(n*k,1,prob = 1/(1+exp(-Y))),ncol = k)
  
  return(list(Y=Y.bin,Y.raw=Y,X=X,beta=beta,epsilon=epsilon,pho=pho))
  
}




# -------------------- PART 2: TEST & METHOD -------------------------



########################## SCORE TEST ############################

score_glm <- function(Y,X,cov, col_binary=NULL,nperm=1e3,one_side=TRUE){
  n <- nrow(X)
  k <- ncol(Y)
  
  #glm_XonZ <- glm(X ~ cov, family = "gaussian")
  if(is.null(cov)){
    X <- X
  }else{
    X <- resid(lm(X ~ cov))
  }

  
  
  # Create Matrix to save output
  X_perm <- cbind(X, replicate(nperm, sample(X)))
  
  if(one_side){
    p_l <- matrix(NA,ncol = k,nrow = nperm+1)
    p_u <- matrix(NA,ncol = k,nrow = nperm+1)
  }else{
    p_l <- matrix(NA,ncol = k,nrow = nperm+1)
  }
  
  # Fit null model before permutation
  if(length(col_binary)>0){
    for (j in col_binary) {
      assign(paste0("fit_cat_",j),glm(Y[,j]~1, family = "binomial"))
    }
    k_cont <- (1:k)[-col_binary]
  }else{
    k_cont <- 1:k
  }
  
  fit.cont <- list()
  for (i in k_cont) {
    fit.cont[[i]] <- glm(Y[,i] ~ 1, family = "gaussian")
  }
  
  ## Permutation start
  score <- matrix(NA, ncol = k, nrow = nperm + 1)
  for (h in 1:k) {
    score[, h] <- glm.scoretest(fit.cont[[h]], X_perm)
  }
  
  p_l <- pnorm(score)
  p_u <- pnorm(score, lower.tail = FALSE)
  
  return(list("p_l"=p_l,"p_u"=p_u))
}

# glm score test for binary traits
score_glm_bin <- function(Y,X,nperm=1e3,one_side=TRUE,overall=FALSE,X_permuted=NULL){
  n <- nrow(X)
  k <- ncol(Y)
  
  #glm_XonZ <- glm(X ~ cov, family = "gaussian")
  
  # Create Matrix to save output
  if(!overall){
    X_perm <- cbind(X, replicate(nperm, sample(X)))
  }else{
    X_perm <- X_permuted
  }
  
  if(one_side){
    p_l <- matrix(NA,ncol = k,nrow = nperm+1)
    p_u <- matrix(NA,ncol = k,nrow = nperm+1)
  }else{
    p_l <- matrix(NA,ncol = k,nrow = nperm+1)
  }
  
  # Fit null model before permutation
  
  k_cont <- 1:k
  
  fit.cont <- list()
  for (i in k_cont) {
    fit.cont[[i]] <- glm(Y[,i] ~ 1, family = "binomial")
  }
  
  ## Permutation start
  score <- matrix(NA, ncol = k, nrow = nperm + 1)
  for (h in 1:k) {
    score[, h] <- glm.scoretest(fit.cont[[h]], X_perm)
  }
  
  p_l <- pnorm(score)
  p_u <- pnorm(score, lower.tail = FALSE)
  
  return(list("p_l"=p_l,"p_u"=p_u))
}


## SCORE TEST: NO COVARIATES

score_test <- function(Y,X,sample_size=1000, nperm=1000,one_side=TRUE,binary=FALSE){
  n <- sample_size
  y <- Y
  x <- X
  
  if(nperm>0){
    X_perm <- cbind(X, replicate(nperm, sample(X,n,replace = FALSE)))
    x <- X_perm
  }
  
  
  xres <- scale(x,scale = FALSE)
  yres <- scale(y,scale = FALSE)
  
  U <- t(x) %*% scale(y,scale = FALSE)
 
  if(binary){
    v <-  colMeans(y) %*% t(colMeans(y)) *  as.numeric((t(xres[,1]) %*% xres[,1]))  # variance is same across permutations
    v <- diag(v)
  }else{
    v <-  cov(y) *  as.numeric((t(xres[,1]) %*% xres[,1]))  # variance is same across permutations
    v <- diag(v)
  }

 
  s <- t(t(U)/sqrt(v))
  
  if(one_side == FALSE){
    pval <- 2*(1-pnorm(abs(s)))
    return(list(U=U,U_cov=v,p_value=pval))
  }else{
    p_lower <- pnorm(s, lower.tail = TRUE)
    p_upper <- pnorm(s, lower.tail = FALSE)
    
    return(list(U=U,U_cov=v,U.standard=s,p_u=p_upper,p_l=p_lower))
  }
  
}

######################### AF METHOD ########################


##### COMPUTE AF P-VALUES
AdaptiveFisher_compute <- function(p,log=FALSE){
  n <- nrow(p)
  k <- ncol(p)
  if(log){
    p <- -p
  }else{
    p <- -log(p)
  }
  
  p_sort <- apply(p, 1, sort,decreasing=TRUE)
  p_cumsum <- apply(p_sort, 2, cumsum)
  
  # permutation p_value
  pval_rank <- apply(-p_cumsum, 1, rank,ties.method="max")
  pval_perm <- pval_rank/n
  AF_stat <- apply(pval_perm, 1, min)    # tuning parameter?
  AF_stat_p <- rank(AF_stat, ties.method = "max")/n
  AF_perm_pval <- AF_stat_p[-1]
  AF_pval <- AF_stat_p[1]
  return(list("AF P_value"= AF_pval,
              "AF_perm_p" = AF_perm_pval,
              "p_vector" = c(AF_pval,AF_perm_pval)
  ))
}



####### COMBINE TWO ONE-SIDE P-VALUES

AdaptiveFisher_combo <- function(pu, pl, log=FALSE){
  
  AF_lower <- AdaptiveFisher_compute(pl,log = log)
  AF_upper <- AdaptiveFisher_compute(pu, log = log)
  
  #p_combo <- rbind(c(AF_lower$`AF P_value`,AF_upper$`AF P_value`),cbind(AF_lower$AF_perm_p,AF_upper$AF_perm_p))
  p_combo <- cbind(AF_lower$p_vector,AF_upper$p_vector)
  AF_comb <- AdaptiveFisher_compute(p_combo)
  
  return(list("AF_pvalue" = AF_comb$"AF P_value","P_vector"=AF_comb$"p_vector"))
  
}



# ----------------PART 3: improve method --------------

##### PCA ####

Y_PCA <- function(Y){
  pca <- prcomp(Y,scale. = TRUE)
  Y.new <- pca$x
  return(Y.new)
}

## residuals of Y regress on Z
Y.resid <- function(Y,Z){
  k <- ncol(Y)
  r <- matrix(NA,ncol = k,nrow = nrow(Y))
  for (i in 1:ncol(Y)) {
    glm.null <- glm(Y[,i] ~ Z,family = "gaussian")
    r[,i] <- glm.null$residuals 
  }
  return(r)
}

# ---------------- PART 4: Adds on -------------

min_p <- function(p){
  n <- nrow(p)
  zz<-apply(p, 1, min)
  p.vector <- rank(zz,ties.method = "max")/n
  p.minp<-mean(zz<=zz[1])
  return(list("p.min" = p.minp,"p.vector"=p.vector))
}
