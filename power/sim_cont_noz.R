

### power_cont_noz

args <- commandArgs(TRUE)
core.i <- as.numeric(args)


# load packages, functions
path="/users/PAS1149/osu10039/PAS1149/qldeng/AF_data/june_simulations"
setwd(paste0(path,"/","code","/","functions"))
source("af_functions_noz.R")


n = 1000 # sample size
betaz = 0 # effect of covariates

t <- 100 # traits
cor_matrix = "cs"  # correlation structure
rho = 0.6  # correlation strength
sparsity="sparse"

rep <- 50   # permutation for each core


# matrix/vectors save results
pval_AF_comb <- matrix(NA,nrow = rep,ncol = 3)
pval_minP <- rep(NA,rep)
pval_MANOVA <- rep(NA,rep)
pval_SPU <- matrix(NA,nrow = rep,ncol = 10)
pval_SPU_ex <- matrix(NA,nrow = rep,ncol = 10)
pval_MSKAT <- matrix(NA,nrow = rep,ncol = 2)
pval_TATES <- rep(NA,rep)
pval_MultiPhen <- rep(NA,rep)



for(i in (1:rep)){
  print(i)
  set.seed(core.i*1000 + i)  # random seed
  
  # beta
  #beta <- c(runif(1,0.15,0.25),rep(0,9))
  #beta <- c(runif(4,0.05,0.15),rep(0,6))
  #beta <- c(runif(1,0.2,0.4),rep(0,49))
  #beta <- c(runif(10,0.05,0.12),rep(0,40))
  beta <- c(runif(2,0.15,0.3),rep(0,98))
  #beta <- c(runif(20,0.02,0.1),rep(0,80))
  
  sim <- SIM_covariate_power_Z(betaz=betaz,constr = cor_matrix,beta,traits = t,cat_ratio =0,cor=rho,n.cov = 2, MAF = 0.3)
  
  # do PCA
  Y.res <- sim$Y
  X.res <- sim$X
  Y.pca <- Y_PCA(Y.res)
  
  ## Score test
  score <- score_test(cbind(Y.pca,Y.res),X.res,n)
    
  ###################### Original AF  ######################
  ## AF w./ PCA
  PCA_res <- AdaptiveFisher_combo(score$p_l[,1:t],score$p_u[,1:t])
  p_perm_pcares <- matrix(PCA_res$P_vector,ncol = 1)
  # AF w.o. PCA
  AF <- AdaptiveFisher_combo(score$p_l[,-(1:t)],score$p_u[,-(1:t)])
  p_perm_af <- matrix(AF$P_vector,ncol = 1)
  # combine 
  combo <- AdaptiveFisher_compute(cbind(p_perm_af,p_perm_pcares))
  pval_AF_comb[i,1] <- combo$`AF P_value`
  # record each method's p-value and min
  pval_AF_comb[i,2] <- PCA_res$AF_pvalue
  pval_AF_comb[i,3] <- AF$AF_pvalue
  
  
  
  #################### Other Methods #######################
  ### minP ###
  p_matrix_2side <- 2*pnorm(abs(qnorm(score$p_l[,-(1:t)])),lower.tail = FALSE)
  pval_minP[i] <- min_p(p_matrix_2side)$p.min
  
  ### SPUs & aSPU ###
  invisible(capture.output(suppressMessages(pan<-GEEaSPU(traits = sim$Y,geno = sim$X,model = "gaussian"))))
  pval_SPU[i,] <- pan[1:10]
  
  invisible(capture.output(suppressMessages(pan<-GEEaSPU(traits = sim$Y,geno = sim$X,model = "gaussian",corstr = "exchangeable"))))
  pval_SPU_ex[i,] <- pan[1:10] 
  

  ###  MSKAT  ###
  X_matrix <- matrix(sim$X,ncol = 1)
  SKAT <- MSKAT(MSKAT.cnull(sim$Y,X=NULL), X_matrix, W.beta=c(1,25))
  pval_MSKAT[i,] <- SKAT
  
  #### MANOVA ####
  man <- manova(sim$Y~sim$X)
  pval_MANOVA[i] <- summary(man)$stats["sim$X","Pr(>F)"]
  
  ### TATES ###
  TATEs <- TATES(sim$Y,score$p_l[,-(1:t)],t,n.snp = 1)
  pval_TATES[i] <- TATEs[t+1]
  
  #### MultiPhen ### 
  X = as.matrix(sim$X)
  dimnames(X) <- list(1:n,1)
  Y = as.matrix(sim$Y)
  dimnames(Y) <- list(1:n,1:t)
  opts = mPhen.options(c("regression","pheno.input"))
  opts$mPhen.scoreTest = TRUE
  invisible(capture.output(suppressMessages(m <- mPhen(X[,1,drop=FALSE],Y,phenotypes = "all", opts = opts)))) 
  pval_MultiPhen[i] <- m$Results[,,,"pvalue"] ["JointModel"]
}


print(paste("power_cont_noz",t,rho,sep = "_"))
colMeans(pval_AF_comb<0.05,na.rm = TRUE)
colMeans(pval_MSKAT<0.05,na.rm = TRUE)
colMeans(pval_SPU<0.05,na.rm = TRUE)
colMeans(pval_SPU_ex<0.05,na.rm = TRUE)
mean(pval_MultiPhen<0.05,na.rm = TRUE)
mean(pval_TATES<0.05,na.rm = TRUE)
mean(pval_MANOVA<0.05,na.rm = TRUE)
mean(pval_minP<0.05,na.rm = TRUE)


setwd(paste0(path,"/","power/continuous/no_covariates/trait",t,"/",sparsity))
save(list = ls(all.names = TRUE), file = paste0("power_cont_noz","_",t,"_",sparsity,"_",rho,"_",core.i,".RData"), envir = .GlobalEnv)

