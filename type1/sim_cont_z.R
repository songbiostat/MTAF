### type 1 error

args <- commandArgs(TRUE)
core.i <- as.numeric(args)


# load packages, functions
path="/fs/project/PAS1149/qldeng/AF_data/june_simulations/"
setwd(paste0(path,"/","code","/","functions"))
source("af_functions_z.R")




n = 1000 # sample size
t <-100  # traits
cor_matrix = "cs"
rho = 0.6
betaZ = 0.5  # effect of covariates

rep <- 50  # permutations per core

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
  set.seed(core.i*10000 + i)  # random seed
  
  # type1: null beta
  beta <- rep(0,t)

  sim <- SIM_covariate_power_Z(betaz=betaZ,constr = cor_matrix,beta,traits = t,cat_ratio = 0,cor=rho,n.cov = 2, MAF = 0.3)
  
  # Regress Y on Z and do PCA
  Y.res <- resid(lm(sim$Y~sim$Z))
  X.res <- resid(lm(sim$X~sim$Z))
  Y.pca <- Y_PCA(Y.res)
  
  ## Score test
  score <- score_glm(cbind(Y.pca,Y.res),sim$X,cov = sim$Z,nperm = 1e3)
  
   
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
  invisible(capture.output(suppressMessages(pan1<-GEEaSPU(traits = sim$Y,geno = sim$X,Z=sim$Z,model = "gaussian"))))
  pval_SPU[i,] <- pan1[1:10]
  
  invisible(capture.output(suppressMessages(pan2<-GEEaSPU(traits = sim$Y,geno = sim$X,Z=sim$Z,model = "gaussian",corstr = "exchangeable"))))
  pval_SPU_ex[i,] <- pan2[1:10] 
  

  ###  MSKAT  ###
  X_matrix <- matrix(sim$X,ncol = 1)
  SKAT <- MSKAT(MSKAT.cnull(sim$Y,sim$Z), X_matrix, W.beta=c(1,25))
  pval_MSKAT[i,] <- SKAT
  
  #### MANOVA ####
  man <- manova(sim$Y~sim$X+sim$Z)
  pval_MANOVA[i] <- summary(man)$stats["sim$X","Pr(>F)"]
  
  ### TATES ###
  TATEs <- TATES(sim$Y,score$p_l[,-(1:t)],t,n.snp = 1)
  pval_TATES[i] <- TATEs[t+1]
  
  #### MultiPhen ### 
  X = as.matrix(sim$X)
  dimnames(X) <- list(1:n,1)
  Y = as.matrix(cbind(sim$Y,sim$Z))
  dimnames(Y) <- list(1:n,c(1:t,"z1","z2"))

  opts = mPhen.options(c("regression","pheno.input"))
  invisible(capture.output(m <- mPhen(X[,1,drop=FALSE],Y,phenotypes = dimnames(Y)[[2]][1:t],covariates = c("z1","z2")))) 
  pval_MultiPhen[i] <- m$Results[,,,"pvalue"] ["JointModel"]
}

print(paste("type1_cont_z_",t,rho,sep = "_"))
colMeans(pval_AF_comb<0.05,na.rm = TRUE)
colMeans(pval_MSKAT<0.05,na.rm = TRUE)
colMeans(pval_SPU<0.05,na.rm = TRUE)
colMeans(pval_SPU_ex<0.05,na.rm = TRUE)
mean(pval_MultiPhen<0.05,na.rm = TRUE)
mean(pval_TATES<0.05,na.rm = TRUE)
mean(pval_MANOVA<0.05,na.rm = TRUE)
mean(pval_minP<0.05,na.rm = TRUE)

setwd(paste0(path,"/","type1/continuous/covariates/trait",t))
save(list = ls(all.names = TRUE), file = paste0("type1_f_cont_z","_",t,"_",rho,"_",core.i,".RData"), envir = .GlobalEnv)




