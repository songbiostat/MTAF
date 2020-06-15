### power,mix

args <- commandArgs(TRUE)
core.i <- as.numeric(args)


# load packages, functions
path="/users/PAS1149/osu10039/PAS1149/qldeng/AF_data/june_simulations"
setwd(paste0(path,"/","code","/","functions"))
source("af_functions_mix.R")

t <- 50
rho = 0.6
cor_matrix = "cs"
n=1000


rep <- 50  # 1000/n. cores

# matrix/vectors save results
pval_AF_comb <- matrix(NA,nrow = rep,ncol = 5)
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
  #beta <- c(runif(4,0.05,0.3),rep(0,6))
  beta <- c(runif(10,0.05,0.25),rep(0,40))
  
  sim <- SIM_AF(beta_z=0.5,cor=rho,beta=beta,traits=t,mix=TRUE)
  
  ## Score test
  col_bin <- seq(1,t,2)
  Y_cont <- sim$Y[,-col_bin]
  Y_bin <- sim$Y[,col_bin]
  
  n_perm=1e3
  X_perm <- cbind(sim$X, replicate(n_perm, sample(sim$X)))
  ######################  AF Bin ######################
  # bin
  score_bin <- score_glm_bin(Y_bin,sim$X,cov=sim$Z,overall=TRUE,X_permuted=X_perm)
  # AF of bin
  AF_bin <- AdaptiveFisher_combo(score_bin$p_l,score_bin$p_u)
  
  ######################  AF Cont ######################
  # Regress Y on Z and do PCA
  Y.res <- resid(lm(Y_cont~sim$Z))
  Y.pca <- Y_PCA(Y.res)
  
  ## Score test
  score <- score_glm(cbind(Y.pca,Y.res),sim$X,cov = sim$Z,nperm = 1e3,overall=TRUE,X_permuted=X_perm)
  t=t-length(col_bin)
  ###################### Original AF  ######################
  ## AF w./ PCA
  PCA_res <- AdaptiveFisher_combo(score$p_l[,1:t],score$p_u[,1:t])
  p_perm_pcares <- matrix(PCA_res$P_vector,ncol = 1)
  # AF w.o. PCA
  AF <- AdaptiveFisher_combo(score$p_l[,-(1:t)],score$p_u[,-(1:t)])
  p_perm_af <- matrix(AF$P_vector,ncol = 1)
  # combine 
  AF_cont <- AdaptiveFisher_compute(cbind(p_perm_af,p_perm_pcares))
  
  
  # combine 
  p_perm_bin <- matrix(AF_bin$P_vector,ncol = 1)
  p_perm_cont <- matrix(AF_cont$p_vector,ncol = 1)
  combo <- AdaptiveFisher_compute(cbind(p_perm_bin,p_perm_cont))
  pval_AF_comb[i,1] <- combo$`AF P_value`
  pval_AF_comb[i,2] <- AF_bin$AF_pvalue
  pval_AF_comb[i,3] <- AF_cont$`AF P_value`
  pval_AF_comb[i,4] <- AF$AF_pvalue
  pval_AF_comb[i,5] <- PCA_res$AF_pvalue
  
  
  t=t+length(col_bin)
  score_all <- matrix(rep(0,t*(n+1)),ncol = t)
  score_all[,col_bin] <- score_bin$p_l
  t=t-length(col_bin)
  score_all[,-col_bin] <- score$p_l[,-(1:t)]
  t=t+length(col_bin)
  #################### Other Methods #######################
  
  ### minP ###
  p_matrix_2side <- 2*pnorm(abs(qnorm(score_all)),lower.tail = FALSE)
  pval_minP[i] <- min_p(p_matrix_2side)$p.min
  
  ### TATES ###
  TATEs <- TATES(sim$Y,score_all,t,n.snp = 1)
  pval_TATES[i] <- TATEs[t+1]
  
  #### MultiPhen ### 
  X = as.matrix(sim$X)
  dimnames(X) <- list(1:n,1)
  Y = as.matrix(sim$Y)
  dimnames(Y) <- list(1:n,1:t)
  opts = mPhen.options(c("regression","pheno.input"))
  invisible(capture.output(m <- mPhen(X[,1,drop=FALSE],Y,phenotypes = dimnames(Y)[[2]][1:t],covariates = c("z1","z2"),opts = opts))) 
  pval_MultiPhen[i] <- m$Results[,,,"pvalue"] ["JointModel"]  
  
}

print(paste("power_mix",t,rho,sep = "_"))
colMeans(pval_AF_comb<0.05,na.rm = TRUE)
mean(pval_MultiPhen<0.05,na.rm = TRUE)
mean(pval_TATES<0.05,na.rm = TRUE)
mean(pval_minP<0.05,na.rm = TRUE)

setwd(paste0(path,"/","power/mix/trait",t))

save(list = ls(all.names = TRUE), file = paste0("power_mix","_",t,"_",rho,"_",core.i,".RData"), envir = .GlobalEnv)

