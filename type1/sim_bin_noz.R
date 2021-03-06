### type 1 error

args <- commandArgs(TRUE)
core.i <- as.numeric(args)

# load packages, functions
path="/fs/project/PAS1149/qldeng/AF_data/june_simulations"
setwd(paste0(path,"/","code","/","functions"))
source("af_functions_noz.R")


n=1000
t <- 50
cor_matrix = "cs"
rho = 0.6

rep <- 50    # 1000/n. cores

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
  set.seed(core.i*30000 + i)  # random seed
  
  # type1: null beta
  beta <- rep(0,t)

  sim <- SIM_binary_power_noZ(constr = cor_matrix, beta=beta,traits = t,cor=rho,n.cov = 2, MAF = 0.3)
  
  # 
  Y.res <- sim$Y
  X.res <- sim$X
  
  ## Score test
  score <- score_test(Y.res,X.res,n, binary = TRUE)
  #score <- score_glm_bin(Y.res,X.res,overall = FALSE)
  
  ###################### Original AF  ######################
  # AF w.o. PCA
  AF <- AdaptiveFisher_combo(score$p_l,score$p_u)
  pval_AF_comb[i,1] <- AF$AF_pvalue
  
  
  #################### Other Methods #######################
  ### minP ###
  p_matrix_2side <- 2*pnorm(abs(qnorm(score$p_l)),lower.tail = FALSE)
  pval_minP[i] <- min_p(p_matrix_2side)$p.min
  
  ### SPUs & aSPU ###
 
  invisible(capture.output(suppressMessages(pan1<-GEEaSPU(sim$Y,geno = sim$X,model = "binomial"))))
  pval_SPU[i,] <- pan1[1:10]
  invisible(capture.output(suppressMessages(pan2<-GEEaSPU(sim$Y,geno = sim$X,model = "binomial",corstr = "exchangeable"))))
  pval_SPU_ex[i,] <- pan2[1:10]

  
  ### TATES ###
  TATEs <- TATES(sim$Y,score$p_l,t,n.snp = 1)
  pval_TATES[i] <- TATEs[t+1]
  
  #### MultiPhen ### 
  X = as.matrix(sim$X)
  dimnames(X) <- list(1:n,1)
  Y = as.matrix(sim$Y)
  dimnames(Y) <- list(1:n,1:t)
  opts = mPhen.options(c("regression","pheno.input"))
  invisible(capture.output(suppressMessages(m <- mPhen(X[,1,drop=FALSE],Y,phenotypes = "all", opts = opts)))) 
  pval_MultiPhen[i] <- m$Results[,,,"pvalue"] ["JointModel"]  
}

print(paste("type1_cont_z_",t,rho,sep = "_"))
colMeans(pval_AF_comb<0.05,na.rm = TRUE)
colMeans(pval_SPU<0.05,na.rm = TRUE)
colMeans(pval_SPU_ex<0.05,na.rm = TRUE)
mean(pval_MultiPhen<0.05,na.rm = TRUE)
mean(pval_TATES<0.05,na.rm = TRUE)
mean(pval_minP<0.05,na.rm = TRUE)

setwd(paste0(path,"/","type1/binary/no_covariates/trait",t))

save(list = ls(all.names = TRUE), file = paste0("type1_f_bin_noz","_",t,"_",rho,"_",core.i,".RData"), envir = .GlobalEnv)




