library(mvtnorm)
library(bigstatsr)
library(genio)

PATH <- "/Users/ritwizkamal/Reproducibility/Reproducibility_Level4/"
outPATH <- "/Users/ritwizkamal/Reproducibility/Reproducibility_Level4/sample_output/"

#================================================================================================

#Change here to reproduce different simulation configurations

#Call Saved Seed
myseed <- read.table('/Users/ritwizkamal/Reproducibility/Simul_Analysis_1/Additional/Seeds/seed_exp1.4.1.seed')
myseed <- myseed$V1
set.seed(myseed)
print(myseed)

#Generate New Seed
# myseed <- as.integer(Sys.time())
# set.seed(myseed)
# print(myseed)

exp_num <- 'exp1.4.1'

Naux1_train <- 10000
Naux2_train <- 10000
Naux3_train <- 10000
Naux4_train <- 10000
Ntar_train <- 1000
Ntar_test <- 1000
Ntar_val <- 1000

rho_matrix <- matrix(data = c(1, 0.8, 0.8, 0.8, 0.8,
                              0.8, 1, 0.8, 0.8, 0.8,
                              0.8, 0.8, 1, 0.8, 0.8,
                              0.8, 0.8, 0.8, 1, 0.8,
                              0.8, 0.8, 0.8, 0.8, 1), nrow = 5)
h2 = 0.7

#================================================================================================

pheno_generation_5Pops = function(N_aux1, N_aux2, N_aux3, N_aux4, N_tar_train, N_tar_test, N_tar_val, 
                                  Z_aux1, Z_aux2, Z_aux3, Z_aux4, Z_tar_train, Z_tar_test, Z_tar_val, 
                                  rho_matrix, h2){
  
  M <- dim(Z_tar_train)[2]  # number of SNPs
  
  beta_aux1 <- rep(0, M)  # coef for aux1
  beta_aux2 <- rep(0, M)  # coef for aux2
  beta_aux3 <- rep(0, M)  # coef for aux3
  beta_aux4 <- rep(0, M)  # coef for aux4
  beta_tar <- rep(0, M)  # coef for target
  
  b <- rmvnorm(M, sigma = rho_matrix, )
  
  beta_aux1 <- b[,1]
  beta_aux2 <- b[,2]
  beta_aux3 <- b[,3]
  beta_aux4 <- b[,4]
  beta_tar <- b[,5]
  
  # phenotype for aux 1
  pheno_aux1 <- ((as.vector(Z_aux1%*%beta_aux1)/sd(Z_aux1%*%beta_aux1))*sqrt(h2)) + (rnorm(N_aux1, 0, 1) * sqrt(1-h2))
  
  # phenotype for aux 2
  pheno_aux2 <- ((as.vector(Z_aux2%*%beta_aux2)/sd(Z_aux2%*%beta_aux2))*sqrt(h2)) + (rnorm(N_aux2, 0, 1) * sqrt(1-h2))
  
  # phenotype for aux 3
  pheno_aux3 <- ((as.vector(Z_aux3%*%beta_aux3)/sd(Z_aux3%*%beta_aux3))*sqrt(h2)) + (rnorm(N_aux3, 0, 1) * sqrt(1-h2))
  
  # phenotype for aux 4
  pheno_aux4 <- ((as.vector(Z_aux4%*%beta_aux4)/sd(Z_aux4%*%beta_aux4))*sqrt(h2)) + (rnorm(N_aux4, 0, 1) * sqrt(1-h2))
  
  # phenotype for target train
  pheno_tar_train <- ((as.vector(Z_tar_train%*%beta_tar)/sd(Z_tar_train%*%beta_tar))*sqrt(h2)) + (rnorm(N_tar_train, 0, 1) * sqrt(1-h2))
  
  # phenotype for target test
  pheno_tar_test <- ((as.vector(Z_tar_test%*%beta_tar)/sd(Z_tar_test%*%beta_tar))*sqrt(h2)) + (rnorm(N_tar_test, 0, 1) * sqrt(1-h2))
  
  # phenotype for target val
  pheno_tar_val <- ((as.vector(Z_tar_val%*%beta_tar)/sd(Z_tar_val%*%beta_tar))*sqrt(h2)) + (rnorm(N_tar_val, 0, 1) * sqrt(1-h2)) 
  
  return(list(pheno_aux1=pheno_aux1, pheno_aux2=pheno_aux2, pheno_aux3=pheno_aux3, pheno_aux4=pheno_aux4,
              pheno_tar_train=pheno_tar_train, pheno_tar_val=pheno_tar_val, pheno_tar_test=pheno_tar_test,
              beta_aux1=beta_aux1, beta_aux2=beta_aux2, beta_aux3=beta_aux3, beta_aux4=beta_aux4, 
              beta_tar=beta_tar))
}

#================================================================================================

aux1_train_file <- 'E0_EUR.geno'
aux1info_file <- 'eur_info.info'
aux1_full <- read.table(paste0(PATH,aux1_train_file), sep=" ")
aux1_samples <- sample(1:dim(aux1_full)[2], Naux1_train, replace=FALSE, prob=NULL)
aux1_train <- as.matrix(aux1_full[,aux1_samples])
aux1_train_raw <- aux1_train
aux1_train <- t(aux1_train)
aux1_train <- scale(aux1_train, center = TRUE, scale = TRUE)
rm(aux1_full)

aux2_train_file <- 'E0_EAS.geno'
aux2info_file <- 'eas_info.info'
aux2_full <- read.table(paste0(PATH,aux2_train_file), sep=" ")
aux2_samples <- sample(1:dim(aux2_full)[2], Naux2_train, replace=FALSE, prob=NULL)
aux2_train <- as.matrix(aux2_full[,aux2_samples])
aux2_train_raw <- aux2_train
aux2_train <- t(aux2_train)
aux2_train <- scale(aux2_train, center = TRUE, scale = TRUE)
rm(aux2_full)

aux3_train_file <- 'E0_AMR.geno'
aux3info_file <- 'amr_info.info'
aux3_full <- read.table(paste0(PATH,aux3_train_file), sep=" ")
aux3_samples <- sample(1:dim(aux3_full)[2], Naux3_train, replace=FALSE, prob=NULL)
aux3_train <- as.matrix(aux3_full[,aux3_samples])
aux3_train_raw <- aux3_train
aux3_train <- t(aux3_train)
aux3_train <- scale(aux3_train, center = TRUE, scale = TRUE)
rm(aux3_full)

aux4_train_file <- 'E0_AFR.geno'
aux4info_file <- 'afr_info.info'
aux4_full <- read.table(paste0(PATH,aux4_train_file), sep=" ")
aux4_samples <- sample(1:dim(aux4_full)[2], Naux4_train, replace=FALSE, prob=NULL)
aux4_train <- as.matrix(aux4_full[,aux4_samples])
aux4_train_raw <- aux4_train
aux4_train <- t(aux4_train)
aux4_train <- scale(aux4_train, center = TRUE, scale = TRUE)
rm(aux4_full)

sas_train_file <- 'E0_SAS.geno'
snplist_file <- 'all_5_pops_common_snplist.snps'
sasinfo_file <- 'sas_info.info'
sas_full <- read.table(paste0(PATH,sas_train_file), sep=" ")
sas_samples <- sample(1:dim(sas_full)[2], Ntar_train+Ntar_val+Ntar_test, replace=FALSE, prob=NULL)
sas_train_samples <- sas_samples[1:Ntar_train]
sas_test_samples <- sas_samples[(Ntar_train+1):(Ntar_train+Ntar_test)]
sas_val_samples <- sas_samples[(Ntar_train+Ntar_test+1):(Ntar_train+Ntar_test+Ntar_val)]
sas_train <- as.matrix(sas_full[,sas_train_samples])
sas_test <- as.matrix(sas_full[,sas_test_samples])
sas_val <- as.matrix(sas_full[,sas_val_samples])
sas_train_raw <- sas_train
sas_test_raw <- sas_test
sas_val_raw <- sas_val
sas_train <- t(sas_train)
sas_test <- t(sas_test)
sas_val <- t(sas_val)
sas_train <- scale(sas_train, center = TRUE, scale = TRUE)
sas_test <- scale(sas_test, center = TRUE, scale = TRUE)
sas_val <- scale(sas_val, center = TRUE, scale = TRUE)
rm(sas_full)
#==============================================================================

pheno <- pheno_generation_5Pops(Naux1_train, Naux2_train, Naux3_train, Naux4_train, Ntar_train, Ntar_test, Ntar_val, 
                                aux1_train, aux2_train, aux3_train, aux4_train, sas_train, sas_test, sas_val, 
                                rho_matrix, h2 = h2)

beta_aux1 <- pheno$beta_aux1
beta_aux2 <- pheno$beta_aux2
beta_aux3 <- pheno$beta_aux3
beta_aux4 <- pheno$beta_aux4
beta_tar <- pheno$beta_tar
pheno_aux1 <- pheno$pheno_aux1
pheno_aux2 <- pheno$pheno_aux2
pheno_aux3 <- pheno$pheno_aux3
pheno_aux4 <- pheno$pheno_aux4
pheno_tar_train <- pheno$pheno_tar_train
pheno_tar_test <- pheno$pheno_tar_test
pheno_tar_val <- pheno$pheno_tar_val

ss_aux1 <- read.table(paste0(PATH, aux1info_file), sep=" ", header = TRUE)
temp_ss_aux1 <- big_univLinReg(X = as_FBM(aux1_train), pheno_aux1)
names(temp_ss_aux1) <- c("beta", "se", "t")
ss_aux1$beta <- temp_ss_aux1$beta
ss_aux1$se <- temp_ss_aux1$se
ss_aux1$t <- temp_ss_aux1$t
ss_aux1$pvalue <- 2*pnorm(-abs(temp_ss_aux1$t))
ss_aux1$n <- Naux1_train
rm(temp_ss_aux1)

ss_aux2 <- read.table(paste0(PATH, aux2info_file), sep=" ", header = TRUE)
temp_ss_aux2 <- big_univLinReg(X = as_FBM(aux2_train), pheno_aux2)
names(temp_ss_aux2) <- c("beta", "se", "t")
ss_aux2$beta <- temp_ss_aux2$beta
ss_aux2$se <- temp_ss_aux2$se
ss_aux2$t <- temp_ss_aux2$t
ss_aux2$pvalue <- 2*pnorm(-abs(temp_ss_aux2$t))
ss_aux2$n <- Naux2_train
rm(temp_ss_aux2)

ss_aux3 <- read.table(paste0(PATH, aux3info_file), sep=" ", header = TRUE)
temp_ss_aux3 <- big_univLinReg(X = as_FBM(aux3_train), pheno_aux3)
names(temp_ss_aux3) <- c("beta", "se", "t")
ss_aux3$beta <- temp_ss_aux3$beta
ss_aux3$se <- temp_ss_aux3$se
ss_aux3$t <- temp_ss_aux3$t
ss_aux3$pvalue <- 2*pnorm(-abs(temp_ss_aux3$t))
ss_aux3$n <- Naux3_train
rm(temp_ss_aux3)

ss_aux4 <- read.table(paste0(PATH, aux4info_file), sep=" ", header = TRUE)
temp_ss_aux4 <- big_univLinReg(X = as_FBM(aux4_train), pheno_aux4)
names(temp_ss_aux4) <- c("beta", "se", "t")
ss_aux4$beta <- temp_ss_aux4$beta
ss_aux4$se <- temp_ss_aux4$se
ss_aux4$t <- temp_ss_aux4$t
ss_aux4$pvalue <- 2*pnorm(-abs(temp_ss_aux4$t))
ss_aux4$n <- Naux4_train
rm(temp_ss_aux4)

ss_tar <- read.table(paste0(PATH, sasinfo_file), sep=" ", header = TRUE)
temp_ss_tar <- big_univLinReg(X = as_FBM(sas_train), pheno_tar_train)
names(temp_ss_tar) <- c("beta", "se", "t")
ss_tar$beta <- temp_ss_tar$beta
ss_tar$se <- temp_ss_tar$se
ss_tar$t <- temp_ss_tar$t
ss_tar$pvalue <- 2*pnorm(-abs(temp_ss_tar$t))
ss_tar$n <- Ntar_train
rm(temp_ss_tar)

write.table(ss_tar, paste0(outPATH,'target',exp_num,'.sumstat'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(ss_aux1, paste0(outPATH,'aux_one',exp_num,'.sumstat'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(ss_aux2, paste0(outPATH,'aux_two',exp_num,'.sumstat'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(ss_aux3, paste0(outPATH,'aux_three',exp_num,'.sumstat'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(ss_aux4, paste0(outPATH,'aux_four',exp_num,'.sumstat'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)

write.table(beta_aux1, paste0(outPATH,'beta_aux_one',exp_num,'.truebetas'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(beta_aux2, paste0(outPATH,'beta_aux_two',exp_num,'.truebetas'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(beta_aux3, paste0(outPATH,'beta_aux_three',exp_num,'.truebetas'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(beta_aux4, paste0(outPATH,'beta_aux_four',exp_num,'.truebetas'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(beta_tar, paste0(outPATH,'beta_tar',exp_num,'.truebetas'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)

write.table(myseed, paste0(outPATH,'seed',exp_num,'.seed'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)

write.table(pheno_aux1, paste0(outPATH,'pheno_aux_one',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(pheno_aux2, paste0(outPATH,'pheno_aux_two',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(pheno_aux3, paste0(outPATH,'pheno_aux_three',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(pheno_aux4, paste0(outPATH,'pheno_aux_four',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)

write.table(pheno_tar_train, paste0(outPATH,'pheno_tar_train',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(pheno_tar_test, paste0(outPATH,'pheno_tar_test',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(pheno_tar_val, paste0(outPATH,'pheno_tar_val',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)

write_bed(paste0(outPATH,'Target_train',exp_num,'.bed'), sas_train_raw)
write_bed(paste0(outPATH,'Target_test',exp_num,'.bed'), sas_test_raw)
write_bed(paste0(outPATH,'Target_val',exp_num,'.bed'), sas_val_raw)

write.table(sas_train, paste0(outPATH,'Target_train',exp_num,'.geno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(sas_test, paste0(outPATH,'Target_test',exp_num,'.geno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(sas_val, paste0(outPATH,'Target_val',exp_num,'.geno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)

