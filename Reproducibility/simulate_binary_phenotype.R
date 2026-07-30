library(mvtnorm)
library(bigstatsr)
library(genio)

PATH <- "<INPUT PATH>"
outPATH <- "<OUTPUT PATH>"

#================================================================================================

#Generate New Seed
myseed <- as.integer(Sys.time())
set.seed(myseed)
print(myseed)

exp_num <- 'exp1.1.1'

Naux1_train <- 1000
Naux2_train <- 1000
Naux3_train <- 1000
Naux4_train <- 1000
Ntar_train <- 1000
Ntar_test <- 1000
Ntar_val <- 1000

rho_matrix <- matrix(data = c(1, 0.8, 0.8, 0.8, 0.8,
                              0.8, 1, 0.8, 0.8, 0.8,
                              0.8, 0.8, 1, 0.8, 0.8,
                              0.8, 0.8, 0.8, 1, 0.8,
                              0.8, 0.8, 0.8, 0.8, 1), nrow = 5)
h2 = 0.7
prevalence = 0.2

#================================================================================================

impute_geno = function(in_matrix){
  zero_sd_rows <- which(apply(in_matrix, 1, sd) == 0)
  for (row in zero_sd_rows)
  {
    temp_row <- in_matrix[row,]
    temp_sample <- sample(1:length(temp_row), 1)
    val <- temp_row[temp_sample]
    temp_row[temp_sample] <- ifelse(val == 0, 1,
                                    ifelse(val == 1, 0,
                                           ifelse(val == 2, 1, val)))
    in_matrix[row,] <- temp_row
  }
  
  return(in_matrix)
}

pheno_generation_binary_5Pops <- function(N_aux1, N_aux2, N_aux3, N_aux4, N_tar_train, N_tar_test, N_tar_val, 
                                          Z_aux1, Z_aux2, Z_aux3, Z_aux4, Z_tar_train, Z_tar_test, Z_tar_val, 
                                          rho_matrix, h2, prevalence) {
  
  M <- dim(Z_tar_train)[2]  # number of SNPs
  
  beta_aux1 <- rep(0, M)
  beta_aux2 <- rep(0, M)
  beta_aux3 <- rep(0, M)
  beta_aux4 <- rep(0, M)
  beta_tar <- rep(0, M)
  
  b <- rmvnorm(M, sigma = rho_matrix)
  
  beta_aux1 <- b[,1]
  beta_aux2 <- b[,2]
  beta_aux3 <- b[,3]
  beta_aux4 <- b[,4]
  beta_tar <- b[,5]
  
  # --- 2. Calculate Continuous Liability  ---
  
  # liability for aux 1
  lia_aux1 <- ((as.vector(Z_aux1 %*% beta_aux1) / sd(Z_aux1 %*% beta_aux1)) * sqrt(h2)) + (rnorm(N_aux1, 0, 1) * sqrt(1-h2))
  
  # liability for aux 2
  lia_aux2 <- ((as.vector(Z_aux2 %*% beta_aux2) / sd(Z_aux2 %*% beta_aux2)) * sqrt(h2)) + (rnorm(N_aux2, 0, 1) * sqrt(1-h2))
  
  # liability for aux 3
  lia_aux3 <- ((as.vector(Z_aux3 %*% beta_aux3) / sd(Z_aux3 %*% beta_aux3)) * sqrt(h2)) + (rnorm(N_aux3, 0, 1) * sqrt(1-h2))
  
  # liability for aux 4
  lia_aux4 <- ((as.vector(Z_aux4 %*% beta_aux4) / sd(Z_aux4 %*% beta_aux4)) * sqrt(h2)) + (rnorm(N_aux4, 0, 1) * sqrt(1-h2))
  
  # liability for target train
  lia_tar_train <- ((as.vector(Z_tar_train %*% beta_tar) / sd(Z_tar_train %*% beta_tar)) * sqrt(h2)) + (rnorm(N_tar_train, 0, 1) * sqrt(1-h2))
  
  # liability for target test
  lia_tar_test <- ((as.vector(Z_tar_test %*% beta_tar) / sd(Z_tar_test %*% beta_tar)) * sqrt(h2)) + (rnorm(N_tar_test, 0, 1) * sqrt(1-h2))
  
  # liability for target val
  lia_tar_val <- ((as.vector(Z_tar_val %*% beta_tar) / sd(Z_tar_val %*% beta_tar)) * sqrt(h2)) + (rnorm(N_tar_val, 0, 1) * sqrt(1-h2))
  
  # --- 3. NEW: Apply Threshold to get Binary Phenotypes ---
  
  # Calculate the threshold from the standard normal distribution
  # e.g., for prevalence = 0.1, threshold is qnorm(0.9) approx 1.28
  threshold <- qnorm(prevalence, lower.tail = FALSE)
  
  # Convert liabilities to binary 0/1 (Control/Case)
  pheno_aux1 <- as.numeric(lia_aux1 > threshold)
  pheno_aux2 <- as.numeric(lia_aux2 > threshold)
  pheno_aux3 <- as.numeric(lia_aux3 > threshold)
  pheno_aux4 <- as.numeric(lia_aux4 > threshold)
  pheno_tar_train <- as.numeric(lia_tar_train > threshold)
  pheno_tar_test <- as.numeric(lia_tar_test > threshold)
  pheno_tar_val <- as.numeric(lia_tar_val > threshold)
  
  # --- 4. Return All Results ---
  
  return(list(
    # Binary Phenotypes
    pheno_aux1 = pheno_aux1, pheno_aux2 = pheno_aux2, pheno_aux3 = pheno_aux3, pheno_aux4 = pheno_aux4,
    pheno_tar_train = pheno_tar_train, pheno_tar_val = pheno_tar_val, pheno_tar_test = pheno_tar_test,
    
    # Underlying Liabilities
    lia_aux1 = lia_aux1, lia_aux2 = lia_aux2, lia_aux3 = lia_aux3, lia_aux4 = lia_aux4,
    lia_tar_train = lia_tar_train, lia_tar_val = lia_tar_val, lia_tar_test = lia_tar_test,
    
    # Betas
    beta_aux1 = beta_aux1, beta_aux2 = beta_aux2, beta_aux3 = beta_aux3, beta_aux4 = beta_aux4, 
    beta_tar = beta_tar
  ))
}

#================================================================================================

aux1_train_file <- 'E0_EUR.geno'
aux1_info_file_path <- 'eur_info.info'
aux1_full <- read.table(paste0(PATH,aux1_train_file), sep=" ") #This is a SNPs x Samples matrix
aux1_samples <- sample(1:dim(aux1_full)[2], Naux1_train, replace=FALSE, prob=NULL)
aux1_train <- as.matrix(aux1_full[,aux1_samples])
if(length(which(apply(aux1_train, 1, sd) == 0)) > 0){
  aux1_train <- impute_geno(aux1_train)
}
aux1_train_raw <- aux1_train
aux1_train <- t(aux1_train)
aux1_train <- scale(aux1_train, center = TRUE, scale = TRUE)
rm(aux1_full)

aux2_train_file <- 'E0_EAS.geno'
aux2_info_file_path <- 'eas_info.info'
aux2_full <- read.table(paste0(PATH,aux2_train_file), sep=" ") #This is a SNPs x Samples matrix
aux2_samples <- sample(1:dim(aux2_full)[2], Naux2_train, replace=FALSE, prob=NULL)
aux2_train <- as.matrix(aux2_full[,aux2_samples])
if(length(which(apply(aux2_train, 1, sd) == 0)) > 0){
  aux2_train <- impute_geno(aux2_train)
}
aux2_train_raw <- aux2_train
aux2_train <- t(aux2_train)
aux2_train <- scale(aux2_train, center = TRUE, scale = TRUE)
rm(aux2_full)

aux3_train_file <- 'E0_AMR.geno'
aux3_info_file_path <- 'amr_info.info'
aux3_full <- read.table(paste0(PATH,aux3_train_file), sep=" ") #This is a SNPs x Samples matrix
aux3_samples <- sample(1:dim(aux3_full)[2], Naux3_train, replace=FALSE, prob=NULL)
aux3_train <- as.matrix(aux3_full[,aux3_samples])
if(length(which(apply(aux3_train, 1, sd) == 0)) > 0){
  aux3_train <- impute_geno(aux3_train)
}
aux3_train_raw <- aux3_train
aux3_train <- t(aux3_train)
aux3_train <- scale(aux3_train, center = TRUE, scale = TRUE)
rm(aux3_full)

aux4_train_file <- 'E0_AFR.geno'
aux4_info_file_path <- 'afr_info.info'
aux4_full <- read.table(paste0(PATH,aux4_train_file), sep=" ") #This is a SNPs x Samples matrix
aux4_samples <- sample(1:dim(aux4_full)[2], Naux4_train, replace=FALSE, prob=NULL)
aux4_train <- as.matrix(aux4_full[,aux4_samples])
if(length(which(apply(aux4_train, 1, sd) == 0)) > 0){
  aux4_train <- impute_geno(aux4_train)
}
aux4_train_raw <- aux4_train
aux4_train <- t(aux4_train)
aux4_train <- scale(aux4_train, center = TRUE, scale = TRUE)
rm(aux4_full)


tar_train_file <- 'E0_SAS.geno'
tar_info_file_path <- 'sas_info.info'
tar_full <- read.table(paste0(PATH,tar_train_file), sep=" ") #This is a SNPs x Samples matrix
tar_all_samples <- sample(1:dim(tar_full)[2], Ntar_train+Ntar_val+Ntar_test, replace=FALSE, prob=NULL)
tar_train_samples <- tar_all_samples[1:Ntar_train]
tar_val_samples <- tar_all_samples[(Ntar_train+1):(Ntar_train+Ntar_val)]
tar_test_samples <- tar_all_samples[(Ntar_train+Ntar_val+1):(Ntar_train+Ntar_val+Ntar_test)]
tar_train <- as.matrix(tar_full[,tar_train_samples])
tar_test <- as.matrix(tar_full[,tar_test_samples])
tar_val <- as.matrix(tar_full[,tar_val_samples])
if(length(which(apply(tar_train, 1, sd) == 0)) > 0){
  tar_train <- impute_geno(tar_train)
}
if(length(which(apply(tar_val, 1, sd) == 0)) > 0){
  tar_val <- impute_geno(tar_val)
}
if(length(which(apply(tar_test, 1, sd) == 0)) > 0){
  tar_test <- impute_geno(tar_test)
}
tar_train_raw <- tar_train
tar_test_raw <- tar_test
tar_val_raw <- tar_val
tar_train <- t(tar_train)
tar_test <- t(tar_test)
tar_val <- t(tar_val)
tar_train <- scale(tar_train, center = TRUE, scale = TRUE)
tar_test <- scale(tar_test, center = TRUE, scale = TRUE)
tar_val <- scale(tar_val, center = TRUE, scale = TRUE)
rm(tar_full)
#==============================================================================

pheno <- pheno_generation_binary_5Pops(Naux1_train, Naux2_train, Naux3_train, Naux4_train, Ntar_train, 
                                       Ntar_test, Ntar_val, aux1_train, aux2_train, aux3_train, aux4_train, 
                                       tar_train, tar_test, tar_val, rho_matrix, h2 = h2, prevalence=prevalence)

beta_aux1 <- pheno$beta_aux1
beta_aux2 <- pheno$beta_aux2
beta_aux3 <- pheno$beta_aux3
beta_aux4 <- pheno$beta_aux4
beta_tar <- pheno$beta_tar

lia_aux1 <- pheno$lia_aux1
lia_aux2 <- pheno$lia_aux2
lia_aux3 <- pheno$lia_aux3
lia_aux4 <- pheno$lia_aux4
lia_tar_train <- pheno$lia_tar_train
lia_tar_test <- pheno$lia_tar_test
lia_tar_val <- pheno$lia_tar_val

pheno_aux1 <- pheno$pheno_aux1
pheno_aux2 <- pheno$pheno_aux2
pheno_aux3 <- pheno$pheno_aux3
pheno_aux4 <- pheno$pheno_aux4
pheno_tar_train <- pheno$pheno_tar_train
pheno_tar_test <- pheno$pheno_tar_test
pheno_tar_val <- pheno$pheno_tar_val

# Write Ground Truth Betas
write.table(beta_aux1, paste0(outPATH,'beta_aux_one_',exp_num,'.truebetas'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(beta_aux2, paste0(outPATH,'beta_aux_two_',exp_num,'.truebetas'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(beta_aux3, paste0(outPATH,'beta_aux_three_',exp_num,'.truebetas'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(beta_aux4, paste0(outPATH,'beta_aux_four_',exp_num,'.truebetas'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(beta_tar, paste0(outPATH,'beta_tar_',exp_num,'.truebetas'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)

# Write Seed for this simulation
write.table(myseed, paste0(outPATH,'seed_',exp_num,'.seed'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)

# Write Ground Truth Phenotypes
write.table(pheno_aux1, paste0(outPATH,'pheno_aux_one_',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(pheno_aux2, paste0(outPATH,'pheno_aux_two_',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(pheno_aux3, paste0(outPATH,'pheno_aux_three_',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(pheno_aux4, paste0(outPATH,'pheno_aux_four_',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)

write.table(pheno_tar_train, paste0(outPATH,'pheno_tar_train_',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(pheno_tar_test, paste0(outPATH,'pheno_tar_test_',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(pheno_tar_val, paste0(outPATH,'pheno_tar_val_',exp_num,'.truepheno'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)

# Write Ground Truth Liabilities
write.table(lia_aux1, paste0(outPATH,'lia_aux_one_',exp_num,'.trueliability'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(lia_aux2, paste0(outPATH,'lia_aux_two_',exp_num,'.trueliability'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(lia_aux3, paste0(outPATH,'lia_aux_three_',exp_num,'.trueliability'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(lia_aux4, paste0(outPATH,'lia_aux_four_',exp_num,'.trueliability'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)

write.table(lia_tar_train, paste0(outPATH,'lia_tar_train_',exp_num,'.trueliability'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(lia_tar_test, paste0(outPATH,'lia_tar_test_',exp_num,'.trueliability'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(lia_tar_val, paste0(outPATH,'lia_tar_val_',exp_num,'.trueliability'), sep=" ", row.names = FALSE, col.names = FALSE, quote=FALSE)

# Write Genotypes

#Aux 1
write_bed(paste0(outPATH,'aux1_train_',exp_num,'.bed'), aux1_train_raw) # bed file
aux1info_file <- read.table(paste0(PATH, aux1_info_file_path), sep=" ", header=TRUE)
aux1_bim_data <- as.data.frame(aux1info_file$chr)
aux1_bim_data$snp <- aux1info_file$rsid
aux1_bim_data$zeros <- rep(0, dim(aux1_train_raw)[1])
aux1_bim_data$pos <- aux1info_file$bp
aux1_bim_data$eff <- aux1info_file$effAllele
aux1_bim_data$ref <- aux1info_file$refAllele
write.table(aux1_bim_data, paste0(outPATH,'aux1_train_',exp_num,'.bim'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #bim file
aux1_fam_data <- as.data.frame(1:dim(aux1_train_raw)[2])
aux1_fam_data$col2 <- aux1_fam_data[,1]
aux1_fam_data$col3 <- rep(0, length(aux1_fam_data[,1]))
aux1_fam_data$col4 <- rep(0, length(aux1_fam_data[,1]))
aux1_fam_data$col5 <- rep(0, length(aux1_fam_data[,1]))
aux1_fam_data$col6 <- rep(-9, length(aux1_fam_data[,1]))
write.table(aux1_fam_data, paste0(outPATH,'aux1_train_',exp_num,'.fam'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #fam file

#Aux 2
write_bed(paste0(outPATH,'aux2_train_',exp_num,'.bed'), aux2_train_raw) # bed file
aux2info_file <- read.table(paste0(PATH, aux2_info_file_path), sep=" ", header=TRUE)
aux2_bim_data <- as.data.frame(aux2info_file$chr)
aux2_bim_data$snp <- aux2info_file$rsid
aux2_bim_data$zeros <- rep(0, dim(aux2_train_raw)[1])
aux2_bim_data$pos <- aux2info_file$bp
aux2_bim_data$eff <- aux2info_file$effAllele
aux2_bim_data$ref <- aux2info_file$refAllele
write.table(aux2_bim_data, paste0(outPATH,'aux2_train_',exp_num,'.bim'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #bim file
aux2_fam_data <- as.data.frame(1:dim(aux2_train_raw)[2])
aux2_fam_data$col2 <- aux2_fam_data[,1]
aux2_fam_data$col3 <- rep(0, length(aux2_fam_data[,1]))
aux2_fam_data$col4 <- rep(0, length(aux2_fam_data[,1]))
aux2_fam_data$col5 <- rep(0, length(aux2_fam_data[,1]))
aux2_fam_data$col6 <- rep(-9, length(aux2_fam_data[,1]))
write.table(aux2_fam_data, paste0(outPATH,'aux2_train_',exp_num,'.fam'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #fam file

#Aux 3
write_bed(paste0(outPATH,'aux3_train_',exp_num,'.bed'), aux3_train_raw) # bed file
aux3info_file <- read.table(paste0(PATH, aux3_info_file_path), sep=" ", header=TRUE)
aux3_bim_data <- as.data.frame(aux3info_file$chr)
aux3_bim_data$snp <- aux3info_file$rsid
aux3_bim_data$zeros <- rep(0, dim(aux3_train_raw)[1])
aux3_bim_data$pos <- aux3info_file$bp
aux3_bim_data$eff <- aux3info_file$effAllele
aux3_bim_data$ref <- aux3info_file$refAllele
write.table(aux3_bim_data, paste0(outPATH,'aux3_train_',exp_num,'.bim'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #bim file
aux3_fam_data <- as.data.frame(1:dim(aux3_train_raw)[2])
aux3_fam_data$col2 <- aux3_fam_data[,1]
aux3_fam_data$col3 <- rep(0, length(aux3_fam_data[,1]))
aux3_fam_data$col4 <- rep(0, length(aux3_fam_data[,1]))
aux3_fam_data$col5 <- rep(0, length(aux3_fam_data[,1]))
aux3_fam_data$col6 <- rep(-9, length(aux3_fam_data[,1]))
write.table(aux3_fam_data, paste0(outPATH,'aux3_train_',exp_num,'.fam'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #fam file

#Aux 4
write_bed(paste0(outPATH,'aux4_train_',exp_num,'.bed'), aux4_train_raw) # bed file
aux4info_file <- read.table(paste0(PATH, aux4_info_file_path), sep=" ", header=TRUE)
aux4_bim_data <- as.data.frame(aux4info_file$chr)
aux4_bim_data$snp <- aux4info_file$rsid
aux4_bim_data$zeros <- rep(0, dim(aux4_train_raw)[1])
aux4_bim_data$pos <- aux4info_file$bp
aux4_bim_data$eff <- aux4info_file$effAllele
aux4_bim_data$ref <- aux4info_file$refAllele
write.table(aux4_bim_data, paste0(outPATH,'aux4_train_',exp_num,'.bim'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #bim file
aux4_fam_data <- as.data.frame(1:dim(aux4_train_raw)[2])
aux4_fam_data$col2 <- aux4_fam_data[,1]
aux4_fam_data$col3 <- rep(0, length(aux4_fam_data[,1]))
aux4_fam_data$col4 <- rep(0, length(aux4_fam_data[,1]))
aux4_fam_data$col5 <- rep(0, length(aux4_fam_data[,1]))
aux4_fam_data$col6 <- rep(-9, length(aux4_fam_data[,1]))
write.table(aux4_fam_data, paste0(outPATH,'aux4_train_',exp_num,'.fam'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #fam file

#Tar Train
write_bed(paste0(outPATH,'tar_train_',exp_num,'.bed'), tar_train_raw) # bed file
tarinfo_file <- read.table(paste0(PATH, tar_info_file_path), sep=" ", header=TRUE)
tar_train_bim_data <- as.data.frame(tarinfo_file$chr)
tar_train_bim_data$snp <- tarinfo_file$rsid
tar_train_bim_data$zeros <- rep(0, dim(tar_train_raw)[1])
tar_train_bim_data$pos <- tarinfo_file$bp
tar_train_bim_data$eff <- tarinfo_file$effAllele
tar_train_bim_data$ref <- tarinfo_file$refAllele
write.table(tar_train_bim_data, paste0(outPATH,'tar_train_',exp_num,'.bim'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #bim file
tar_train_fam_data <- as.data.frame(1:dim(tar_train_raw)[2])
tar_train_fam_data$col2 <- tar_train_fam_data[,1]
tar_train_fam_data$col3 <- rep(0, length(tar_train_fam_data[,1]))
tar_train_fam_data$col4 <- rep(0, length(tar_train_fam_data[,1]))
tar_train_fam_data$col5 <- rep(0, length(tar_train_fam_data[,1]))
tar_train_fam_data$col6 <- rep(-9, length(tar_train_fam_data[,1]))
write.table(tar_train_fam_data, paste0(outPATH,'tar_train_',exp_num,'.fam'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #fam file

#Tar Val
write_bed(paste0(outPATH,'tar_val_',exp_num,'.bed'), tar_val_raw) # bed file
tarinfo_file <- read.table(paste0(PATH, tar_info_file_path), sep=" ", header=TRUE)
tar_val_bim_data <- as.data.frame(tarinfo_file$chr)
tar_val_bim_data$snp <- tarinfo_file$rsid
tar_val_bim_data$zeros <- rep(0, dim(tar_val_raw)[1])
tar_val_bim_data$pos <- tarinfo_file$bp
tar_val_bim_data$eff <- tarinfo_file$effAllele
tar_val_bim_data$ref <- tarinfo_file$refAllele
write.table(tar_val_bim_data, paste0(outPATH,'tar_val_',exp_num,'.bim'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #bim file
tar_val_fam_data <- as.data.frame(1:dim(tar_val_raw)[2])
tar_val_fam_data$col2 <- tar_val_fam_data[,1]
tar_val_fam_data$col3 <- rep(0, length(tar_val_fam_data[,1]))
tar_val_fam_data$col4 <- rep(0, length(tar_val_fam_data[,1]))
tar_val_fam_data$col5 <- rep(0, length(tar_val_fam_data[,1]))
tar_val_fam_data$col6 <- rep(-9, length(tar_val_fam_data[,1]))
write.table(tar_val_fam_data, paste0(outPATH,'tar_val_',exp_num,'.fam'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #fam file

#Tar Test
write_bed(paste0(outPATH,'tar_test_',exp_num,'.bed'), tar_test_raw) # bed file
tarinfo_file <- read.table(paste0(PATH, tar_info_file_path), sep=" ", header=TRUE)
tar_test_bim_data <- as.data.frame(tarinfo_file$chr)
tar_test_bim_data$snp <- tarinfo_file$rsid
tar_test_bim_data$zeros <- rep(0, dim(tar_test_raw)[1])
tar_test_bim_data$pos <- tarinfo_file$bp
tar_test_bim_data$eff <- tarinfo_file$effAllele
tar_test_bim_data$ref <- tarinfo_file$refAllele
write.table(tar_test_bim_data, paste0(outPATH,'tar_test_',exp_num,'.bim'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #bim file
tar_test_fam_data <- as.data.frame(1:dim(tar_test_raw)[2])
tar_test_fam_data$col2 <- tar_test_fam_data[,1]
tar_test_fam_data$col3 <- rep(0, length(tar_test_fam_data[,1]))
tar_test_fam_data$col4 <- rep(0, length(tar_test_fam_data[,1]))
tar_test_fam_data$col5 <- rep(0, length(tar_test_fam_data[,1]))
tar_test_fam_data$col6 <- rep(-9, length(tar_test_fam_data[,1]))
write.table(tar_test_fam_data, paste0(outPATH,'tar_test_',exp_num,'.fam'), sep='\t', 
            row.names = FALSE, col.names = FALSE, quote = FALSE) #fam file

gc()
gc()
