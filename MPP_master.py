# LOAD ALL NECESSARY LIBRARIES =======================================================================
import os
import csv
import sys
import torch
import random
import warnings
import argparse
import subprocess
import scipy as sp
import numpy as np
import pandas as pd
from tqdm import tqdm
import concurrent.futures
import statsmodels.api as sm
from functools import partial
from datetime import datetime
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
from statsmodels.formula.api import ols
from pandas_plink import read_plink1_bin, read_plink
from sklearn.metrics import precision_recall_curve, auc
# =====================================================================================================

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected (True/False).')

def check_chr_files(prefix, chrs, extensions):
    """
    Helper function to verify that all per-chromosome files exist in a directory.
    Takes a prefix (e.g., 'dir/tar_train_chr'), appends the chromosome number and extension.
    """
    missing = []
    for c in chrs:
        for ext in extensions:
            file_path = f"{prefix}{c}{ext}"
            if not os.path.isfile(file_path):
                missing.append(file_path)
    return missing

def setup_parser():
    parser = argparse.ArgumentParser(
        description="Master script for MultiPopPred (MPP) pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Version
    parser.add_argument('--version', type=str, default='MPP-PRS+', 
                        choices=['MPP-PRS+', 'MPP-PRS', 'MPP-PRS-TarSS', 'MPP-GWAS', 'MPP-GWAS-TarSS', 'MPP-GWAS-Admix'],
                        help="Specify which version of MultiPopPred to use.")

    # Chromosomes
    parser.add_argument('--chr', type=str, default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22",
                        help="Comma-separated integers (1-22) specifying chromosomes to process.")

    # Target Population Training & Files 
    parser.add_argument('--bfileTrain', type=str, help="Prefix path to target population's training genotype data (e.g., /path/tar_train_chr).")
    parser.add_argument('--bfileExtLD', type=str, help="Prefix path to target population's external LD reference file (e.g., /path/ext_ld_chr).")
    parser.add_argument('--phenoTrain', type=str, help="Path to the target population's ground-truth training phenotype.")
    parser.add_argument('--ssTar', type=str, help="Path to the target population's summary statistics file (.txt).")
    parser.add_argument('--prsTar', type=str, help="Path to the target population's single ancestry PRS file (.txt).")
    parser.add_argument('--covTrain', type=str, help="Prefix path to target population's training set covariates (Optional).")

    # Auxiliary Populations
    parser.add_argument('--auxPops', type=str, help="Comma-separated auxiliary populations (e.g., EUR,EAS,AMR,AFR).")
    parser.add_argument('--ssAux', type=str, help="Comma-separated paths to auxiliary populations' summary statistics.")
    parser.add_argument('--prsAux', type=str, help="Comma-separated paths to auxiliary populations' single ancestry PRS (.txt).")
    parser.add_argument('--AdmixWeights', type=str, help="Path to file containing admixture proportions (for MPP-GWAS-Admix).")

    parser.add_argument('--trait_type', type=str, choices=['binary','continuous'], help="Specify whether the trait is continuous or binary.")
    parser.add_argument('--Neff', type=str, help="Comma-separated effective sample size for auxiliary as well as target population. Target population Neff should be at the end.")

    # Validation
    parser.add_argument('--validate', type=str2bool, default=True, help="True/False to enable validation. Default is True.")
    parser.add_argument('--bfileVal', type=str, help="Prefix path to validation genotype data (e.g., /path/tar_val_chr).")
    parser.add_argument('--phenoVal', type=str, help="Path to the validation ground-truth phenotype.")
    parser.add_argument('--covVal', type=str, help="Prefix path to validation covariates (Optional).")

    # Penalties
    parser.add_argument('--L1', type=str, default="0.0001,0.001,0.01,0.1,1,10,100,1000,10000", help="Comma-separated L1 penalty floats.")
    parser.add_argument('--L2', type=str, default="0.0001,0.001,0.01,0.1,0.5,0.9,0.99", help="Comma-separated L2 penalty floats.")

    # LBFGS & System Params
    parser.add_argument('--smoothing', type=float, default=0.1, help="Smoothing parameter.")
    parser.add_argument('--max_iter', type=int, default=10000, help="Maximum number of iterations for LBFGS.")
    parser.add_argument('--tol', type=float, default=0.0001, help="Tolerance for LBFGS.")
    parser.add_argument('--max_fun', type=int, default=10, help="Maximum number of function evaluations for LBFGS.")
    parser.add_argument('--seed', type=int, default=123, help="Seed for reproducibility.")
    parser.add_argument('--Nthreads', type=str, default='all', help="Number of threads (integer or 'all').")
    parser.add_argument('--out', type=str, required=True, help="Path to the output directory.")
    
    # Verbose Output
    parser.add_argument('--verbose', type=str2bool, default=True, help="True/False to display the verbose outputs.")

    return parser

def crosscheck_args(args):
    """Handles all custom dependency and constraint cross checking."""
    
    # --- Chromosome Validation ---
    try:
        chrs = [int(c.strip()) for c in args.chr.split(',')]
        if any(c < 1 or c > 22 for c in chrs):
            raise ValueError
    except ValueError:
        sys.exit("Error: --chr must contain only comma-separated integers between 1 and 22.")
    
    args.chr_list = chrs 

    # --- Version Group Definitions ---
    v_train_group = ['MPP-PRS+', 'MPP-PRS', 'MPP-GWAS', 'MPP-GWAS-Admix']
    v_extld_group = ['MPP-PRS-TarSS', 'MPP-GWAS-TarSS']
    v_ssaux_group = ['MPP-GWAS', 'MPP-GWAS-TarSS', 'MPP-GWAS-Admix']
    v_prsaux_group = ['MPP-PRS+', 'MPP-PRS', 'MPP-PRS-TarSS']
    v_l2_group = ['MPP-GWAS-TarSS', 'MPP-PRS-TarSS']

    plink_exts = ['.bed', '.bim', '.fam']
    cov_exts = ['.txt'] 

    # --- Argument Validation Logic ---
    
    if args.version not in ['MPP-PRS+','MPP-PRS','MPP-PRS-TarSS','MPP-GWAS','MPP-GWAS-TarSS','MPP-GWAS-Admix']:
        sys.exit(f"Error: Invalid --version {args.version} entered.")
        
    if args.version in v_extld_group:
      if args.trait_type == 'binary':
        sys.exit(f"Error: Version {args.version} is not supported for binary phenotypes.")
    
    if args.version in v_train_group:
        if not args.bfileTrain or not args.phenoTrain:
            sys.exit(f"Error: --bfileTrain and --phenoTrain are required for version {args.version}.")
        
        missing_train = check_chr_files(args.bfileTrain, chrs, plink_exts)
        if missing_train:
            sys.exit(f"Error: Missing --bfileTrain PLINK files for specified chromosomes. Could not find:\n" + "\n".join(missing_train[:3]) + "\n...etc")

        if args.covTrain:
            missing_cov = check_chr_files(args.covTrain, chrs, cov_exts)
            if missing_cov:
                sys.exit(f"Error: Missing --covTrain files for specified chromosomes. Could not find:\n" + "\n".join(missing_cov[:3]) + "\n...etc")
    else:
        if args.bfileTrain or args.phenoTrain or args.covTrain:
            sys.exit(f"Error: --bfileTrain, --phenoTrain, and --covTrain cannot be used with version {args.version}.")

    if args.version in v_extld_group:
        if not args.bfileExtLD:
            sys.exit(f"Error: --bfileExtLD is required for version {args.version}.")
        missing_extld = check_chr_files(args.bfileExtLD, chrs, plink_exts)
        if missing_extld:
            sys.exit(f"Error: Missing --bfileExtLD PLINK files for specified chromosomes. Could not find:\n" + "\n".join(missing_extld[:3]) + "\n...etc")
    elif args.bfileExtLD:
        sys.exit(f"Error: --bfileExtLD cannot be used with version {args.version}.")

    if args.version == 'MPP-GWAS-TarSS':
        if not args.ssTar: sys.exit("Error: --ssTar is required for version MPP-GWAS-TarSS.")
    elif args.ssTar: sys.exit(f"Error: --ssTar cannot be used with version {args.version}.")

    if args.version == 'MPP-PRS-TarSS':
        if not args.prsTar: sys.exit("Error: --prsTar is required for version MPP-PRS-TarSS.")
    elif args.prsTar: sys.exit(f"Error: --prsTar cannot be used with version {args.version}.")

    num_aux = len(args.auxPops.split(',')) if args.auxPops else 0

    if args.version in v_ssaux_group:
        if not args.ssAux: sys.exit(f"Error: --ssAux is required for version {args.version}.")
        if len(args.ssAux.split(',')) != num_aux: sys.exit("Error: --ssAux paths must match --auxPops count.")
    elif args.ssAux: sys.exit(f"Error: --ssAux cannot be used with version {args.version}.")

    if args.version in v_prsaux_group:
        if not args.prsAux: sys.exit(f"Error: --prsAux is required for version {args.version}.")
        if len(args.prsAux.split(',')) != num_aux: sys.exit("Error: --prsAux paths must match --auxPops count.")
    elif args.prsAux: sys.exit(f"Error: --prsAux cannot be used with version {args.version}.")

    if args.version == 'MPP-GWAS-Admix':
        if not args.AdmixWeights: sys.exit("Error: --AdmixWeights is required for version MPP-GWAS-Admix.")
    elif args.AdmixWeights: sys.exit(f"Error: --AdmixWeights cannot be used with version {args.version}.")
    
    if args.validate:
        if not args.bfileVal or not args.phenoVal:
            sys.exit("Error: --bfileVal and --phenoVal are required when --validate is True.")
        
        missing_val = check_chr_files(args.bfileVal, chrs, plink_exts)
        if missing_val: sys.exit(f"Error: Missing --bfileVal PLINK files.")

        if args.covVal:
            missing_cov_val = check_chr_files(args.covVal, chrs, cov_exts)
            if missing_cov_val: sys.exit(f"Error: Missing --covVal files.")
    else:
        if args.bfileVal or args.phenoVal or args.covVal:
             sys.exit("Error: --bfileVal, --phenoVal, and --covVal cannot be used when --validate is False.")
        
        if len(args.L1.split(',')) > 1: sys.exit("Error: --L1 must be a single float when --validate is False.")
        if args.version in v_l2_group and len(args.L2.split(',')) > 1:
            sys.exit("Error: --L2 must be a single float when --validate is False.")

    return args

def print_summary(args):
    """Prints a clean, version-specific summary of all parsed inputs."""
    
    v_train_group = ['MPP-PRS+', 'MPP-PRS', 'MPP-GWAS', 'MPP-GWAS-Admix']
    v_extld_group = ['MPP-PRS-TarSS', 'MPP-GWAS-TarSS']
    v_ssaux_group = ['MPP-GWAS', 'MPP-GWAS-TarSS', 'MPP-GWAS-Admix']
    v_prsaux_group = ['MPP-PRS+', 'MPP-PRS', 'MPP-PRS-TarSS']
    v_l2_group = ['MPP-GWAS-TarSS', 'MPP-PRS-TarSS']

    print("\n" + "="*60)
    print(f" MultiPopPred (MPP) Run Summary")
    print("="*60)
    
    print(f"Version:           {args.version}")
    print(f"Chromosomes:       {args.chr}")
    print(f"Output Directory:  {args.out}")
    print(f"Seed:              {args.seed}")
    print(f"Threads:           {args.Nthreads}")
    print("-" * 60)
    
    print("TARGET POPULATION INPUTS:")
    if args.version in v_train_group:
        print(f"  Train Genotypes: {args.bfileTrain}")
        print(f"  Train Phenotype: {args.phenoTrain}")
        print(f"  Train Covariate: {args.covTrain if args.covTrain else 'Not specified (No adjustment)'}")
    
    if args.version in v_extld_group:
        print(f"  External LD Ref: {args.bfileExtLD}")
        
    if args.version == 'MPP-GWAS-TarSS':
        print(f"  Target SumStats: {args.ssTar}")
        
    if args.version == 'MPP-PRS-TarSS':
        print(f"  Target PRS:      {args.prsTar}")

    print("-" * 60)
    
    if args.version in v_ssaux_group or args.version in v_prsaux_group:
        print("AUXILIARY POPULATION INPUTS:")
        print(f"  Populations:     {args.auxPops if args.auxPops else 'Not specified'}")
        if args.version in v_ssaux_group:
            print(f"  Aux SumStats:    {args.ssAux}")
        if args.version in v_prsaux_group:
            print(f"  Aux PRS:         {args.prsAux}")
        if args.version == 'MPP-GWAS-Admix':
            print(f"  Admix Weights:   {args.AdmixWeights}")
        print("-" * 60)
        
    print("VALIDATION SETTINGS:")
    print(f"  Perform Validate: {args.validate}")
    if args.validate:
        print(f"  Val Genotypes:    {args.bfileVal}")
        print(f"  Val Phenotype:    {args.phenoVal}")
        print(f"  Val Covariate:    {args.covVal if args.covVal else 'Not specified (No adjustment)'}")
    print("-" * 60)
    
    print("MODEL PARAMETERS:")
    print(f"  L1 Penalty:       {args.L1}")
    if args.version in v_l2_group:
        print(f"  L2 Penalty:       {args.L2}")
    print(f"  Smoothing:        {args.smoothing}")
    print(f"  LBFGS Max Iter:   {args.max_iter}")
    print(f"  LBFGS Tolerance:  {args.tol}")
    print(f"  LBFGS Max Fun:    {args.max_fun}")
    
    print("="*60 + "\n")

    
def pre_process(chr_num, args):
  
  if args.verbose:
    print("="*60)
    print(f"Processing chromosome {chr_num}...")
    print("="*60 + "\n")
  
  if args.version in ['MPP-PRS+','MPP-PRS','MPP-GWAS','MPP-GWAS-Admix']:
      
      """Loads target training genotype data - true LD """
      bfileTrain_path = args.bfileTrain+str(chr_num)
      G_tar_train  = read_plink1_bin(bfileTrain_path + '.bed',\
                               bfileTrain_path + '.bim',\
                               bfileTrain_path + '.fam',\
                               verbose=False)
      X_tar_train = G_tar_train.values
      X_tar_train = np.where(X_tar_train == 2, 0, np.where(X_tar_train == 0, 2, X_tar_train))
      X_tar_train = np.array(X_tar_train, dtype=float)
      all_snps_tar_train = G_tar_train.snp.values
      
      if args.verbose:
        print("Loaded "+str(X_tar_train.shape[0])+" samples and "+str(X_tar_train.shape[1])+" SNPs from target training genotype file...")
      
      """Loads target training phenotype data - true LD """
      Y_train = pd.read_csv(args.phenoTrain, delimiter='\t')
      if Y_train.shape[0] != X_tar_train.shape[0]:
        sys.exit(f"Error: Target training data - X and Y dimensions do not match!")
      
      if args.covTrain:
        """Loads target training covariates """
        Cov_train = pd.read_csv(args.covTrain, delimiter='\t')
        if Cov_train.shape[0] != Y_train.shape[0]:
          sys.exit(f"Target training data - Cov and Y dimensions do not match!")
        if args.verbose:
          print("Loaded "+str(Cov_train.shape[0])+" samples and "+str(Cov_train.shape[1]-2)+" covariates from target training covariates file...")
  
  # ===========================================================================
  # ===========================================================================
      
  if args.version in ['MPP-PRS-TarSS','MPP-GWAS-TarSS']:
      
      """Loads target training genotype data - external LD """
      bfileExtLD_path = args.bfileExtLD+str(chr_num)
      G_tar_train  = read_plink1_bin(bfileExtLD_path + '.bed',\
                               bfileExtLD_path + '.bim',\
                               bfileExtLD_path + '.fam',\
                               verbose=False)
      X_tar_train = G_tar_train.values
      X_tar_train = np.where(X_tar_train == 2, 0, np.where(X_tar_train == 0, 2, X_tar_train))
      X_tar_train = np.array(X_tar_train, dtype=float)
      all_snps_tar_train = G_tar_train.snp.values
      
      if args.verbose:
        print("Loaded "+str(X_tar_train.shape[0])+" samples and "+str(X_tar_train.shape[1])+" SNPs from target training external genotype file...")
      
      if args.version == 'MPP-GWAS-TarSS':
        """Loads target summary statistics """
        tar_ss = pd.read_csv(args.ssTar, delimiter='\t')
      elif args.version == 'MPP-PRS-TarSS':
        """Loads target PRS """
        tar_prs = pd.read_csv(args.prsTar, delimiter='\t')
    
  # ===========================================================================
  # ===========================================================================  
    
  if args.version in ['MPP-PRS+','MPP-PRS','MPP-PRS-TarSS']:
      
      """Loads auxiliary PRS data into a dictionary."""
      # 1. Split the comma-separated strings into lists
      pops = args.auxPops.split(',')
      prs_paths = args.prsAux.split(',')

      # 2. Initialize an empty dictionary to hold your K datasets
      aux_prs_data = {}

      # 3. Loop through the populations and paths simultaneously using zip()
      for pop, path in zip(pops, prs_paths):
          if args.verbose:
              print(f"Loading auxiliary PRS data for {pop} from {path}...")

          # Read the file (assuming tab-separated, adjust 'sep' if it's a CSV or space-delimited)
          aux_prs_data[pop] = pd.read_csv(path+str(chr_num)+'.txt', delimiter='\t')
          aux_prs_data[pop] = aux_prs_data[pop].dropna()
          if args.verbose:
              print(str(aux_prs_data[pop].shape[0])+' SNPs loaded from '+str(pop)+' PRS ...')
              
      if args.version in ['MPP-PRS+','MPP-PRS']:
        snp_sets = [set(df['SNP']) for df in aux_prs_data.values()]
        snp_sets.append(set(all_snps_tar_train))
        common_snps = set.intersection(*snp_sets)
      elif args.version == 'MPP-PRS-TarSS':
        snp_sets = [set(df['SNP']) for df in aux_prs_data.values()]
        snp_sets.append(set(all_snps_tar_train))
        snp_sets.append(set(tar_prs['SNP']))
        common_snps = set.intersection(*snp_sets)
      
      if args.verbose:
        print(str(len(common_snps))+' SNPs common across auxiliary and target populations ...')
      
      for pop in pops:
        aux_prs_data[pop] = aux_prs_data[pop][aux_prs_data[pop]['SNP'].isin(common_snps)].reset_index(drop=True)
      
      """Aggregate auxiliary BETAs using equal weighting"""  
      aux_betas_list = [np.array(df['BETA'], dtype=float) for df in aux_prs_data.values()]
      aux_betas = np.mean(np.stack(aux_betas_list), axis=0)
      
      snp_order = [list(df['SNP']) for df in aux_prs_data.values()][0]
      
      if args.version in ['MPP-PRS+','MPP-PRS']:
        X_tar_bim = pd.read_csv(args.bfileTrain+str(chr_num)+'.bim', sep='\t', header=None, names=['CHR','SNP','zeros','POS','A1','A2'])
        X_tar_bim = X_tar_bim[X_tar_bim['SNP'].isin(common_snps)].reset_index(drop=True)
        tar_snps = list(X_tar_bim['SNP'])
        tar_chr = list(X_tar_bim['CHR'])
        tar_pos = list(X_tar_bim['POS'])
        tar_a1 = list(X_tar_bim['A1'])
        tar_a2 = list(X_tar_bim['A2'])
      elif args.version == 'MPP-PRS-TarSS':
        tar_prs = tar_prs[tar_prs['SNP'].isin(common_snps)].reset_index(drop=True)
        tar_snps = list(tar_prs['SNP'])
        tar_chr = list(tar_prs['CHR'])
        tar_pos = list(tar_prs['POS'])
        tar_a1 = list(tar_prs['A1'])
        tar_a2 = list(tar_prs['A2'])
      
      """Filter genotype data to retain only common SNPs"""
      temp_indices_tar = [i for i, e in enumerate(all_snps_tar_train) if e in common_snps]
      all_snps_tar_train = [i for i in all_snps_tar_train if i in common_snps]
      X_tar_train = X_tar_train[:, temp_indices_tar]
      if not (all_snps_tar_train == snp_order):
        sys.exit(f"Error: SNP order in the target genotype matrix and aggregated auxiliary betas does not match!")
      
  if args.version in ['MPP-GWAS','MPP-GWAS-TarSS','MPP-GWAS-Admix']:
    
    """Loads auxiliary Summary Statistics data into a dictionary."""

    # 1. Split the comma-separated strings into lists
    pops = args.auxPops.split(',')
    ss_paths = args.ssAux.split(',')

    # 2. Initialize an empty dictionary to hold your K datasets
    aux_ss_data = {}

    # 3. Loop through the populations and paths simultaneously using zip()
    for pop, path in zip(pops, ss_paths):
        if args.verbose:
            print(f"Loading auxiliary Summ Stats data for {pop} from {path}...")

        # Read the file (assuming tab-separated, adjust 'sep' if it's a CSV or space-delimited)
        aux_ss_data[pop] = pd.read_csv(path+str(chr_num)+'.txt', delimiter='\t')
        aux_ss_data[pop] = aux_ss_data[pop].dropna()
        if args.verbose:
            print(str(aux_ss_data[pop].shape[0])+' SNPs loaded from '+str(pop)+' Summ Stats ...')
            
    
    if args.version in ['MPP-GWAS','MPP-GWAS-Admix']:
      snp_sets = [set(df['SNP']) for df in aux_ss_data.values()]
      snp_sets.append(set(all_snps_tar_train))
      common_snps = set.intersection(*snp_sets)
    elif args.version == 'MPP-GWAS-TarSS':
      snp_sets = [set(df['SNP']) for df in aux_ss_data.values()]
      snp_sets.append(set(all_snps_tar_train))
      snp_sets.append(set(tar_ss['SNP']))
      common_snps = set.intersection(*snp_sets)
    
    if args.verbose:
      print(str(len(common_snps))+' SNPs common across auxiliary and target populations ...')
      
    for pop in pops:
        aux_ss_data[pop] = aux_ss_data[pop][aux_ss_data[pop]['SNP'].isin(common_snps)].reset_index(drop=True)
      
    if args.version in ['MPP-GWAS','MPP-GWAS-TarSS']:
      """Aggregate auxiliary BETAs using equal weighting"""
      aux_betas_list = [np.array(df['BETA'], dtype=float) for df in aux_ss_data.values()]
      aux_betas = np.mean(np.stack(aux_betas_list), axis=0)
      
    elif args.version == 'MPP-GWAS-Admix':
      """Aggregate auxiliary BETAs using admix weighting"""
      aux_betas_list = [np.array(df['BETA'], dtype=float) for df in aux_ss_data.values()]
      aux_betas_matrix = np.stack(aux_betas_list)
      admix_wts = pd.read_csv(args.AdmixWeights, delimiter=' ', header=None)
      admix_wts = admix_wts.to_numpy()
      temp_wts = admix_wts @ aux_betas_matrix
      aux_betas = np.median(temp_wts, axis=0)
      
    snp_order = [list(df['SNP']) for df in aux_ss_data.values()][0]
    
    if args.version in ['MPP-GWAS','MPP-GWAS-Admix']:
      X_tar_bim = pd.read_csv(args.bfileTrain+str(chr_num)+'.bim', sep='\t', header=None, names=['CHR','SNP','zeros','POS','A1','A2'])
      X_tar_bim = X_tar_bim[X_tar_bim['SNP'].isin(common_snps)].reset_index(drop=True)
      tar_snps = list(X_tar_bim['SNP'])
      tar_chr = list(X_tar_bim['CHR'])
      tar_pos = list(X_tar_bim['POS'])
      tar_a1 = list(X_tar_bim['A1'])
      tar_a2 = list(X_tar_bim['A2'])
    elif args.version == 'MPP-GWAS-TarSS':
      tar_ss = tar_ss[tar_ss['SNP'].isin(common_snps)].reset_index(drop=True)
      tar_snps = list(tar_ss['SNP'])
      tar_chr = list(tar_ss['CHR'])
      tar_pos = list(tar_ss['POS'])
      tar_a1 = list(tar_ss['A1'])
      tar_a2 = list(tar_ss['A2'])
    
    """Filter genotype data to retain only common SNPs"""
    temp_indices_tar = [i for i, e in enumerate(all_snps_tar_train) if e in common_snps]
    all_snps_tar_train = [i for i in all_snps_tar_train if i in common_snps]
    X_tar_train = X_tar_train[:, temp_indices_tar]
    if not (all_snps_tar_train == snp_order):
      sys.exit(f"Error: SNP order in the target genotype matrix and aggregated auxiliary betas does not match!")
    
  # ===========================================================================
  # ===========================================================================  
    
  """Create a dictionary to return all processed data together"""
  if args.version in ['MPP-PRS+','MPP-PRS','MPP-GWAS','MPP-GWAS-Admix']:
    my_dict = {
      "X_tar_train": X_tar_train,
      "Y_tar_train": np.array(Y_train['Pheno'], dtype=float),
      "aux_betas": aux_betas,
      "tar_snps": tar_snps,
      "tar_chr": tar_chr,
      "tar_pos": tar_pos,
      "tar_a1": tar_a1,
      "tar_a2": tar_a2
    }
    if args.covTrain:
      my_dict.update({
        "Cov_train": Cov_train
      })
        
  elif args.version == 'MPP-PRS-TarSS':
    my_dict = {
      "X_tar_train": X_tar_train,
      "aux_betas": aux_betas,
      "tar_prs": tar_prs,
      "tar_snps": tar_snps,
      "tar_chr": tar_chr,
      "tar_pos": tar_pos,
      "tar_a1": tar_a1,
      "tar_a2": tar_a2
    }
      
  elif args.version == 'MPP-GWAS-TarSS':
    my_dict = {
      "X_tar_train": X_tar_train,
      "aux_betas": aux_betas,
      "tar_ss": tar_ss,
      "tar_snps": tar_snps,
      "tar_chr": tar_chr,
      "tar_pos": tar_pos,
      "tar_a1": tar_a1,
      "tar_a2": tar_a2
    }

  return my_dict

def doMPP(chr_num, L1_penalty, L2_penalty, args):
  
  processed_data = pre_process(chr_num, args)
  
  # ==============================================================================================
  # ==============================================================================================
  
  if args.trait_type == 'continuous':
  
    if args.version in ['MPP-PRS+','MPP-PRS','MPP-GWAS','MPP-GWAS-Admix']:
      X_tar_train = processed_data['X_tar_train']
      Y_tar_train = processed_data['Y_tar_train']
      aux_betas = processed_data['aux_betas']
      tar_snps = processed_data['tar_snps']
      tar_chr = processed_data['tar_chr']
      tar_pos = processed_data['tar_pos']
      tar_a1 = processed_data['tar_a1']
      tar_a2 = processed_data['tar_a2']
      
      if args.covTrain:
        Cov_train = processed_data['Cov_train']
        temp_df = pd.concat([Y_tar_train[['Pheno']],Cov_train.iloc[:,2:]],axis=1)
        model = ols('Pheno ~ .', data=temp_df).fit()
        temp_y_hat = model.predict(temp_df)
        Y_tar_train = np.array(Y_tar_train['Pheno'], dtype=float) - temp_y_hat
      
      Beta_hat_init = np.random.normal(0, 0.00000000001, len(tar_snps))
      
      X_tensor = torch.from_numpy(X_tar_train)
      XtX = torch.matmul(X_tensor.T, X_tensor)
      XtX_np = XtX.numpy()
  
      y_tensor = torch.from_numpy(Y_tar_train)
      Xty = torch.matmul(X_tensor.T, y_tensor)
      Xty_np = Xty.numpy()
      
      mu = float(args.smoothing)
      L1_penalty = L1_penalty
      L2_penalty = L2_penalty # Not used in these versions
      max_iter = int(args.max_iter)
      tol = float(args.tol)
      max_fun = int(args.max_fun)
      
      def costfunc_lin(Bs, L1_penalty):
        t1 = torch.matmul(X_tensor, torch.from_numpy(Bs)).numpy()
        temp = Bs - aux_betas
        nesterov = mu*np.log((0.5*np.exp(-temp/mu))+(0.5*np.exp(temp/mu)))
        term2 = L1_penalty*sum(nesterov)
        val = sum((Y_tar_train - t1)**2) + term2
        return val
      
      def gradfunc_lin(Bs, L1_penalty):
          term1 = 2*torch.matmul(XtX, torch.from_numpy(Bs)).numpy()
          term2 = 2*Xty_np
          temp = Bs - aux_betas
          nesterov = np.divide((-np.exp(-temp/mu) + np.exp(temp/mu)),(np.exp(-temp/mu) + np.exp(temp/mu)))
          term3 = L1_penalty*nesterov
          return term1-term2+term3
  
      ans = sp.optimize.minimize(costfunc_lin, jac=gradfunc_lin, x0=Beta_hat_init, args=(L1_penalty),\
                                 method='L-BFGS-B', options={'maxiter':max_iter, 'ftol':tol,'maxfun':max_fun})
                                 
      final_snps = tar_snps
      final_betas = ans.x
      final_chr = tar_chr
      final_pos = tar_pos
      final_a1 = tar_a1
      final_results_df = pd.DataFrame({'CHR':final_chr, 'SNP':final_snps, 'POS':final_pos,\
                                       'A1':final_a1, 'BETA':final_betas,})
      
    if args.version == 'MPP-PRS-TarSS':
      X_tar_train = processed_data['X_tar_train']
      aux_betas = processed_data['aux_betas']
      tar_prs = processed_data['tar_prs']
      tar_snps = processed_data['tar_snps']
      tar_chr = processed_data['tar_chr']
      tar_pos = processed_data['tar_pos']
      tar_a1 = processed_data['tar_a1']
      tar_a2 = processed_data['tar_a2']
      
      Neff = int(args.Neff.split(',')[-1])
      tar_corr = np.array(tar_prs['BETA'], dtype=float)*Neff
      
      Beta_hat_init = np.random.normal(0, 0.00000000001, len(tar_snps))
      
      mu = float(args.smoothing)
      L1_penalty = L1_penalty
      L2_penalty = L2_penalty
      max_iter = int(args.max_iter)
      tol = float(args.tol)
      max_fun = int(args.max_fun)
      
      X_tensor = torch.from_numpy(X_tar_train)
      XtX = torch.matmul(X_tensor.T, X_tensor)
      XtX_np = XtX.numpy()
      
      XtX_np = abs(1-L2_penalty)*XtX_np
      r = tar_corr
      
      def costfunc_lin_tarss(Bs, L1_penalty, L2_penalty):
          t1 = torch.matmul(torch.from_numpy(Bs).T, torch.from_numpy(XtX_np))
          term1 = torch.matmul(t1, torch.from_numpy(Bs)).numpy()
          term2 = 2*torch.matmul(torch.from_numpy(Bs).T, torch.from_numpy(r)).numpy()
          term3 = L2_penalty*np.matmul(Bs.T,Bs)
          temp = Bs - aux_betas
          nesterov = mu*np.log((0.5*np.exp(-temp/mu))+(0.5*np.exp(temp/mu)))
          term4 = L1_penalty*sum(nesterov)
          return term1-term2+term3+term4
  
      def gradfunc_lin_tarss(Bs, L1_penalty, L2_penalty):
          term1 = 2*torch.matmul(torch.from_numpy(XtX_np),torch.from_numpy(Bs)).numpy()
          term2 = 2*r
          term3 = L2_penalty*2*Bs
          temp = Bs - aux_betas
          nesterov = np.divide((-np.exp(-temp/mu) + np.exp(temp/mu)),(np.exp(-temp/mu) + np.exp(temp/mu)))
          term4 = L1_penalty*nesterov
          return term1-term2+term3+term4
        
      ans = sp.optimize.minimize(costfunc_lin_tarss, jac=gradfunc_lin_tarss, x0=Beta_hat_init,\
                                 args=(L1_penalty, L2_penalty), method='L-BFGS-B',\
                                 options={'maxiter':max_iter, 'ftol':tol,'maxfun':max_fun})
                                 
      final_snps = tar_snps
      final_betas = ans.x
      final_chr = tar_chr
      final_pos = tar_pos
      final_a1 = tar_a1
      final_results_df = pd.DataFrame({'CHR':final_chr, 'SNP':final_snps, 'POS':final_pos,\
                                       'A1':final_a1, 'BETA':final_betas,})
                                       
    if args.version == 'MPP-GWAS-TarSS':
      X_tar_train = processed_data['X_tar_train']
      aux_betas = processed_data['aux_betas']
      tar_ss = processed_data['tar_ss']
      tar_snps = processed_data['tar_snps']
      tar_chr = processed_data['tar_chr']
      tar_pos = processed_data['tar_pos']
      tar_a1 = processed_data['tar_a1']
      tar_a2 = processed_data['tar_a2']
      
      Neff = int(args.Neff.split(',')[-1])
      tar_corr = np.array(tar_ss['BETA'], dtype=float)*Neff
      
      Beta_hat_init = np.random.normal(0, 0.00000000001, len(tar_snps))
      
      mu = float(args.smoothing)
      L1_penalty = L1_penalty
      L2_penalty = L2_penalty
      max_iter = int(args.max_iter)
      tol = float(args.tol)
      max_fun = int(args.max_fun)
      
      X_tensor = torch.from_numpy(X_tar_train)
      XtX = torch.matmul(X_tensor.T, X_tensor)
      XtX_np = XtX.numpy()
      
      XtX_np = abs(1-L2_penalty)*XtX_np
      r = tar_corr
      
      def costfunc_lin_tarss(Bs, L1_penalty, L2_penalty):
          t1 = torch.matmul(torch.from_numpy(Bs).T, torch.from_numpy(XtX_np))
          term1 = torch.matmul(t1, torch.from_numpy(Bs)).numpy()
          term2 = 2*torch.matmul(torch.from_numpy(Bs).T, torch.from_numpy(r)).numpy()
          term3 = L2_penalty*np.matmul(Bs.T,Bs)
          temp = Bs - aux_betas
          nesterov = mu*np.log((0.5*np.exp(-temp/mu))+(0.5*np.exp(temp/mu)))
          term4 = L1_penalty*sum(nesterov)
          return term1-term2+term3+term4
  
      def gradfunc_lin_tarss(Bs, L1_penalty, L2_penalty):
          term1 = 2*torch.matmul(torch.from_numpy(XtX_np),torch.from_numpy(Bs)).numpy()
          term2 = 2*r
          term3 = L2_penalty*2*Bs
          temp = Bs - aux_betas
          nesterov = np.divide((-np.exp(-temp/mu) + np.exp(temp/mu)),(np.exp(-temp/mu) + np.exp(temp/mu)))
          term4 = L1_penalty*nesterov
          return term1-term2+term3+term4
        
      ans = sp.optimize.minimize(costfunc_lin_tarss, jac=gradfunc_lin_tarss, x0=Beta_hat_init,\
                                 args=(L1_penalty, L2_penalty), method='L-BFGS-B',\
                                 options={'maxiter':max_iter, 'ftol':tol,'maxfun':max_fun})
                                 
      final_snps = tar_snps
      final_betas = ans.x
      final_chr = tar_chr
      final_pos = tar_pos
      final_a1 = tar_a1
      final_results_df = pd.DataFrame({'CHR':final_chr, 'SNP':final_snps, 'POS':final_pos,\
                                       'A1':final_a1, 'BETA':final_betas,})
      
  # ==============================================================================================
  # ==============================================================================================  
  
  if args.trait_type == 'binary':
    
    if args.version in ['MPP-PRS+','MPP-PRS','MPP-GWAS','MPP-GWAS-Admix']:
      X_tar_train = processed_data['X_tar_train']
      Y_tar_train = processed_data['Y_tar_train']
      aux_betas = processed_data['aux_betas']
      tar_snps = processed_data['tar_snps']
      tar_chr = processed_data['tar_chr']
      tar_pos = processed_data['tar_pos']
      tar_a1 = processed_data['tar_a1']
      tar_a2 = processed_data['tar_a2']
      
      N_cases = sum(Y_tar_train == 1)
      N_controls = sum(Y_tar_train == 0)
      Neff = 4/((1/N_cases) + (1/N_controls))
      
      M = X_tar_train.shape[1]
      
      if args.covTrain:
        Cov_train = processed_data['Cov_train']
        Cov_train = Cov_train.iloc[:,2:Cov_train.shape[1]]
        Cov_train = Cov_train.values
        numCov = Cov_train.shape[1]
        X_tar_train = np.concatenate((Cov_train, X_tar_train), axis = 1)
        
      #Initialize Beta_hats
      Beta_hat_init = np.random.normal(0, 0.00000000001, X_tar_train.shape[1])
      
      X_tensor = torch.from_numpy(X_tar_train)

      def sigmoid(x):
        return 1 / (1 + np.exp(-x))
      
      mu = float(args.smoothing)
      L1_penalty = L1_penalty
      L2_penalty = L2_penalty # Not used in these versions
      max_iter = int(args.max_iter)
      tol = float(args.tol)
      max_fun = int(args.max_fun)
      
      if args.covTrain:
        def costfunc_log(Bs, L1_penalty):
            t1 = torch.matmul(X_tensor, torch.from_numpy(Bs)).numpy()
            p = sigmoid(t1)
            term1 = np.dot(Y_tar_train, np.log(p)) + np.dot((1-Y_tar_train),np.log(1-p))
            temp = Bs[numCov:] - aux_betas
            nesterov = mu*np.log((0.5*np.exp(-temp/mu))+(0.5*np.exp(temp/mu)))
            term2 = L1_penalty*sum(nesterov)
            val = (-1/Neff)*term1 + term2
            return val
    
        def gradfunc_log(Bs, L1_penalty):
            t1 = torch.matmul(X_tensor, torch.from_numpy(Bs)).numpy()
            p = sigmoid(t1)
            t2 = (Y_tar_train - p).T
            term1 = torch.matmul(torch.from_numpy(t2), X_tensor).numpy()
            temp = Bs[numCov:] - aux_betas
            nesterov = np.divide((-np.exp(-temp/mu) + np.exp(temp/mu)),(np.exp(-temp/mu) + np.exp(temp/mu)))
            term2 = L1_penalty*nesterov
            term3 = np.concatenate((np.zeros((numCov,)),term2))
            val = (-1/Neff)*term1 + term3
            return val
      else:
        def costfunc_log(Bs, L1_penalty):
            t1 = torch.matmul(X_tensor, torch.from_numpy(Bs)).numpy()
            p = sigmoid(t1)
            term1 = np.dot(Y_tar_train, np.log(p)) + np.dot((1-Y_tar_train),np.log(1-p))
            temp = Bs - aux_betas
            nesterov = mu*np.log((0.5*np.exp(-temp/mu))+(0.5*np.exp(temp/mu)))
            term2 = L1_penalty*sum(nesterov)
            val = (-1/Neff)*term1 + term2
            return val
    
        def gradfunc_log(Bs, L1_penalty):
            t1 = torch.matmul(X_tensor, torch.from_numpy(Bs)).numpy()
            p = sigmoid(t1)
            t2 = (Y_tar_train - p).T
            term1 = torch.matmul(torch.from_numpy(t2), X_tensor).numpy()
            temp = Bs - aux_betas
            nesterov = np.divide((-np.exp(-temp/mu) + np.exp(temp/mu)),(np.exp(-temp/mu) + np.exp(temp/mu)))
            term2 = L1_penalty*nesterov
            term3 = term2
            val = (-1/Neff)*term1 + term3
            return val
          
      ans = sp.optimize.minimize(costfunc_log, jac=gradfunc_log, x0=Beta_hat_init, args=(L1_penalty),\
                                 method='L-BFGS-B', options={'maxiter':max_iter, 'ftol':tol,'maxfun':max_fun})
                                 
      final_snps = tar_snps
      if args.covTrain:
        final_betas = ans.x[numCov:(M+numCov)]
      else:
        final_betas = ans.x
      final_chr = tar_chr
      final_pos = tar_pos
      final_a1 = tar_a1
      final_results_df = pd.DataFrame({'CHR':final_chr, 'SNP':final_snps, 'POS':final_pos,\
                                       'A1':final_a1, 'BETA':final_betas,})
                                 
  return final_results_df

def doValidation(args):
  
  if not args.validate:
    sys.exit(f"Error: Validation cannot be peformed when --validate is False.")
    
  else:
    if args.trait_type == 'continuous':
      if args.version in ['MPP-PRS+','MPP-PRS','MPP-GWAS','MPP-GWAS-Admix']:
        L1_vals = args.L1.split(',')
        L1_vals = [float(i) for i in L1_vals]
        temp_path = args.out + '/temp'
        os.system(f"mkdir -p {temp_path}")
        best_L1 = 0
        best_R2 = 0
        for L1_penalty in L1_vals:
          num_chrs = args.chr.split(',')
          num_chrs = [int(i) for i in num_chrs]
          for chr_num in num_chrs:
            bfile_path = args.bfileVal + str(chr_num)
            score_path = args.out + '/mpp_predbetas_chr'+ str(chr_num)+'_L1_'+str(L1_penalty)+'.txt'
            out_path = temp_path +'/y_predbetas_chr'+ str(chr_num)+'_L1_'+str(L1_penalty)
            os.system(f"./plink2 --bfile {bfile_path} --score {score_path} 2 4 5 header --out {out_path} --silent")
            
          if args.trait_type == 'continuous':
            Y_tar_val = pd.read_csv(args.phenoVal, delimiter='\t')
            if args.covVal:
              Cov_val = pd.read_csv(args.covVal, delimiter='\t')
              temp_df = pd.concat([Y_tar_val[['Pheno']],Cov_val.iloc[:,2:]],axis=1)
              model = ols('Pheno ~ .', data=temp_df).fit()
              temp_y_hat = model.predict(temp_df)
              Y_tar_val = np.array(Y_tar_val['Pheno'], dtype=float) - temp_y_hat
            else:
              Y_tar_val = np.array(Y_tar_val['Pheno'], dtype=float)
            
            for i in range(len(num_chrs)):
              if i == 0:
                out_path = temp_path +'/y_predbetas_chr'+ str(num_chrs[i])+'_L1_'+str(L1_penalty)+'.sscore'
                Yhat_val_df = pd.read_csv(out_path, delimiter='\t')
                Yhat_val = np.array(Yhat_val_df['SCORE1_AVG'], dtype=float) * np.array(Yhat_val_df['ALLELE_CT'], dtype=float)
              else:
                out_path = temp_path +'/y_predbetas_chr'+ str(num_chrs[i])+'_L1_'+str(L1_penalty)+'.sscore'
                temp = pd.read_csv(out_path, delimiter='\t')
                temp = np.array(temp['SCORE1_AVG'], dtype=float) * np.array(temp['ALLELE_CT'], dtype=float)
                Yhat_val += temp
                
            this_R2 = (sp.stats.pearsonr(Y_tar_val, Yhat_val)[0])**2
            if this_R2 > best_R2:
              best_L1 = L1_penalty
              best_R2 = this_R2
        
        print("Validation complete! Best L1: "+str(best_L1)+" and Best R2: "+str(best_R2))
      
      elif args.version in ['MPP-PRS-TarSS','MPP-GWAS-TarSS']:
        L1_vals = args.L1.split(',')
        L1_vals = [float(i) for i in L1_vals]
        L2_vals = args.L2.split(',')
        L2_vals = [float(i) for i in L2_vals]
        temp_path = args.out + '/temp'
        os.system(f"mkdir -p {temp_path}")
        best_L1 = 0
        best_L2 = 0
        best_R2 = 0
        for L1_penalty in L1_vals:
          for L2_penalty in L2_vals:
            num_chrs = args.chr.split(',')
            num_chrs = [int(i) for i in num_chrs]
            for chr_num in num_chrs:
              bfile_path = args.bfileVal + str(chr_num)
              score_path = args.out + '/mpp_predbetas_chr'+ str(chr_num)+'_L1_'+str(L1_penalty)+'.txt'
              out_path = temp_path +'/y_predbetas_chr'+ str(chr_num)+'_L1_'+str(L1_penalty)
              os.system(f"./plink2 --bfile {bfile_path} --score {score_path} 2 4 5 header --out {out_path} --silent")
              
            if args.trait_type == 'continuous':
              Y_tar_val = pd.read_csv(args.phenoVal, delimiter='\t')
              if args.covVal:
                Cov_val = pd.read_csv(args.covVal, delimiter='\t')
                temp_df = pd.concat([Y_tar_val[['Pheno']],Cov_val.iloc[:,2:]],axis=1)
                model = ols('Pheno ~ .', data=temp_df).fit()
                temp_y_hat = model.predict(temp_df)
                Y_tar_val = np.array(Y_tar_val['Pheno'], dtype=float) - temp_y_hat
              else:
                Y_tar_val = np.array(Y_tar_val['Pheno'], dtype=float)
              
              for i in range(len(num_chrs)):
                if i == 0:
                  out_path = temp_path +'/y_predbetas_chr'+ str(num_chrs[i])+'_L1_'+str(L1_penalty)+'.sscore'
                  Yhat_val_df = pd.read_csv(out_path, delimiter='\t')
                  Yhat_val = np.array(Yhat_val_df['SCORE1_AVG'], dtype=float) * np.array(Yhat_val_df['ALLELE_CT'], dtype=float)
                else:
                  out_path = temp_path +'/y_predbetas_chr'+ str(num_chrs[i])+'_L1_'+str(L1_penalty)+'.sscore'
                  temp = pd.read_csv(out_path, delimiter='\t')
                  temp = np.array(temp['SCORE1_AVG'], dtype=float) * np.array(temp['ALLELE_CT'], dtype=float)
                  Yhat_val += temp
                  
              this_R2 = (sp.stats.pearsonr(Y_tar_val, Yhat_val)[0])**2
              if this_R2 > best_R2:
                best_L1 = L1_penalty
                best_L2 = L2_penalty
                best_R2 = this_R2
        
        print("Validation complete! Best L1: "+str(best_L1)+", Best L2:"+str(best_L2)+" and Best R2: "+str(best_R2))
    
    if args.trait_type == 'binary':
      if args.version in ['MPP-PRS+','MPP-PRS','MPP-GWAS','MPP-GWAS-Admix']:
        L1_vals = args.L1.split(',')
        L1_vals = [float(i) for i in L1_vals]
        temp_path = args.out + '/temp'
        os.system(f"mkdir -p {temp_path}")
        best_L1 = 0
        best_PRAUC = 0
        for L1_penalty in L1_vals:
          num_chrs = args.chr.split(',')
          num_chrs = [int(i) for i in num_chrs]
          for chr_num in num_chrs:
            bfile_path = args.bfileVal + str(chr_num)
            score_path = args.out + '/mpp_predbetas_chr'+ str(chr_num)+'_L1_'+str(L1_penalty)+'.txt'
            out_path = temp_path +'/y_predbetas_chr'+ str(chr_num)+'_L1_'+str(L1_penalty)
            os.system(f"./plink2 --bfile {bfile_path} --score {score_path} 2 4 5 header --out {out_path} --silent")
            
          if args.trait_type == 'binary':
            Y_tar_val = pd.read_csv(args.phenoVal, delimiter='\t')
            Y_tar_val = np.array(Y_tar_val['Pheno'], dtype=float)
            if args.covVal:
              Cov_val = pd.read_csv(args.covVal, delimiter='\t')
              Cov_val = Cov_val.iloc[:,2:]
            
            for i in range(len(num_chrs)):
              if i == 0:
                out_path = temp_path +'/y_predbetas_chr'+ str(num_chrs[i])+'_L1_'+str(L1_penalty)+'.sscore'
                Yhat_val_df = pd.read_csv(out_path, delimiter='\t')
                Yhat_val = np.array(Yhat_val_df['SCORE1_AVG'], dtype=float) * np.array(Yhat_val_df['ALLELE_CT'], dtype=float)
              else:
                out_path = temp_path +'/y_predbetas_chr'+ str(num_chrs[i])+'_L1_'+str(L1_penalty)+'.sscore'
                temp = pd.read_csv(out_path, delimiter='\t')
                temp = np.array(temp['SCORE1_AVG'], dtype=float) * np.array(temp['ALLELE_CT'], dtype=float)
                Yhat_val += temp
                
            if args.covVal:
              temp_df = np.concatenate((Cov_val.values,Yhat_val.reshape(-1, 1)), axis=1)
            else:
              temp_df = Yhat_val.reshape(-1, 1)
            temp_df = sm.add_constant(temp_df)
            glm_model = sm.GLM(Y_tar_val, temp_df, family=sm.families.Binomial())
            glm_result = glm_model.fit()
            y_final_pred = glm_result.predict(temp_df, linear=True)
            precision, recall, thresholds = precision_recall_curve(Y_tar_val, y_final_pred)
            this_PRAUC = auc(recall, precision)
            if this_PRAUC > best_PRAUC:
              best_L1 = L1_penalty
              best_PRAUC = this_PRAUC
        
        print("Validation complete! Best L1: "+str(best_L1)+" and Best PR-AUC: "+str(best_PRAUC))
        
  return 

def main():
    parser = setup_parser()
    args = parser.parse_args()
    
    # Run strict cross checking
    args = crosscheck_args(args)
    
    print_summary(args)
    print(f"All arguments are valid. Starting {args.version} pipeline...\n")
    
    np.random.seed(args.seed)
    
    if (args.Nthreads == False) or (args.Nthreads == 'all'):
      Nthreads = None
    else:
      Nthreads = int(args.Nthreads)
      
    if not args.validate:
      if args.version in ['MPP-PRS+','MPP-PRS','MPP-GWAS','MPP-GWAS-Admix']:
        L1_penalty = float(args.L1)
        L2_penalty = 0
      elif args.version in ['MPP-PRS-TarSS','MPP-GWAS-TarSS']:
        L1_penalty = float(args.L1)
        L2_penalty = float(args.L2)
    
      # Proceed with method execution logic
      num_chrs = args.chr.split(',')
      num_chrs = [int(i) for i in num_chrs]
      
      # List to hold the results
      results = [None] * len(num_chrs)
      process_with_args = partial(doMPP, L1_penalty=L1_penalty, L2_penalty=L2_penalty, args=args)
      
      print('Start Time: '+ str(datetime.now()))
      
      # Use ProcessPoolExecutor for parallelism
      with concurrent.futures.ProcessPoolExecutor(max_workers=Nthreads) as executor:
          # Submit tasks
          future_to_index = {
              executor.submit(process_with_args, chr_num): i for i, chr_num in enumerate(num_chrs)
          }
      
          # Collect results as they complete
          for future in tqdm(concurrent.futures.as_completed(future_to_index), total=len(num_chrs), desc="Processing Chromosomes"):
              index = future_to_index[future]
              try:
                  results[index] = future.result()
              except Exception as exc:
                  print(f"Iteration {index} generated an exception: {exc}")
      
      print('End Time: '+ str(datetime.now()))
      
      # Store output files
      os.system(f"mkdir -p {args.out}")
      for i in range(len(results)):
        results[i].to_csv(args.out + '/mpp_predbetas_chr'+ str(num_chrs[i])+'.txt', sep='\t', index=False, header=True)
        
    elif args.validate:
      L1_vals = args.L1.split(',')
      L1_vals = [float(i) for i in L1_vals]
      if args.version in ['MPP-PRS+','MPP-PRS','MPP-GWAS','MPP-GWAS-Admix']:
        L2_penalty = 0
        for L1_penalty in L1_vals:
          # Proceed with method execution logic
          num_chrs = args.chr.split(',')
          num_chrs = [int(i) for i in num_chrs]
          
          # List to hold the results
          results = [None] * len(num_chrs)
          process_with_args = partial(doMPP, L1_penalty=L1_penalty, L2_penalty=L2_penalty, args=args)
          
          print('Performing Validation using L1: '+ str(L1_penalty))
          
          # Use ProcessPoolExecutor for parallelism
          with concurrent.futures.ProcessPoolExecutor(max_workers=Nthreads) as executor:
              # Submit tasks
              future_to_index = {
                  executor.submit(process_with_args, chr_num): i for i, chr_num in enumerate(num_chrs)
              }
          
              # Collect results as they complete
              for future in tqdm(concurrent.futures.as_completed(future_to_index), total=len(num_chrs), desc="Processing Chromosomes"):
                  index = future_to_index[future]
                  try:
                      results[index] = future.result()
                  except Exception as exc:
                      print(f"Iteration {index} generated an exception: {exc}")
          
          # Store output files
          os.system(f"mkdir -p {args.out}")
          for i in range(len(results)):
            results[i].to_csv(args.out + '/mpp_predbetas_chr'+ str(num_chrs[i])+'_L1_'+str(L1_penalty)+'.txt',\
                              sep='\t', index=False, header=True)
                              
        doValidation(args)
                              
      elif args.version in ['MPP-PRS-TarSS','MPP-GWAS-TarSS']:
        for L1_penalty in L1_vals:
          for L2_penalty in L2_vals:
            # Proceed with method execution logic
            num_chrs = args.chr.split(',')
            num_chrs = [int(i) for i in num_chrs]
            
            # List to hold the results
            results = [None] * len(num_chrs)
            process_with_args = partial(doMPP, L1_penalty=L1_penalty, L2_penalty=L2_penalty, args=args)
            
            print('Performing Validation using L1: '+ str(L1_penalty)+' and L2: '+str(L2_penalty))
            
            # Use ProcessPoolExecutor for parallelism
            with concurrent.futures.ProcessPoolExecutor(max_workers=Nthreads) as executor:
                # Submit tasks
                future_to_index = {
                    executor.submit(process_with_args, chr_num): i for i, chr_num in enumerate(num_chrs)
                }
            
                # Collect results as they complete
                for future in tqdm(concurrent.futures.as_completed(future_to_index), total=len(num_chrs), desc="Processing Chromosomes"):
                    index = future_to_index[future]
                    try:
                        results[index] = future.result()
                    except Exception as exc:
                        print(f"Iteration {index} generated an exception: {exc}")
            
            # Store output files
            os.system(f"mkdir -p {args.out}")
            for i in range(len(results)):
              results[i].to_csv(args.out + '/mpp_predbetas_chr'+ str(num_chrs[i])+'_L1_'+str(L1_penalty)+'_L2_'+str(L2_penalty)+'.txt',\
                                sep='\t', index=False, header=True)
                                
        doValidation(args)

    print("MultiPopPred pipeline is complete!")

if __name__ == "__main__":
    main()
