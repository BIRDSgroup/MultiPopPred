# MultiPopPred

This is the official repository of the manuscript "MultiPopPred: A Trans-Ethnic Disease Risk Prediction Method, and its Application to the South Asian Population" by Ritwiz Kamal and Manikandan Narayanan.

- [License Preamble](#license-preamble)
- [Getting Started](#section-1-multipoppred---getting-started)
- [Using MultiPopPred](#section-2-using-multipoppred)
  - [Running Jupyter Notebooks](#running-jupyter-notebooks)
  - [Input Requirements](#input-requirements)
  - [Expected Outputs](#multipoppred-output)
  - [Example Data](#running-multipoppred-with-example-data)
- [Data Availability and Reproducibility](#section-4-data-availability-and-reproducibility)
- [MultiPopPred - Five Versions](#section-5-multipoppred---five-versions)
- [Support](#section-6-support)
- [Citation](#section-7-citation)

## License Preamble
Copyright 2024 BIRDS Group, IIT Madras

MultiPopPred is a free software: you can redistribute it and modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

MultiPopPred is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please take a look at the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with MultiPopPred. If not, see [https://www.gnu.org/licenses/](https://www.gnu.org/licenses/).

## Section 1: MultiPopPred - Getting Started

- The MultiPopPred Github repository can be downloaded using
   ```
   git clone https://github.com/BIRDSgroup/MultiPopPred.git
   ```
- The basic code for the five versions of MultiPopPred are provided as Jupyter notebooks in [Scripts](./Scripts).
- Example data for running the scripts is made available here: [Link to example data](https://1drv.ms/f/c/1d2cade3bfb64a9a/EtX4VV-0h1dJrXaW2kA4gWoBQxZX04sE25UcJ6V5jXccNA?e=LCRCHO).
- Hardware Requirements
  - MultiPopPred requires only a standard computer with enough RAM to support the in-memory operations.
- Software Requirements
   - OS Requirements
     - MultiPopPred is supported for MacOS and Linux. We have tested MultiPopPred on the following systems:
       - Linux: Ubuntu 18.04, 24.04 and CentOS 7.9.2009
       - MacOS: Ventura 13.5
   - Python Dependencies
     - Basic Python packages required to use MultiPopPred are available in [requirements.txt](./requirements.txt). Typical install time for each package on a "normal" desktop computer is ~2-5 minutes.
     - Optional: The complete conda environment can be reproduced using [environment.yml](./environment.yml).
       ```
       conda env create -f environment.yml
       ``` 
- Additional information such as Supplementary Data files associated with the manuscript are provided in [Supplementary Data](https://bit.ly/4nemWud) .

## Section 2: Using MultiPopPred

### Coming Soon

A command line version of MultiPopPred (executable as shown below) will be made available through this Github soon. Please watch this space for more updates.

```
python MultiPopPred-master.py --version MPP-PRS+ --aux_pops EUR,EAS,AMR,AFR --tar_pop SAS --aux_ss eur_ss.txt,eas_ss.txt,amr.txt,afr.txt --tar_ss sas_ss.txt --tar_geno training_geno --tar_pheno training_pheno --penalty 10 --out SAS_MPP_out.txt
```

### Running Jupyter Notebooks

The basic code for all five versions of MultiPopPred has been provided as Jupyter Notebooks (.ipynb files). To run these notebooks, one must have Jupyter installed on their systems which can easily be done following the instructions [over here](https://jupyter.org/install). Once Jupyter has been installed properly, Jupyter Notebook can be initialized as follows:

(a) For initializing Jupyter Notebook on local system
```
jupyter-notebook
```
(b) For initializing Jupyter Notebook on remote server

Step 1: From ternimal 1 enter the command
```
jupyter-notebook --no-browser
```
Step 2: From terminal 2 login to remote server using 
```
ssh -L 8080:localhost:xxxx username@hostname
```
(Here xxxx comes from the port number that you get after Step 1)

Step 3: From your browser, go to [http://localhost:8080/](http://localhost:8080/) and enter the token (Obtained after Step 1) for logging into jupyter notebook.

### Input Requirements

#### MPP-PRS+

The code for MPP-PRS+ is available at [MPP-PRS+.ipynb](./Scripts/MPP-PRS+.ipynb)

MPP-PRS+ requires the following inputs:
1. Outputs (from single-ancestry PRS analysis) from Lassosum-TrueLD for each auxiliary population (and optionally the target population).
   (A detailed account of the code, input and output requirements for Lassosum-TrueLD are provided later in this section)
   ```
     CHR	SNP	        POS	        A1	BETA
     22	rs5747999	17075353.0	A	0.00039767037558578256
     22	rs1807512	17221495.0	C	-0.0008640544683478455
   ```
2. Target Population's Genotype files for training, validation and testing datasets in [PLINK's](https://www.cog-genomics.org/plink/) bed/bim/fam format.
   
   2.1. The .bed file contains the binary encoded genotypes

   2.2. The .bim file contains information on SNPs in the following format
   ```
   22 rs5747999 0 17075353 A C
   22 rs1807512 0 17221495 C T
   ```

   2.3. The .fam file contains information on samples/individuals in the following format
   ```
   1001 1001 0 0 0 -9
   1002 1002 0 0 0 -9
   ```
   
3. Target Population's Phenotype files for training, validation and testing datasets in the following format
   ```
   -0.570771708057521
   -0.503465336753098
   0.0783526847346506
   1.08188259542337
   ```
4. A L1 penalty value.

#### Lassosum-TrueLD

The code for Lassosum-TrueLD is available at [Lassosum_TrueLD.ipynb](./Scripts/Lassosum_TrueLD.ipynb)

Lassosum-TrueLD works on a single population/anestry at a time and requires the following inputs:

1. The population in consideration's Genotype files for training, validation and testing datasets in [PLINK's](https://www.cog-genomics.org/plink/) bed/bim/fam format.
   
   1.1. The .bed file contains the binary encoded genotypes.

   1.2. The .bim file contains information on SNPs in the same format as depicted in point 2.2 under MPP-PRS+ input requirements.

   1.3. The .fam file contains information on samples/individuals in the same format as depicted in point 2.3 under MPP-PRS+ input requirements.
   
2. The population in consideration's Phenotype files for training, validation and testing datasets in the same format as depicted in point 3 under MPP-PRS+ input requirements.
   
3. A L1 penalty value.

#### MPP-PRS

The code for MPP-PRS is available at [MPP-PRS.ipynb](./Scripts/MPP-PRS.ipynb)

MPP-PRS requires the following inputs:

1. Outputs from Lassosum-ExtLD (or [Lassosum2 (Zhang et al., Nat. Comms., 2024)](https://github.com/Jingning-Zhang/PROSPER)) for each auxiliary population (and optionally the target population).
   (A detailed account of the code, input and output requirements for Lassosum-ExtLD are provided later in this section)
   ```
   rsid	        a1      a0      weight
   rs5747999	A	C	-0.011848444003009
   rs1807512	C	T	-0.0131474125666559
   ```
2. Target Population's Genotype files for training, validation and testing datasets in [PLINK's](https://www.cog-genomics.org/plink/) bed/bim/fam format.
   
   2.1. The .bed file contains the binary encoded genotypes

   2.2. The .bim file contains information on SNPs in the same format as depicted in point 2.2 under MPP-PRS+ input requirements.

   2.3. The .fam file contains information on samples/individuals in the same format as depicted in point 2.3 under MPP-PRS+ input requirements.
   
3. Target Population's Phenotype files for training, validation and testing datasets in the same format as depicted in point 3 under MPP-PRS+ input requirements.
   
4. A L1 penalty value.

#### Lassosum-ExtLD

The code for Lassosum-ExtLD (or Lassosum2) can be found at the Github repository of PROSPER ([Zhang et al., Nat. Comms., 2024](https://github.com/Jingning-Zhang/PROSPER)).

Lassosum-ExtLD works on a single population/anestry at a time and requires the following inputs (A more detailed account of Lassosum2's code, input and output requirements can be found on [PROSPER's Github](https://github.com/Jingning-Zhang/PROSPER)):

1. The population in consideration's GWAS summary statistics in the following format:
   ```
   rsid	        chr	a1	a0	beta	                beta_se	                n_eff
   rs5747999	22	A	C	-0.0046297381055207	0.0308228833408787	1000
   rs1807512	22	C	T	-0.0015843367342689	0.0308231909399726	1000
   ```
   
2. The population in consideration's Genotype files for training, validation and testing datasets in [PLINK's](https://www.cog-genomics.org/plink/) bed/bim/fam format.

   2.1. The .bed file contains the binary encoded genotypes.

   2.2. The .bim file contains information on SNPs in the same format as depicted in point 2.2 under MPP-PRS+ input requirements.

   2.3. The .fam file contains information on samples/individuals in the same format as depicted in point 2.3 under MPP-PRS+ input requirements.
       
3. The population in consideration's Phenotype files for training, validation and testing datasets merged into a single file in the following format:
   ```
   1001 1001 -0.570771708057521
   1002 1002 -0.503465336753098
   1003 1003 0.0783526847346506
   1004 1004 1.08188259542337
   ```
4. The population in consideration's external LD reference panel which can be obtained from [PROSPER's Github](https://github.com/Jingning-Zhang/PROSPER).

#### MPP-GWAS

The code for MPP-GWAS can be found at [MPP-GWAS.ipynb](./Scripts/MPP-GWAS.ipynb)

MPP-GWAS requires the following inputs:

1. GWAS Summary Statistics for each auxiliary population (and optionally the target population) in the following format:
   ```
   CHR   SNP          GENETIC.DIST       POS         A1   A2   BETA                   SE                    T                    P                 N
   22    rs5747999    0.0127800389155    17075353    A    C    -0.0046297381055207    0.0308228833408787    -0.15020457542272    0.880603216185275 1000
   22    rs1807512    0.153130271278     17221495    C    T    -0.0015843367342689    0.0308231909399726    -0.0514008019920572  0.959006145722472 1000
   22    rs4819535    0.154243623837     17278762    C    T    0.0316657163724304     0.0308069292177538    1.02787642833878     0.304007958825193 1000
   ```
2. Target Population's Genotype files for training, validation and testing datasets in [PLINK's](https://www.cog-genomics.org/plink/) bed/bim/fam format.
   
   2.1. The .bed file contains the binary encoded genotypes

   2.2. The .bim file contains information on SNPs in the following format
   ```
   22 rs5747999 0 17075353 A C
   22 rs1807512 0 17221495 C T
   ```

   2.3. The .fam file contains information on samples/individuals in the following format
   ```
   1001 1001 0 0 0 -9
   1002 1002 0 0 0 -9
   ```
   
3. Target Population's Phenotype files for training, validation and testing datasets in the following format
   ```
   -0.570771708057521
   -0.503465336753098
   0.0783526847346506
   1.08188259542337
   ```
4. A L1 penalty value.

#### MPP-GWAS-TarSS

The code for MPP-GWAS-TarSS can be found at [MPP-GWAS-TarSS.ipynb](./Scripts/MPP-GWAS-TarSS.ipynb)

MPP-GWAS-TarSS requires the following inputs:

1. GWAS Summary Statistics for each auxiliary population as well as the target population in the following format:
   ```
   CHR   SNP          GENETIC.DIST       POS         A1   A2   BETA                   SE                    T                    P                 N
   22    rs5747999    0.0127800389155    17075353    A    C    -0.0046297381055207    0.0308228833408787    -0.15020457542272    0.880603216185275 1000
   22    rs1807512    0.153130271278     17221495    C    T    -0.0015843367342689    0.0308231909399726    -0.0514008019920572  0.959006145722472 1000
   22    rs4819535    0.154243623837     17278762    C    T    0.0316657163724304     0.0308069292177538    1.02787642833878     0.304007958825193 1000
   ```
2. Target Population's external LD reference panel in [PLINK's](https://www.cog-genomics.org/plink/) bed/bim/fam format. The reference panel can be obtained from an appropriate source such as [1000 Genomes](https://cncr.nl/research/magma/). The computation of LD happens during runtime of the code.
   
3. L1 and L2 penalty values.

#### MPP-GWAS-Admixture

The code for MPP-GWAS-Admixture can be found at [MPP-GWAS-Admixture.ipynb](./Scripts/MPP-GWAS-Admixture.ipynb)

MPP-GWAS-Admixture requires the following inputs:

1. GWAS Summary Statistics for each auxiliary population (and optionally the target population) in the following format:
   ```
   CHR   SNP          GENETIC.DIST       POS         A1   A2   BETA                   SE                    T                    P                 N
   22    rs5747999    0.0127800389155    17075353    A    C    -0.0046297381055207    0.0308228833408787    -0.15020457542272    0.880603216185275 1000
   22    rs1807512    0.153130271278     17221495    C    T    -0.0015843367342689    0.0308231909399726    -0.0514008019920572  0.959006145722472 1000
   22    rs4819535    0.154243623837     17278762    C    T    0.0316657163724304     0.0308069292177538    1.02787642833878     0.304007958825193 1000
   ```
2. Target Population's Genotype files for training, validation and testing datasets in [PLINK's](https://www.cog-genomics.org/plink/) bed/bim/fam format.
   
   2.1. The .bed file contains the binary encoded genotypes

   2.2. The .bim file contains information on SNPs in the following format
   ```
   22 rs5747999 0 17075353 A C
   22 rs1807512 0 17221495 C T
   ```

   2.3. The .fam file contains information on samples/individuals in the following format
   ```
   1001 1001 0 0 0 -9
   1002 1002 0 0 0 -9
   ```
   
3. Target Population's Phenotype files for training, validation and testing datasets in the following format
   ```
   -0.570771708057521
   -0.503465336753098
   0.0783526847346506
   1.08188259542337
   ```
4. The admixture proportions for individuals in the target population in the following format:
   (The admixture proportions can be obtained with the help of the tool [Admixture](https://dalexander.github.io/admixture/))
   ```
   AUX1.PROP   AUX2.PROP   AUX3.PROP   AUX4.PROP
   0.663683    0.000010    0.243133    0.093174
   0.726047    0.000010    0.000010    0.273933
   0.512244    0.000010    0.000010    0.487736
   0.616460    0.021437    0.248714    0.113389
   ```
5. A L1 penalty value.

### MultiPopPred Output

All five versions of MultiPopPred produce output in the following format:
```
CHR	SNP	        POS	        A1	BETA
22	rs5747999	17075353.0	A	0.008025369947048135
22	rs1807512	17221495.0	C	0.0011822995349338257   
```

The PRS can then be computed with the improved target betas using either the evaluation function provided in the jupyter notebooks or using PLINK's score flag.

### Running MultiPopPred with example data

To run MultiPopPred on the provided [example data](https://1drv.ms/f/c/1d2cade3bfb64a9a/EtX4VV-0h1dJrXaW2kA4gWoBQxZX04sE25UcJ6V5jXccNA?e=LCRCHO), the following steps need to be followed:

1. Download example_data.tar.gz (~50 MB) from the link given above and decompress it using
   ```
   tar -xvzf example_data.tar.gz
   ```
2. Download MultiPopPred scripts from [Scripts](./Scripts).
3. Place the MultiPopPred script(s) within your example_data directory.
4. Modify the PATH variables within the script(s) as per your system.
5. Run the jupyter notebook script(s).

Typical expected runtime on a "normal" desktop computer is ~1-3 seconds for each of the five versions of MultiPopPred applied to the provided example data, for a given single hyperparameter configuration.

## Section 4: Data Availability and Reproducibility

### Data Availability

Genotype-Phenotype data generated in this work is available through the following links:

- Complete simulated genotype data for EUR, EAS, AMR, AFR, SAS populations are available as follows:
  - EUR: Coming Soon
  - EAS: Coming Soon
  - AMR: Coming Soon
  - AFR: Coming Soon
  - SAS: Coming Soon

- Complete simulated genotype-phenotype data for each of our simulation analyses are available as follows:
  - Simulation Analyses 1 (Varying heritability of simulated trait): Coming Soon
  - Simulation Analyses 2 (Varying number of samples in auxiliary populations): Coming Soon
  - Simulation Analyses 3 (Varying number of sample in target population): Coming Soon

### Reproducibility - Level 1: Reproducing Figures from Final Results

Data and code for reproducing the main text figures from final results are available in [Figures](./Figures)

### Reproducibility - Level 2: Reproducing Final Results from Raw Outputs

Coming Soon

### Reproducibility - Level 3: Reproducing Raw Outputs from Input Genotype-Phenotype Data

Coming Soon

### Reproducibility - Level 4: Reproducing Input Genotype-Phenotype Data

Coming Soon

## Section 5: MultiPopPred - five versions
The code for five versions of MultiPopPred, as illustrated in the figure below, are provided.

![Methodology Overview](https://github.com/BIRDSgroup/MultiPopPred/blob/main/Figures/Plots/method_github.png)

## Section 6: Support

For any queries regarding MultiPopPred, please contact Ritwiz Kamal (ritwiz@cse.iitm.ac.in) and/or Manikandan Narayanan (nmanik@cse.iitm.ac.in)

## Section 7: Citation

Kamal, R. and Narayanan, M., 2024. MultiPopPred: A Trans-Ethnic Disease Risk Prediction Method, and its Application to the South Asian Population. bioRxiv, pp.2024-11. [https://doi.org/10.1101/2024.11.26.625410](https://doi.org/10.1101/2024.11.26.625410)






