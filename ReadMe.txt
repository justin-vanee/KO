Description of data files and code associated with "Data reconciliation in multi-trait experiments with kinship ordination."

The algorithms contains five .stan model statements for our kinship ordination model (KO) as well as unconstrained ordination (UO), constrained ordination (CO), concurrent ordination (CCO), and the multivariate mixed model (MMM). In addition, there is also functions.R which is a collection of functions referenced in code.RMD for simulating data, fitting models, and extracting performance. 

Outputs is the default save location for all temporary R objects produced in code.RMD. The cleaned data for the Brooms tectorum analysis is preloaded in this folder for convenience. 

Figures is the daft save location for all figures created in code.RMD. 

Data contains the raw data for these analyses. Note that these files are also housed on the GitHub of Peter B. Adler and Megan L. Vahsen, and we recommend checking those pages for all data associated with BromeCast. 

The main code file for reproducing the results described in "Data reconciliation in multi-trait experiments with kinship ordination" is code.RMD. The code is delineated into various sections and subsections 

# Set up 

Lists all packages used, sets working directory, and load user defined functions in functions.R 

# Toy Example 

This simulates a dataset from the multivariate mixed model defined in MMM.stan. The simulated dataset can have up to 284 unique genotypes based on Q1 cross in mice described in more detail in the qtl2 package. Code is provided for fitting all five described models to the simulated dataset and plotting performance.   

# Heatmap of kinship 

Creates a Heatmap of kinship matrix for genotypes in mouse cross. 

# Prior Predictive Checks Heritability 

Plots implied prior distributions on heritability for a varying number of latent variables in our proposed kinship ordination model. 

# Data Wrangling (Cheatgrass: Growth Chamber + Common Gardens + Physiology, for viewing only can provide raw data upon request)

This calls various raw data files in Data and combines them into "data_list.RData". The list includes the scaled, centered, and potentially log-transformed trait data as well as the environmental covariates and kinship matrix.

# Parallel Simulation Study 

Broken into three subsections "Simulate Datasets," "Run Simulation Study in Parallel," and "Make Simulation Study Plots," which are self-explanatory. 

# Real Data Analysis (Cheatgrass: Growth Chamber + Common Gardens + Physiology)

This contains all the code for reproducing our data reconciliation analysis of Bromus tectorum (cheatgrass). 

## Data Preparation (need to run for all scripts below)

This create the various lists for fitting the described models using Stan. The multivariate mixed model cannot be fit to the full-dataset so separate lists are created for each experiment. As indicated, must sections below this one depend on these lists or the data frames themselves, so need to run it first. 

## Fit Models  

Fits the KO and MMM models and checks convergence with histograms of Rhat values (potential scale reduction factor). Note that it can take around 1-2 days to fit the MMM to the common garden experiment on an Apple M3 Max machine with 128 GB ram. It takes around 90 minutes to fit the KO model.   

## Heritability Plot 

Combines the heritability results for the four model fits (KO fit to all data, MMM fit separately to growth chamber, common gardens, physiology experiments) and compares estimates in terms of precision and congruence. Create one plot of all heritability estimates and intervals for both methods. 

## Trait Selection Plot

Calls the posterior predictive distributions for trait means of the 93 genotypes planted in the common gardens calculates the association (OLS estimate) between various traits and fitness for each MCMC iteration. Then creates a plot of selection coefficients. 

## 10-fold cross-validation

Again, split into three self-explanatory sections. The study takes a long time (even in parallel) because of the MMM. Without the MMM, it can finish in around 2 hours. The folds for cross validation are stratified by experiment to prevent all physiology traits, for example, from appearing in one fold. 