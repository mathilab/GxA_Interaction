# Notes 

This directory contains results from **Simulation Study 2**.

## Local Ancestry Average Correlation of Effect Sizes

The result files are in the zipped file `laacor.zip`. Unzip `laacor.zip` to obtain six CSV files corresponding to the six phenotypes studied. Details of each file are described below.

- Each file is named `$Phenotype_all_LAACor_results.csv`, where `$Phenotype` is specified
- Each file has $25,000$ rows ($25$ choices of `Rho`; $10$ choice of `R2_Causal`; $100$ effect vector reps) 
- Quantity reported per column:
	- `R2_Causal`: Value of $r^2$ used in simulations.
	- `Rho`: Value of $\rho$ used in simulations.
	- `Global_LAACor_prime_hat`: Empirical correlations of local ancestry average causal effects under the Global Model. Quantity is rep-specific.
	- `Mean_Global_LAACor_prime_j_True`: Analytical approximation of the local ancestry average causal effect parameter under the Global Model, $\overline{\text{LAACor}'_\text{Glo}}$ (not dependent on the rep). See Supplementary Material Box C.
	- `Local_LAACor_prime_hat`: Empirical correlations of local ancestry average causal effects under the Local Model. Quantity is rep-specific.
	- `Mean_Local_LAACor_prime_j`: Analytical approximation of the local ancestry average causal effect correlation parameter under the Local Model, $\overline{\text{LAACor}'_\text{Loc}}$ (not dependent on the rep). See Supplementary Material Box C.
	- `Global_LAACor_hat`: Empirical correlations of local ancestry average tagging effects under the Global Model. Quantity is rep-specific.        
	- `Mean_Global_LAACor_j_True`: Analytical approximation of the local ancestry average tagging effect correlation parameter under the Global Model, $\overline{\text{LAACor}_\text{Glo}}$ (not dependent on the rep). See Supplementary Material Box D.
	- `Local_LAACor_hat`: Empirical correlations of local ancestry average tagging effects under the Local Model. Quantity is rep-specific.   
	- `Mean_Local_LAACor_j`: Analytical approximation of the local ancestry average tagging effect correlation parameter under the Local Model, $\overline{\text{LAACor}_\text{Loc}}$ (not dependent on the rep). See Supplementary Material Box D.
	- `a_bar_tag`: Sample average global ancestry computed across all tagging variants.
	- `a_bar_causal`: Sample average global ancestry computed across all causal variants.
	- `nvars`: Number of variants included in trait simulation.

## Polygenic Score Performance 

The result files are in the zipped file `batch_sim.zip`. Unzip `batch_sim.zip` to obtain four subdirectories of files. Each subdirectory comprises a batch of simulations:
- `batch_sim/batch1`: Simulations with $\rho\in\{0.2,0.35,0.5,0.65,0.8\}$
- `batch_sim/batch2`: Simulations with $\rho\in\{0.92,0.94,0.96,0.98,1\}$
- `batch_sim/batch3`: Simulations with $\rho\in\{0.91,0.93,0.95,0.97,0.99\}$
- `batch_sim/batch4`: Simulations with $\rho\in\{0.25,0.3,0.4,0.45,0.55,0.6,0.7,0.75,0.85,0.9\}$

Under each batch subdirectory, results for each phenotype and choice of Local or Global model are contained in a csv file ($6\times2=12$ files in total). Details of each file are provided below.

## Phenotype and GxA Model Result

- This is `$Phenotype_$Model_model_all_diagnostic_results.csv`, where `$Phenotype` and `$Model` are specified
- All `batch4` files have $20,000$ rows ($10\times10$ choices of `R2_Causal` and `Rho`; $200$ reps)
- All other batch subdirectory files have $10,000$ rows ($10\times5$ choices of `R2_Causal` and `Rho`; $200$ reps)
- Quantity reported per column:
	- `var_y`: Empirical variance of phenotype across sample. Quantity is rep-specific.
	- `var_totPGS_demean`: Empirical variance of Total PGS across sample. Quantity is rep-specific.
	- `var_parPGS_demean`: Empirical variance of Partial PGS across sample. Quantity is rep-specific.
	- `cor2_y_G`: Squared empirical correlation of phenotype and true genotype risk score. Quantity is rep-specific. 
	- `cor2_y_totPGS_demean`: Squared empirical correlation of phenotype and Total PGS. Quantity is rep-specific. 
	- `cor2_y_parPGS_demean`: Squared empirical correlation of phenotype and Partial PGS. Quantity is rep-specific.
	- `E_cor2_y_totPGS_demean`: Analytical approximation of the expected squared correlation of phenotype and Total PGS (not dependent on the rep). 
	- `E_cor2_y_parPGS_demean`: Analytical approximation of the expected squared correlation of phenotype and Partial PGS (not dependent on the rep). 
	- `R2_Causal`: Value of $r^2$ used in simulations.
	- `Rho`: Value of $\rho$ used in simulations.
	- `a_bar_tag`: Sample average global ancestry computed across all tagging variants.
	- `a_bar_causal`: Sample average global ancestry computed across all causal variants.
	- `cor_beta_causal_eur_tag`: Empirical correlation of European-ancestry causal effects and European-ancestry tagging effects. Quantity is rep-specific. 
	- `cor_beta_causal_afr_tag`: Empirical correlation of African-ancestry causal effects and African-ancestry tagging effects. Quantity is rep-specific. 
	- `cor_eur_tag_afr_tag`: Empirical correlation of African-ancestry tagging effects and European-ancestry tagging effects. Quantity is rep-specific. 
	- `nvars`: Number of variants included in trait simulation.