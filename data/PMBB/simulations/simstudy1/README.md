# Notes 

This directory contains results from **Simulation Study 1**.

Unzip the object `correlations.zip` to obtain two files, whose details are provided below.

## Correlation between partial PGS and phenotype 

- This is `cor_parPGS_y.csv`. It has $720,000$ rows ($6\times6$ choices of `r2` and `rho`; $50$ seeds; $100$ reps per seed; $4$ quantiles)
- Quantity reported per column:
	- `a_bar`: Average global ancestry in PMBB quantile. This is seed-specific.
	- `r2`: Value of $r^2$ used in simulations.
	- `rho`: Value of $\rho$ used in simulations.
	- `Quantile`: PMBB quantile.
	- `Rep`: Rep number for a particular simulation seed.
	- `Seed`: Simulation seed.
	- `Global_Expected`: Expected correlation under the Global Model, as worked out in detail in the Supplement. Quantity is determined fully by `Seed`, `Quantile`, `r2` and `rho`. 
	- `Global_Observed`: Observed correlation under the Global Model. Quantity is rep-specific. 
	- `Local_Expected`: Expected correlation under the Local Model, as worked out in detail in the Supplement. Quantity is determined fully by `Seed`, `Quantile`, `r2` and `rho`.
	- `Local_Observed`: Observed correlation under the Local Model. Quantity is rep-specific. 
	- `Local_Analytical`: Analytical correlation under the Local Model, as presented in the Main Text and derived in the Supplement.
	- `Global_Analytical`: Analytical correlation under the Global Model, as presented in the Main Text and derived in the Supplement.

## Correlation between total PGS and phenotype 

- This is `cor_totPGS_y.csv`. It has $720,000$ rows ($6\times6$ choices of `r2` and `rho`; $50$ seeds; $100$ reps per seed; $4$ quantiles)
- Quantity reported per column:
	- `a_bar`: Average global ancestry in PMBB quantile. This is seed-specific.
	- `r2`: Value of $r^2$ used in simulations.
	- `rho`: Value of $\rho$ used in simulations.
	- `Quantile`: PMBB quantile.
	- `Rep`: Rep number for a particular simulation seed.
	- `Seed`: Simulation seed.
	- `Global_Expected`: Expected correlation under the Global Model, as worked out in detail in the Supplement. Quantity is determined fully by `Seed`, `Quantile`, `r2` and `rho`. 
	- `Global_Observed`: Observed correlation under the Global Model. Quantity is rep-specific. 
	- `Local_Expected`: Expected correlation under the Local Model, as worked out in detail in the Supplement. Quantity is determined fully by `Seed`, `Quantile`, `r2` and `rho`.
	- `Local_Observed`: Observed correlation under the Local Model. Quantity is rep-specific. 
	- `Analytical`: Analytical correlation under either the Global or Local Model, as presented in the Main Text and derived in the Supplement. Quantity is invariant to the model.
