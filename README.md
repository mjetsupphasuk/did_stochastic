# Repository for "Difference-in-differences with stochastic policy shifts of continuous treatments"

This repository provides the code for the manuscript "Difference-in-differences with stochastic policy shifts of a continuous treatment" by Michael Jetsupphasuk, Chenwei Fang, Didong Li and Michael Hudgens [1]. 

## Directory structure

This Github repo contains the code and materials to re-create the results presented in the manuscript. This repo uses RStudio Projects and the `renv` package for reproducibility. Before beginning, run the file `code/preliminaries.R` to create the directory structure used in the rest of the code. 

`/code/`: contains the R scripts used for analysis and the data application.

`/data/`: contains simulation data (created by files in `/code/`) and the data for the main data application (see [2] to obtain the data).

`/results/`: results from the simulation and data application are stored here.

## Simulation

In order to re-create the simulation results found in the manuscript, run the files in the below list in order. Note that these files are meant to run on a computing cluster rather than a personal computer. In particular, `simulation.R` may be submitted with a job array to a computing cluster where each job corresponds to analyzing one simulated dataset. We have found that 5gb of memory and 20 minutes per job was sufficient. All R code files are found in `/code/`. The simulation with multiple time periods follows a similar structure where the relevant files are suffixed with "mtp", e.g., `sims_gen_mtp.R` generates the datasets for the multiple time periods simulation experiment. Parts of these code files are adapted from the code corresponding with [Schindl et al. (2026)](http://arxiv.org/abs/2409.11967) [2] and [Callaway and Sant'Anna (2021)](https://doi.org/10.1016/j.jeconom.2020.12.001) [3]. 

1. `sims_gen.R`: creates the simulated datasets.
2. `simulation.R`: estimates the target parameter in each simulated dataset. 
3. `sims_truth.R`: creates a very large simulated dataset to compute the true target parameter.
4. `sims_eval.R`: aggregates the results of `simulation.R` and reads in the true parameter from `sims_truth.R` to compute statistics (e.g., bias, coverage) evaluating the simulations.

The following are code files that are not directly run but are sourced by other files:

- `sims_scenarios.R`: specifies simulation parameters for the different scenarios.
- `utils.R`: functions used in other files.

The files used to create the illustrations of parallel trends and various stochastic interventions are:

- `exp_tilt_example.R`: illustrates various stochastic interventions.
- `parallel_trends.R`: demonstrates an example of parallel trends assumptions.

Other files:

- `sims_forloop.R`: simple for loop in order to run and evaluate sims without submitting to a computing cluster; practically only used to check performance of oracle nuisance function estimators or parametric model estimators due to computation time.

## Data application

Make sure the data is downloaded from Bartik et al. [4]\ (http://doi.org/10.1257/app.20170487) and placed in the appropriate directory (`data/fracking/`). In particular, the files `county_clean_long.dta`, `county_flat.dta`, and `hdpi_allwells_shale.dta` are used. In order to re-create the results of the paper, run the files in the below list in order. The file `fracking_desc.R` performs data cleaning, data exploration, and saves the dataset used in the analysis. Then, `fracking.R` is run to compute the analysis. Similar to the simulation, this file is meant to be run on a computing cluster with a job array passed in corresponding to a different "case", which defines the target parameter. Then, `fracking_eval.R` evaluates the results from `fracking.R`. 

1. `fracking_desc.R`: data cleaning and data exploration. This file calls in data from the US Census / ACS. Near the top, replace "INSERT YOUR KEY HERE" with your personal Census API key. You may request a key [here](https://api.census.gov/data/key_signup.html): https://api.census.gov/data/key_signup.html. 
2. `fracking.R`: performs the main analysis.
3. `fracking_eval.R`: aggregates and evaluates results from `fracking.R`.

## References

[1] Jetsupphasuk M, Fang C, Li D, Hudgens MG. Difference-in-differences with stochastic policy shifts of continuous treatments. arXiv:2512.00296. doi:10.48550/arXiv.2512.00296

[2] Schindl K, Shen S, Kennedy EH. Incremental effects for continuous exposures. arXiv. Preprint posted online January 27, 2026:arXiv:2409.11967. doi:10.48550/arXiv.2409.11967

[3] Callaway B, Sant’Anna PHC. Difference-in-Differences with multiple time periods. Journal of Econometrics. 2021;225(2):200-230. doi:10.1016/j.jeconom.2020.12.001

[4] Bartik AW, Currie J, Greenstone M, Knittel CR. The Local Economic and Welfare Consequences of Hydraulic Fracturing. American Economic Journal: Applied Economics. 2019;11(4):105-155. doi:10.1257/app.20170487
