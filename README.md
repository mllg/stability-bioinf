# Description of the content of this repository

All R files must be executed with working directory set to the place they are stored.
The experiments were run on a parallel computer system and create intermediate files used in the later analysis.
If you want to avoid re-computation, we are happy to provide our computations upon request via email (the result files consume several GB and are therefore not hosted here).


## Overview of the workflow:
1. Creation and description of the data sets (see folder Data)
2. Conduction of experiments (see folders Experiments)
3. Evaluation of the experiments (see folders Evaluation)


## Data:
- data sets used for the analysis
- data_preprocessing.R: creation of the three data sets as *.RData files
- data_properties: PCA plots used for displaying the properties of the data

## Measures comparison:
### Experiments:
- measures_comparison_experiments.R: experiments for the empirical comparison of the stability measures
- measures_comparison_results.R: collection of the results of the experiments for the empirical comparison of the stability measures

### Evaluation:
- evaluation of the experiments of the empirical comparison of the stability measures
- measures_comparison_evaluation.R: main file for the evaluation of the experiments
- ggplot_functions.R: contains plotting functions
- pareto_optimal.R: selects Pareto optimal models for each data set

## Model finding:
### Experiments:
- model_finding_experiments.R: experiments for the random search for desirable configurations
- model_finding_results.R: collection of the results of the experiments for the random search for desirable configurations

### Evalutaion:
- evaluation of the experiments of the random search for desirable configurations
- model_finding_evaluation.R: main file for the evaluation of the experiments
