# Description of the content of this repository

All R files must be executed with working directory set to the place they are stored.
The experiments were run on a parallel computer system and create intermediate files used in the later analysis.
If you want to avoid re-computation, we are happy to provide our computations upon request via email (the result files consume several GB and are therefore not hosted here).


## Overview of the workflow:
1. Creation and description of the data sets (see folder Data)
2. Conduction of experiments (see folder Experiments)
3. Evaluation of the experiments (see folders Descriptive results, Descriptive other and Random results)
4. For the resulting plots of which some are used in the paper see folder Plots


# Data:
- data sets used for the analysis
- data_preprocessing.R: creation of the three data sets as *.RData files
- data_properties: PCA plots used for displaying the properties of the data

## Descriptive other:
- evaluation of the experiments of the descriptive comparison of the stability measures based on the AP_Breast_Ovary and Stomach (results not shown in the paper)

## Descriptive results:
- evaluation of the experiments of the descriptive comparison of the stability measures based on AP_Colon_Kidney (results shown in the paper)
- stability_descriptive_results_evaluation.R: main file for the evaluation of the experiments
- ggplot_functions.R: contains plotting functions
- pareto_front_selection.R: epsilon constraint selection
- pareto_optimal.R: selects Pareto optimal models for each data set

## Experiments:
- experiments for the descriptive comparison of stability measures and the random search for a good model
- filter_methods.R: filter methods
- measures_of_stability.R: stability measures
- stability_descriptive.R: experiments for the descriptive comparison of stability measures
- stability_descriptive_results.R: collection of the results of the experiments for the descriptive comparison of stability measures
- stability_random.R: experiments for the random search for a good model
- stability_random_results.R: collection of the results of the experiments for the random search for a good model

## Plots:
- contains all plots that are displayed in the paper and additional plots
- all plots are created by stability_descriptive_results_evaluation.R in Deskriptive results

## Random results:
- evaluation of the experiments of the random search for a good model
- stability_random_results_evaluation.R: main file for the evaluation of the experiments
