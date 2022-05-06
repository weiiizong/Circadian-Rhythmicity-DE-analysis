# Circadian Analysis: Rhythmicity detection and DE analysis
The purpose of this software is to demonstrate rhythmicity detection of circadian patterns in gene expression profiles as well as differential expression (DE) analysis using a linear model with gene specific covariate selection. Full analysis using this code can be found in the paper: Diurnal alterations in gene expression across the human dorsal and ventral striatum in psychosis (Ketchesin et al., in review)

## Required R packages: 
- minpack.lm v1.2-1
- Snowfall v1.84-6.1
- doParallel v1.016
- AWfisher https://github.com/Caleb-Huo/AWFisher
- *This code has been tested on R version 4.0.2 

## Initialization Instructions: 
Download and unzip the folders to the your working directory. 
Save each folder to the desktop
code/: code for analysis demostration
ouput/: results generated from code in code/
Install required packages listed above.

## Clinical Data (Example_clinical.csv):
72 subjects (36 psychosis + 36 matched controls). Clinical variables of interest include:
- pair: sample names matched to data columns.
- CorrectedTOD: corrected time of death (TOD) in ZT. 
- Diagnosis: "PSYCHOSIS"-psychosis subjects; "CONTROL"-matched control subjects.
- Sex.Label: sample sex, used as a covariate candidate in DE analysis. 
- Age: sample sample age, used as a covariate candidate in DE analysis. 
- Race.Label: sample race, used as a covariate candidate in DE analysis. 
- RIN_C: RIN value in caudate, used as a covariate candidate in DE analysis. 

## Expression Data (Example_data.csv)::  
- 15041 genes (all genes that were used in the full analysis of caudate region in the paper)
- Data has already been filtered according to method discussed in the paper 
- Expression units are in log2 CPM 


## Demo:
- fitSinCurve.R and Curve_Drawing_axis.R are external functions required to run the rhythmicity analysis.
- bestSelection.R is an external function required to run the DE analysis
- Run all external functions, so that they are in the R environment
- Run code as is in script.R. It contains two main sections: rhythmicity analysis and DE analysis. Three subsections are included in the rhythmicity analysis: Rhythmicity analysis on psychosis subjects, Rhythmicity analysis on control subjects and Change of rhythmicity analysis.
*NOTE: Each method requires permutation. Permutations were set to 10 to save time and computation space. Permuted files are included in the DE (output/DE_output/DE_null) and Rhythmicity (output/Rhythmicity_null_control, output/Rhythmicity_null_psychosis and output/ChangeRhyth_null_psychosis_control) folders to save reviewers time. If permutation files are recreated by the reviewer, output will be slightly different due to the randomness of permutation.


## Expected output:
- output/observed_para_control.csv: the example output of curve fitting parameters and p-values. 
- output/CorePlots: PDF plots of selective core circadian genes.
- output/ChangeRhyth_psychosisVScontrol.csv: the example output of differential rhythmicity parameters (e.g., R^2 (control) - R^2 (psychosis)) and p-values.
- output/DE_output/DE_psychosisVScontrol.csv: the example DE output including original p-value from the likelihood ratio test, corrected p-value and Benjamini Hochberg corrected p-value.
