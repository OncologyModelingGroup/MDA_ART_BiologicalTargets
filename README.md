# MDA Adaptive RT - Biological Targets
This README provides an overview of the MATLAB scripts and data files used for the classification and spatial analysis of Apparent Diffusion Coefficient (ADC) data. These scripts are designed to identify tumor regions of Increased Diffusion Restriction (IDR) versus Decreased/Maintained Diffusion Restriction (DMDR).

## Setup and Requirements
Dependencies
MATLAB (Recommended R2022b or later).
Statistics and Machine Learning Toolbox (Required for ecdf, paretotails, and confusionmat).

## Data Description
| Filename | Description |
|----------| ------------|
|ADC_Data.mat | MATLAB file containing each subject's voxel-wise measurement of apparent diffusion coefficient within the tumor. |
|SpatialDistribution_BestScenario.mat | MATLAB file containing each tumor voxel's class (IDR vs. DMDR) using the optimal cutoff provided in the manuscript.

## Code Summary
| Filename | Description |
|----------| ------------|
| Script1_VoxelClassification_OptimalCutoffResults.m |Performs voxel-wise classification using pre-determined "optimal" cutoffs. Generates confusion matrices, eCDF plots, and performance metrics (Sensitivity, Specificity, DSC). |
|Script2_VoxelClassification_AllCutoffROC.m | Iterates through multiple LTCP (Lower Tail Cumulative Probability) thresholds to perform ROC analysis and determine the best-performing cutoffs for the model. 
|Script3_SpatialDistributionAnalysis_IDR.m | Analyzes the localization of IDR voxels, specifically comparing the "Periphery-to-Core" (PC) ratio between observed and predicted data. |

## Methodology: Voxel-Wise Tumor Classification and Spatial Analysis
This repository contains a suite of MATLAB scripts for processing ADC imaging data, determining optimal classification thresholds via ROC analysis, and analyzing the spatial distribution of predicted tumor patterns.

#### Statistical Mapping & Classification
The script utilize Empirical Cumulative Distribution Functions (eCDF) and Pareto Tail Analysis to establish thresholds for $\Delta ADC$ ($ADC_{final} - ADC_{baseline}$).

| Class | Description |
|-------| ------------|
|IDR (Class 1) | Areas where $\Delta ADC$ is $\le$ the determined threshold (increased restriction). |
|DMDR (Class 2) | Areas where $\Delta ADC$ is > the threshold (maintained/decreased restriction). |
 
#### ROC Analysis
To find the optimal balance between sensitivity and specificity, the analysis tests LTCP thresholds ranging from 20% to 40%. Observed Cutoff: Selected based on the highest Area Under the Curve (AUC). Model Cutoff: Selected using the Youden’s Index ($J$), where $J = Sensitivity + Specificity - 1$.

#### Performance Evaluation
The scripts calculate the following metrics to validate the model against observed data: Sensitivity (TPR) and Specificity (TNR). Dice-Sorensen Coefficient (DSC): Measures spatial overlap/agreement between predicted and observed classes. Confusion Matrices: Generated for each patient case across four different scenarios.

| Scenario | Tumor Regions | Timepoint |
|-----|-------|-----|
| Scenario 1 | Enhancing + Non-Enhancing Tumor Volume | End of radiotherapy |
| Scenario 2 | Enhancing + Non-Enhancing Tumor Volume | 1-month follow-up visit |
| Scenario 3 | Enhancing Tumor Volume | End of radiotherapy |
| Scenario 4 | Enhancing Tumor Volume | 1-month follow-up visit |

#### Spatial Distribution
This analysis classifies IDR voxels based on their location within the tumor: IDR Edge (Periphery) vs. IDR Core. PC Ratio: A ratio $\ge 2$ is used as a criterion to identify tumors with a preference for peripheral IDR localization.

| Class Number  | Description |
|---------------|-------------|
| 0 | DMDR voxel |
| 1 | IDR voxel located in the tumor periphery |
| 2 | IDR voxel located in the tumor core |
