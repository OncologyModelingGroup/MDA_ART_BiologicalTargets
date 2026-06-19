# MDA Adaptive RT - Biological Targets
This README provides an overview of the MATLAB scripts and data files used for the classification and spatial analysis of Apparent Diffusion Coefficient (ADC) data. These scripts are designed to delineate candidate radiotherapy targets informed by changes in diffusion within the tumor. Voxels are classified as Increased Diffusion Restriction (IDR) or Decreased/Maintained Diffusion Restriction (DMDR) in this framework.

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
|Script3_IDR_SpatialDistributionAnalysis_EnhBorder.m | Analyzes the localization of IDR voxels, specifically comparing the "Periphery-to-Core" (PC) ratio between observed and predicted data. |

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
The scripts calculate the following metrics to compare the model against observed data: Sensitivity (TPR) and Specificity (TNR). Dice-Sorensen Coefficient (DSC): Measures spatial overlap/agreement between predicted and observed classes. Confusion Matrices: Generated for each patient case across four different scenarios.

| Scenario | Tumor Regions | Timepoint |
|-----|-------|-----|
| Scenario 1 | Enhancing + Non-Enhancing Tumor Volume | End of radiotherapy |
| Scenario 2 | Enhancing + Non-Enhancing Tumor Volume | 1-month follow-up visit |
| Scenario 3 | Enhancing Tumor Volume | End of radiotherapy |
| Scenario 4 | Enhancing Tumor Volume | 1-month follow-up visit |

#### Spatial Distribution
This analysis classifies IDR voxels based on their location within the tumor. Clinician delineated border between the enhancing and non-enhancing volume was used as the boundary to define tumor periphery (non-enhancing volume) versus tumor core (enhancing volume). IDR Periphery-to-Core Ratio was calculated across subjects and used to assess preferential localization.

| Class Number  | Description | Clinician Delineated Tumor Volume |
|---------------|-------------|------------------------------------|
| 1 | IDR voxel located in the tumor periphery | Non-enhancing volume |
| 2 | IDR voxel located in the tumor core | contrast enhancing volume |
