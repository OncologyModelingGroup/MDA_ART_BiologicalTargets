%% Housekeeping
clear all; clc; close all; 

%Load data and set params
SpatialData = load("SpatialDistribution_BestScenrio.mat").VoxelLocalization_Linearized; 
n_voxel = 2;
r = 2; 
keys = double(table2array(SpatialData.ClassKeys.IDRLocalization)); 
idredge = keys(2); idrcore = keys(3); 
cases = string(fieldnames(SpatialData)); 
cases = cases(contains(cases, "Case_")); 

%% Output arrays
measured_dist = zeros(length(cases), 5); %[casenum, nvoxels_edge, nvoxels_core, pcratio]
model_dist = zeros(length(cases), 5); 

for c = 1:length(cases)
    casenum = double(strrep(cases(c), "Case_", "")); 
    obs_data = SpatialData.(cases(c)).tumorpattern_observed; 
    pred_data = SpatialData.(cases(c)).tumorpattern_predicted;
    
    %Get IDR tumor regions (obs and pred)
    obs_idr_edge = obs_data(obs_data(:,3) == idredge, 3); 
    obs_idr_core = obs_data(obs_data(:,3) == idrcore, 3);
    pred_idr_edge = pred_data(pred_data(:,3) == idredge,3); 
    pred_idr_core = pred_data(pred_data(:,3) == idrcore,3); 

    %Compute ratio
    o_nedge = length(obs_idr_edge); o_ncore = length(obs_idr_core); 
    p_nedge = length(pred_idr_edge); p_ncore = length(pred_idr_core); 
    obs_PCratio = o_nedge / o_ncore; 
    pred_PCratio = p_nedge / p_ncore; 

    %Subject agreement 
    if obs_PCratio >= r
       o_edgepreference = 1; 
    else
        o_edgepreference = 0; 
    end

    if pred_PCratio >= r
        p_edgepreference = 1; 
    else
        p_edgepreference = 0; 
    end

    %save to array
    measured_dist(c,:) = [casenum, o_nedge, o_ncore, obs_PCratio, o_edgepreference]; 
    model_dist(c,:) = [casenum, p_nedge, p_ncore, pred_PCratio, p_edgepreference]; 

end

%Compare subjects
obs_subjects_poscriteria = measured_dist(measured_dist(:,end) == 1, 1); 
pred_subjects_poscriteria = model_dist(model_dist(:,end) == 1, 1);
subj_agreement = sum(obs_subjects_poscriteria == pred_subjects_poscriteria, "all");
subj_mismatch = sum(obs_subjects_poscriteria ~= pred_subjects_poscriteria, "all"); 

disp("Number of subjects correctly identified with PC ratio >= 2:")
disp(subj_agreement)

disp("Number of subjects incorrectly identified with PC ratio >= 2:")
disp(subj_mismatch)

pct_correct = sum(pred_subjects_poscriteria(:,end)) / sum(obs_subjects_poscriteria(:,end)) * 100; 

%Tabulate and report results
t1 = table(measured_dist(:,1), measured_dist(:,2), measured_dist(:,3), measured_dist(:,4), measured_dist(:,5),...
    VariableNames=["CaseNum", "Num IDR Edge", "Num IDR Core", "IDR Periphery-to-Core Ratio", "IDR PC ratio > 2 (0-No, 1-yes)"]); 
t2 = table(model_dist(:,1), model_dist(:,2), model_dist(:,3), model_dist(:,4), model_dist(:,5), ...
    VariableNames=["CaseNum", "Num IDR Edge", "Num IDR Core", "IDR Periphery-to-Core Ratio", "IDR PC ratio > 2 (0-No, 1-yes)"]); 
SpatialAnalysis.Observed_IDRDistribution = t1; 
SpatialAnalysis.ModelPredicted_IDRDistribution = t2; 
SpatialAnalysis.Subjects_AboveThreshold_Obs= obs_subjects_poscriteria;
SpatialAnalysis.Subjects_AboveThreshold_Pred = pred_subjects_poscriteria;
SpatialAnalysis.Subjects_Agreement_Percent = pct_correct; 


%% save data
save("ReproducedResults_SpatialDistribution_BestScenario.mat", "SpatialAnalysis", "-v7.3")