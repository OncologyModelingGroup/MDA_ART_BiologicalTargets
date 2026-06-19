%% Housekeeping
clear all; clc; close all; 

%Load data and set params
SpatialData = load("SpatialDistribution_BestScenario.mat").VoxelLocation_Linearized; 
keys = double(table2array(SpatialData.ClassKeys.IDRLocation)); 
DiffKey = double(table2array(SpatialData.ClassKeys.DiffusionPattern));
idr_key = DiffKey(1); 
dmdr_key = DiffKey(2); 
cases = string(fieldnames(SpatialData)); 
cases = cases(contains(cases, "Case_")); 

%% Output arrays
spatial_dist = zeros(length(cases), 7); %[casenum, obs_nvoxels_idr_edge, obs_nvoxels_idr_core, pred_nvoxels_idr_edge, pred_nvoxels_idr_core, pcratio_obs, pc_ratio_pred]
check_array = zeros(length(cases), 2); 
rep_array = zeros(length(cases), 9); 


for c = 1:length(cases)
    casenum = double(strrep(cases(c), "Case_", "")); 
    obs_data = double(table2array(SpatialData.(cases(c)).tumorpattern_observed)); 
    pred_data = double(table2array(SpatialData.(cases(c)).tumorpattern_predicted));

    %Get location index
    idx_edge_voxels = SpatialData.(cases(c)).edge_location_index;
    idx_core_voxels = SpatialData.(cases(c)).core_location_index; 

    % Split observed class
    obs_idr = obs_data(obs_data(:,2) == 1, :); 
    obs_dmdr = obs_data(obs_data(:,2) == 2, :); 
    pred_idr = pred_data(pred_data(:,2) == 1, :); 
    pred_dmdr = pred_data(pred_data(:,2) == 2, :); 

    %Stratify data by location
    obs_core = obs_data(idx_core_voxels, :); 
    obs_edge = obs_data(idx_edge_voxels, :);  
    pred_edge = pred_data(idx_edge_voxels,:); 
    pred_core = pred_data(idx_core_voxels,:); 

    %Stratify data by diffusion pattern within location (obs and pred)
    obs_edge_idr = obs_idr(ismember(obs_idr(:,1), idx_edge_voxels), :); 
    obs_core_idr = obs_idr(ismember(obs_idr(:,1), idx_core_voxels), :);  
    pred_edge_idr = pred_idr(ismember(pred_idr(:,1), idx_edge_voxels), :); 
    pred_core_idr = pred_idr(ismember(pred_idr(:,1), idx_core_voxels), :); 

    %Get number of IDR elements by spatial location
    [n_pred_edge_idr, ~] = size(pred_edge_idr); 
    [n_pred_core_idr, ~] = size(pred_core_idr); 
    [n_obs_edge_idr, ~] = size(obs_edge_idr); 
    [n_obs_core_idr, ~] = size(obs_core_idr); 
    [n_obs_total_idr, ~] = size(obs_idr); 
    gt_IDR_count = length(obs_idr); 
    pred_IDR_count = length(pred_idr); 

    %Compute ratio
    pred_PC_IDR = (n_pred_edge_idr /gt_IDR_count) / (n_pred_core_idr /gt_IDR_count);
    obs_PC_IDR = (n_obs_edge_idr / gt_IDR_count) / (n_obs_core_idr / gt_IDR_count); 

    %save to array
    spatial_dist(c,:) = [casenum, n_obs_edge_idr, n_obs_core_idr, n_pred_edge_idr, n_pred_core_idr, obs_PC_IDR, pred_PC_IDR]; 
end

%Compare obs and predicted Periphery to Core Ratio (IDR only)
[ccc] = ccc_calc(spatial_dist(:,6), spatial_dist(:,7)); 
disp("Concordance Correlation Coefficient for Periphery-to-Core Ratio (obs vs pred)")
disp(ccc)

%Tabulate and report results
t1 = table(spatial_dist(:,1), spatial_dist(:,2),spatial_dist(:,3),spatial_dist(:,4), spatial_dist(:,5),spatial_dist(:,6),spatial_dist(:,7), VariableNames=["CaseNum", "Num IDR Edge (Obs)", "Num IDR Core (Obs)", "Num IDR Edge (Pred)", "Num IDR Core (Pred)", "Observed IDR Periphery-to-Core Ratio", "Predicted IDR Periphery-to-Core Ratio"]); 
SpatialAnalysis.IDR_SpatialDistribution = t1; 
SpatialAnalysis.ConcordanceCorrelationCoefficient = ccc; 

%% Create Supplemental Figure
PC_array = spatial_dist(:,[1,6,7]); 

%Scatterplot - Compare Obs vs. Pred Periphery-to-Core (PC) Ratio
f1 = figure; 
figure(f1)
scatter(PC_array(:,2), PC_array(:,3), 90, "filled", "o", "MarkerEdgeColor",[0.3686, 0.3686, 0.3686], "MarkerFaceColor",[0.3686, 0.3686, 0.3686])
hold on
xlabel("Observed", FontSize=18, FontName = "Times New Roman")
ylabel("Predicted", FontSize=18, FontName="Times New Roman")
title("IDR Periphery-to-Core Ratio" ,FontSize=20, FontName="Times New Roman")
xlim([0,2.5])
ylim([0,2.5]) 
plot(xlim, ylim, 'k--', 'LineWidth', 1.5, 'Color',"k") 
fontsize(f1, 24, "points")
fontname(f1, "Times New Roman")

%Bar Chart - Case by case report of PC Ratio
f2 = figure; 
figure(f2)
labels = ["1", "2", "3", "4", "5", "6", "8", "9", "10", "11", "14", "16", "17", "18", "20", "21", "23", "24", "26", "27", "28"]; 
face_color = [0.9686, 0.8941 , .8471; ...
            0.9137, .4431,  0.1961]; 
b = bar(labels, PC_array(:,2:3), "FaceColor","flat"); 
for k = 1:size(PC_array(:,2:3),2)
    b(k).CData =face_color(k,:);
end 
xlabel("Case ID")
ylabel("IDR Periphery-to-Core Ratio")
fontsize(f2, 24, "points")
fontname(f2, "Times New Roman")
legend({'Observed', 'Predicted'}, 'Location', 'northeast');


%% save data
save("ReproducedResults_SpatialDistribution_BestScenario.mat", "SpatialAnalysis", "-v7.3")
saveas(f1, "Scatter_IDR-SpatialDistribution.fig")
close(f1)
saveas(f2, "BarChart_IDR-SpatialDistribution.fig")
close(f2)

%% Functions
function [ccc] = ccc_calc(x,y)
%{Input: array of same dimensions

% ccc = concordance correlation coefficient
R = corrcoef(x,y); 
pcc = R(2); % Pearson correlation coefficient

x = x(:); y = y(:);
n = length(x);
u1 = mean(x);
u2 = mean(y);
sx = 0;
sy = 0;
for i = 1:n
  sx = sx+(x(i)-u1)^2;
  sy = sy+(y(i)-u2)^2;
end
sx = sqrt(sx/n);
sy = sqrt(sy/n);
sig1 = sx;
sig2 = sy;
sig12 = 0;

for i=1:n
    sig12 = sig12+((x(i)-u1)*(y(i)-u2));
end
sig12 = sig12/n;
ccc = 2*pcc*sig12/(sig1^2+sig2^2+(u1-u2)^2);
end