%% Housekeeping
clear all; clc; close all; warning('off','all')
ti_execute = datetime("now"); 
ADC_data = load("ADC_Data.mat").ADC_Data;
subjects = string(fieldnames(ADC_data.Scenario_1));
allsubjects = double(strrep(subjects, "Case_", ""));

%% Classification
scenarios = string(fieldnames(ADC_data));
scenarios = scenarios(contains(scenarios, "Scenario")); 
thresholds = [20, 25, 30, 35, 40]; 

%Compute the dADC for each scenario using optimal cutoffs
all_results = struct(); 
classification_results = struct(); 
pareto_results = struct(); 
ecdf_results = struct(); 
population_performance = struct(); 
subjectlevel_performance = struct(); 
roc_aggregated_vars = zeros(4*length(thresholds)^2,6); % scenario, obs cut, model cut, mean sensitivity, mean specificity, mean FPR;
ctr = 0; 
runtimes = zeros(4*(length(thresholds)^2), 3); %obs cut, model cut, run time (seconds) 

for ss = 1:length(scenarios)
    
        %Get optimal lower tail cumulative probability cutoff determined in ROC
        scene = strrep(scenarios(ss), "_", ""); 
        
        for ocut = 1:length(thresholds)
            ocut_i = thresholds(ocut); 
            ocut_tag = strcat("ObsCutoff_", num2str(ocut_i));
            
            for mcut = 1:length(thresholds)
                mcut_i = thresholds(mcut);
                mcut_tag = strcat("ModelCutoff_", num2str(mcut_i)); 
                ctr = ctr+1; 
            
                dadc_threshold_measured = ocut_i; 
                dadc_threshold_model = mcut_i; 
            
                %Get ADC data for each patient in all scenarios
                scene_cases = string(fieldnames(ADC_data.(scenarios(ss)))); 
                metrics_array = zeros(length(scene_cases), 8);
            
                %Run classification
                for cc = 1:length(scene_cases)
                    ID = scene_cases(cc);
                    id_num = double(strrep(ID, "Case_", "")); 
                    subj_data = table2array(ADC_data.(scenarios(ss)).(scene_cases(cc)).ADC); 
                    baseline_adc = subj_data(:,1); 
                    measured_adc_tf = subj_data(:,2); 
                    predicted_adc_tf = subj_data(:,3);
                
                    %Calculate dADC
                    if size(baseline_adc) == size(measured_adc_tf)
                       dadc_measured = measured_adc_tf - baseline_adc;
                    else
                        disp("Measured ADC at final timepoint and baseline data have different dimensions")
                        disp(ID)
                    end 
                    
                    if size(baseline_adc) == size(predicted_adc_tf)
                       dadc_predicted = predicted_adc_tf - baseline_adc;
                    else
                        disp("Predicted ADC at final timepoint and baseline data have different dimensions")
                        disp(ID)
                    end 
            
                    %eCDF fitting
                    [obs_y, obs_x, obs_LB, obs_UB] = ecdf(dadc_measured); %measured data
                    [model_y, model_x, model_LB, model_UB] = ecdf(dadc_predicted); %model predicted data
                    obs_ecdf = [obs_x,obs_y];
                    model_ecdf = [model_x, model_y]; 
                    obs_ecdf_table = array2table(obs_ecdf, VariableNames = ["X_obs", "Y_Obs"]);
                    model_ecdf_table = array2table(model_ecdf, VariableNames= ["X_model", "Y_model"]); 
                    ecdf_results.(scenarios(ss)).(ID).ecdf_observed = obs_ecdf_table;
                    ecdf_results.(scenarios(ss)).(ID).ecdf_model = model_ecdf_table; 
                    
                    %Pareto tail analysis (observed data)
                    LBo = dadc_threshold_measured/100; 
                    UBo = 1-(dadc_threshold_measured/100); 
                    pdo = paretotails(dadc_measured, LBo, UBo, "ecdf");
                    [pobs, qobs] = boundary(pdo); 
                    qobs_low = min(qobs);%observed threshold at selected prob cutoff
                    pareto_results.(scenarios(ss)).(ID).obs_distribution = pdo; 
                    pareto_results.(scenarios(ss)).(ID).obs_boundaries = [pobs,qobs]; 
                    pareto_results.(scenarios(ss)).(ID).obs_LB_cutoff = dadc_threshold_measured; 
                    pareto_results.(scenarios(ss)).(ID).obs_LB_dADCThreshold = qobs_low; 
            
                    % Pareto tail analysis (model data)
                    LB = dadc_threshold_model/100; 
                    UB = 1-(dadc_threshold_model/100); 
                    pd = paretotails(dadc_predicted, LB, UB, 'ecdf'); 
                    [p,q] = boundary(pd); 
                    q_low = min(q); %model threshold at selected prob cutoff
                    pareto_results.(scenarios(ss)).(ID).model_distribution = pd; 
                    pareto_results.(scenarios(ss)).(ID).model_boundaries = [p,q]; 
                    pareto_results.(scenarios(ss)).(ID).model_LB_cutoff = dadc_threshold_model; 
                    pareto_results.(scenarios(ss)).(ID).model_LB_dADCThreshold = q_low; 
            
                    %Voxel Classification by dADC patterns
                    [rows, cols] = size(subj_data); 
                    subj_data_classes = zeros(rows, cols + 1); 
                    pred_class = zeros(rows, 1); 
                    obs_class = pred_class; 
        
                    %Classify observed dADC as IDR/DMDR
                    for o = 1:length(dadc_measured)
                        targ_val = dadc_measured(o,1); 
                        if targ_val <=  qobs_low
                            obs_class(o,1) = 1; % IDR region (areas of increased diff. restriction)
                        else 
                            obs_class(o,1) = 2; % DMDR region (areas of decreased/maintained diff. restriction)
                        end
                    end
                
                    %Classify model predicted dADC as IDR/DMDR 
                    for m = 1:length(dadc_predicted)
                        targ_val = dadc_predicted(m,1); 
                        if targ_val <= q_low
                            pred_class(m,1) = 1; % IDR region (areas of increased diff. restriction)
                        else 
                            pred_class(m,1) = 2; % DMDR region (areas of decreased/maintained diff. restriction)
                        end
                    end
            
                    
                    %Classification Confusion Matrix
                    CM = confusionmat(obs_class, pred_class); 
                    fCM = reshape(CM', 1, numel(CM));
                    TP = fCM(1); 
                    FN = fCM(2); 
                    FP = fCM(3); 
                    TN = fCM(4);
            
                    %Calculate metrics
                    Sensitivity = TP/(TP+FN); 
                    Specificity = TN/(FP+TN); 
                    FPR = FP/(FP+TN);
            
                   %report classification results 
                   classification_results.(scenarios(ss)).(ID).pred_class = pred_class; 
                   classification_results.(scenarios(ss)).(ID).obs_class = obs_class; 
                   classification_results.(scenarios(ss)).(ID).confusion_matrix_array = CM; 
            
                   %Report metrics
                   if cc ==  length(scene_cases)
                       metrics_array(cc,:) = [id_num, TP, FP, TN, FN, Sensitivity, Specificity, FPR]; 
                       metrics_table = array2table(metrics_array, VariableNames=["CaseNum", "TP", "FP", "TN", "FN", "Sensitivity", "Specificity", "FPR"]); 
                       classification_results.(scenarios(ss)).MetricSummary = metrics_table; 
                       metrics_only_array = metrics_array(:,[6:8]); 
                       global_median = median(metrics_only_array); 
                       global_min = min(metrics_only_array); 
                       global_max = max(metrics_only_array);
                       global_mean = mean(metrics_only_array); 
                       global_sd = std(metrics_only_array); 

                       %save info on aggregated ROC variables
                       roc_aggregated_vars(ctr,1) = ss;
                       roc_aggregated_vars(ctr,2) = ocut_i; 
                       roc_aggregated_vars(ctr,3) = mcut_i;
                       roc_aggregated_vars(ctr,4) = global_mean(1,1); %sens
                       roc_aggregated_vars(ctr,5) = global_mean(1,2); %spec
                       roc_aggregated_vars(ctr,6) = global_mean(1,3); %FPR

                       perf_summary = [global_median; global_min; global_max];
                       roc_vars = [global_mean; global_sd]; 
                       glob_table = array2table(perf_summary, RowNames = ["Population Median", "Population Min", "Population Max"], ...
                       VariableNames=["Sensitivity", "Specificity", "FPR"]);
                       rocvars_table = array2table(roc_vars, RowNames=["Mean", "SD"], VariableNames=["Sensitivity", "Specificity", "FPR"]); 
                  
                       %save population level performance metrics
                       population_performance.(scenarios(ss)).PerformanceMetrics = glob_table;
                       population_performance.(scenarios(ss)).ROCMetrics = rocvars_table;

                   else
                       metrics_array(cc,:) = [id_num, TP, FP, TN, FN, Sensitivity, Specificity, FPR]; 

                   end
                
            
                end
                
                
                %Report Run Time 
                tf_execute = datetime("now"); 
                run_time = tf_execute - ti_execute; 
                run_seconds = seconds(run_time); 
                runtimes(ctr,1) = ocut_i; 
                runtimes(ctr,2) = mcut_i; 
                runtimes(ctr,3) = run_seconds; 

                all_results.(ocut_tag).(mcut_tag).classification_results = classification_results; 
                all_results.(ocut_tag).(mcut_tag).population_performance = population_performance;  
                all_results.(ocut_tag).(mcut_tag).ecdf_results = ecdf_results; 
                all_results.(ocut_tag).(mcut_tag).pareto_results = pareto_results;
                all_results.(ocut_tag).(mcut_tag).subjectlevel_performance = metrics_table;

            end
            
        end 

end

%% Save Data
runtimes(:,4) = runtimes(:,3)./ 60; 
disp(["Avg Run Time per cutoff pair (mins): ", num2str(mean(runtimes(:,4)))])
disp(["Avg Run Time per cutoff pair (sec): ", num2str(mean(runtimes(:,3)))])
runtime_table = array2table(runtimes, VariableNames=["Obs Cutoff", "Model Cutoff", "Execution Time (sec)", "Execution Time (mins)"]);
roc_agg_table = array2table(roc_aggregated_vars, VariableNames=["Scenario #", "Obs Cutff", "Model Cutoff", "Mean Sensitivity", "Mean Specificity", "Mean FPR"]); 

save("ExecutionTime.mat", "runtime_table", "-v7.3")
save("ROC_Data.mat", "roc_aggregated_vars", "roc_agg_table", "-v7.3")
save("ReproducedResults_AllCutoff.mat", "all_results", "-v7.3")


%% ROC Analysis
color_palette = ["#333E48", "#00393F", "#00535B", "#006D77", "#00A9B7"]; 
marker = ["-","--", ":", "-.", "-" ]; 
auc_array = zeros(length(thresholds)*4, 3); 
roc_fig = figure; 
tiledlayout(roc_fig,"flow")
ctr = 0; 

for scene = 1:4
    nexttile;
    fig_title = strcat("Scenario ", num2str(scene)); 

    %Set Scenario legend
    if scene == 1
        targ_region = "All regions"; 
        targ_scan = "RT_{end}"; 
    elseif scene == 2
        targ_region = "All regions";
        targ_scan = "RT_{1-month}"; 
    elseif scene == 3
        targ_region = "Enh regions"; 
        targ_scan = "RT_{end}";
    elseif scene == 4
        targ_region = "Enh regions"; 
        targ_scan = "RT_{1-month}";
    end
    
    %Determine Observed Data Cutoff - AUC metric
    for ooo = 1:length(thresholds)
        ctr = ctr + 1; 
        obs_metrics = roc_aggregated_vars(roc_aggregated_vars(:,1) == scene, :); 
        obs_metrics = obs_metrics(obs_metrics(:,2) == thresholds(ooo), :); 
        x = [0; obs_metrics(:,6); 1]; %fpr
        y = [0; obs_metrics(:,4); 1]; %sensitivity
        auc = trapz(x,y);
        auc_array(ctr,1) = scene; 
        auc_array(ctr,2) = thresholds(ooo); 
        auc_array(ctr,3) = auc;

        %plot
        figure(roc_fig)
        plot(x, y, "Color", color_palette(ooo), Marker = "o", MarkerFaceColor = color_palette(ooo), MarkerEdgeColor = color_palette(ooo), LineStyle = marker(ooo),LineWidth=1.25)
        hold on
    
        %plot settings
        ax = gca; 
        ax.XAxis.FontSize = 16; 
        ax.YAxis.FontSize = 16;
        ax.FontName = "Times New Roman"; 
        xlim([0,1.02])
        ylim([0,1.02])
        yticks([0:0.1:1])
        xticks([0:0.1:1])
        xlabel("False Positive Rate", FontSize=20)
        ylabel("True Positive Rate", FontSize=20)
        title(fig_title, FontSize = 24)
        subtitle (strcat("Region: ", targ_region, ", Timepoint: ", targ_scan), FontSize=14)
        plot([0,1], [0,1], '--', "Color", "black")
        hold on
    
        %plot legend (modify)
        if ooo == length(thresholds)
            leg_names = ["Lower 20", "Lower 25", "Lower 30", "Lower 35", "Lower 40"]; 
            lgd = legend(leg_names, 'Location', 'southeastoutside'); 
            title(lgd, "Observed Cutoff")
            fontsize(lgd, 16, "points")
        end
    end
end

figure(roc_fig)
sgtitle("ROC Analysis", FontSize = 32, FontName="Times New Roman")
%ssavefig(roc_fig, "ROC_Plot.fig")

%% Cutoff Selection 
auc_table = array2table(auc_array, VariableNames=["Scenario", "Model Cutoff", "AUC"]); 
best_auc = zeros(4,3); 
best_yi = zeros(4,3); 

%Select Observed Cutoff for all scenarios
for sc = 1:4
    sc_auc = auc_array(auc_array(:,1) == sc, :);
    [max_val, max_idx] = max(sc_auc(:,3)); 
    best_auc(sc,1) = sc; 
    best_auc(sc,2) = max_val; %Best AUC
    best_auc(sc,3) = auc_array(max_idx, 2); %LTCP Cutoff associated with best AUC

end

% Select Model Cutoff for all scenarios
for auc_i = 1:4
    LTCP_obs = best_auc(auc_i, 3); 
    best_curve_sensitivity = roc_aggregated_vars(roc_aggregated_vars(:,1) == auc_i, :);
    best_curve_sensitivity = best_curve_sensitivity(best_curve_sensitivity(:,2) == LTCP_obs,4); 
    best_curve_specificity = roc_aggregated_vars(roc_aggregated_vars(:,1) == auc_i, :); 
    best_curve_specificity = best_curve_specificity(best_curve_specificity(:,2) == LTCP_obs,5); 
    j_stats = round((best_curve_sensitivity + best_curve_specificity - 1),3);
    mod_cutoffs = roc_aggregated_vars(roc_aggregated_vars(:,1) == auc_i, :);
    mod_cutoffs = mod_cutoffs(mod_cutoffs(:,2) == LTCP_obs, 3); 
    yi_data = [mod_cutoffs, j_stats]; 

    %Get best YI
    [max_yi, max_yi_idx] = max(yi_data(:,2));
    yi_cutoff = yi_data(yi_data(:,2)== max_yi,1);
    sens_array_tie = best_curve_sensitivity(yi_data(:,2)== max_yi,1); 

    if length(yi_cutoff) > 1   
       [tie_val, tie_idx] = max(sens_array_tie(:,1)); 
       max_yi = yi_data(best_curve_sensitivity==tie_val, 2); 
       max_yi_idx = find(best_curve_sensitivity==tie_val); 
       max_yi_cutoff = yi_data(max_yi_idx,1); 
       
    elseif length(yi_cutoff) == 1
       max_yi_cutoff = yi_cutoff; 
    end

    best_yi(auc_i,1) = auc_i; %scenario
    best_yi(auc_i,2) = max_yi; %best YI
    best_yi(auc_i,3) = max_yi_cutoff; %model LTCP cutoff
end

%% Report best cutoffs per scenario (Table 1)
scene_description = ["All Regions, RT end", "All Regions, 1month", "Enh Regions, RT end", "Enh Regions, 1month"]; 
best_combination = [best_auc(:,1), scene_description', best_auc(:,3), best_yi(:,3), best_auc(:,2), best_yi(:,2)]; 
best_combination_table = array2table(best_combination, VariableNames=["Scenario", "Description", "LTCP Observed", "LTCP Model", ...
    "Best AUC", "Best Youden's Index"]); 

