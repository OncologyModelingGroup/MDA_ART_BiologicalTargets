%% Housekeeping
clear all; clc; close all; warning('off','all')
ti_execute = datetime("now"); 
ADC_data = load("ADC_Data.mat").ADC_Data;
subjects = string(fieldnames(ADC_data.Scenario_1));
allsubjects = double(strrep(subjects, "Case_", ""));

%% Get Scenarios
scenarios = string(fieldnames(ADC_data));
scenarios = scenarios(contains(scenarios, "Scenario")); 

%Compute the dADC for each scenario using optimal cutoffs
classification_results = struct(); 
pareto_results = struct(); 
ecdf_results = struct(); 
population_performance = struct(); 
subjectlevel_performance = struct(); 

for ss = 1:length(scenarios)
    
    %Get optimal lower tail cumulative probability cutoff determined in ROC
    scene = strrep(scenarios(ss), "_", ""); 
    dadc_threshold_measured = ADC_data.Descriptions.(scene).BestLTCP_Cutoff_Observed; %can be changed to [20, 25, 30, 35, 40]
    dadc_threshold_model = ADC_data.Descriptions.(scene).BestLTCP_Cutoff_Model; %can be changed to [20, 25, 30, 35, 40]

    %Get ADC data for each patient in all scenarios
    scene_cases = string(fieldnames(ADC_data.(scenarios(ss)))); 
    metrics_array = zeros(length(scene_cases), 9);

    %Prep figures
    ecdf_fig = figure;
    tiledlayout(ecdf_fig,"flow")
    cm_fig = figure; 
    tiledlayout(cm_fig, "flow")

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
    
        %Plot Statistical Mapping Strategy
        figure(ecdf_fig)
        nexttile; 
        plot(obs_x, obs_y, '-', LineWidth=2.5, Color="#333E48")
        hold on
        plot(model_x, model_y, '-', LineWidth=2.5, Color = "#BF5700")
        hold on
        xline(q_low, "--")
        hold on
        legend(["Measured", "Model Predicted", "Threshold"], Location="northeastoutside")
        xlabel("\DeltaADC")
        ylabel("eCDF")
        title(strrep(ID, "_", " "))
        fontsize(16,"points")
        fontname(ecdf_fig,"Times New Roman")
        hold off

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
        figure(cm_fig)
        nexttile; 
        cm = confusionchart(CM); 
        cm.Title = strrep(scene_cases(cc), "_", " "); 
        cm.FontSize = 16; 
        cm.FontName = 'Times New Roman'; 
        cm.DiagonalColor = "#00A9B7"; 
        cm.OffDiagonalColor = "#BF5700"; 
        fCM = reshape(CM', 1, numel(CM));
        TP = fCM(1); 
        FN = fCM(2); 
        FP = fCM(3); 
        TN = fCM(4);

        %Calculate metrics
        Sensitivity = TP/(TP+FN); 
        Specificity = TN/(FP+TN); 
        FPR = FP/(FP+TN);
        DSC  = calc_DICE(obs_class, pred_class); 

        % NOTE: the DSC line of code was modified to use linearized ADC data.
        % ADC values in the observed and model predicted data are
        % spatially mapped. Reported Dice-Sorensen coefficients presented in the publication
        % were calculated using 3D arrays of tumor volume.

       %report classification results 
       classification_results.(scenarios(ss)).(ID).pred_class = pred_class; 
       classification_results.(scenarios(ss)).(ID).obs_class = obs_class; 
       classification_results.(scenarios(ss)).(ID).confusion_matrix_array = CM; 


       %Report metrics
       if cc ==  length(scene_cases)
            metrics_array(cc,:) = [id_num, TP, FP, TN, FN, Sensitivity, Specificity, FPR, DSC]; 
            metrics_table = array2table(metrics_array, VariableNames=["CaseNum", "TP", "FP", "TN", "FN", "Sensitivity", "Specificity", "FPR", "DSC"]); 
           classification_results.(scenarios(ss)).MetricSummary = metrics_table; 
       else
           metrics_array(cc,:) = [id_num, TP, FP, TN, FN, Sensitivity, Specificity, FPR, DSC]; 
       end
    

    end
    
    %Calculate population level metrics
    if cc == length(scene_cases)
       metrics_only_array = metrics_array(:,[6,7,9]); 
       global_median = median(metrics_only_array); 
       global_min = min(metrics_only_array); 
       global_max = max(metrics_only_array); 
       perf_summary = [global_median; global_min; global_max]; 
       glob_table = array2table(perf_summary, RowNames = ["Population Median", "Population Min", "Population Max"], ...
           VariableNames=["Sensitivity", "Specificity", "DSC"]);
      
       %save population level performance metrics
       population_performance.(scenarios(ss)).PerformanceMetrics = glob_table; 

    end

%Save figures (pooled subjects)iif
figure(cm_fig)
sgtitle(strrep(scenarios(ss), "_", " "))
fontsize(16, "points")
savefig(cm_fig, strcat("ConfusionMatrices_", scenarios(ss),".fig"))
close(cm_fig)

figure(ecdf_fig)
sgtitle(strrep(scenarios(ss), "_", " "))
fontsize(16, "points")
savefig(ecdf_fig, strcat("ECDF_ParetoTail_", scenarios(ss),".fig"))
close(ecdf_fig)
end 

%Report Run Time 
tf_execute = datetime("now"); 
run_time = tf_execute - ti_execute; 
disp("Total Run Time: ")
disp(run_time)

%% Summary Metrics
sensitivity_array = nan(28,4); 
specificity_array = nan(28,4); 
dsc_array = nan(28,4); 

for ss = 1:length(scenarios)
    si = scenarios(ss); 
    s_metrics = table2array(classification_results.(si).MetricSummary);
    cases = s_metrics(:,1); 

    for cc = 1:length(cases)
        case_i = cases(cc); 
        sens = s_metrics(cc,6); 
        spec = s_metrics(cc,7); 
        dsc = s_metrics(cc,9); 

        %compile arrays
        sensitivity_array(case_i,ss) = sens; 
        specificity_array(case_i, ss) = spec; 
        dsc_array(case_i, ss) = dsc;
        
    end 
end

%Prep table data
sensitivity_array = sensitivity_array(~all(isnan(sensitivity_array), 2), :); 
specificity_array = specificity_array(~all(isnan(specificity_array), 2), :); 
dsc_array = dsc_array(~all(isnan(dsc_array), 2), :); 

sensitivity_table = array2table(sensitivity_array, RowNames = subjects, VariableNames=["Scenario1", "Scenario2", "Scenario3", "Scenario4"]); 
specificity_table = array2table(specificity_array, RowNames = subjects, VariableNames=["Scenario1", "Scenario2", "Scenario3", "Scenario4"]); 
dsc_table = array2table(dsc_array, RowNames = subjects, VariableNames=["Scenario1", "Scenario2", "Scenario3", "Scenario4"]); 

subjectlevel_performance.Sensitivity = sensitivity_table; 
subjectlevel_performance.Specificity = specificity_table; 
subjectlevel_performance.DSC = dsc_table; 

%% Save Generated Data for All Scenarios
save("ReproducedResults_OptimalCutoff.mat", "classification_results", "population_performance", "subjectlevel_performance", "ecdf_results", "-v7.3")

%% functions

function [score] = calc_DICE(array1, array2)

%Get dimensions
[sz1_x, ~, ~] = size(array1); 
[sz2_x, ~, ~] = size(array2); 

%find common elements
for i = 1:sz1_x
    a_elem = array1(i); 
    b_elem = array2(i); 

    if a_elem == b_elem
        similar_elements(i) = 1; 
    else
        continue
    end
end

%Calculate score
num_common = nnz(similar_elements); 
score = 2 * num_common / (sz1_x + sz2_x); 

end
