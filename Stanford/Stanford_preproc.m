% Learning Gait Patterns
% nrokh 2021

clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   input: Stanford University dataset from SimTK (Meta and SubjectData)
%   output: toe-in heuristics to fit FDA, baseline KAM 
%   utils: find_FPA.m, find_PelvisCenter.m, interp100.m,
%          wrtPelvis.m, subjectAvgTrajectory.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EDIT as needed:
destFolder = "D:\sKAM Classifier (2021)\FDA\";
%% 0. prep structure for all subjects

load('subjectHeightsWeights.mat') % (sub_num, weight [kg], height [cm])
load('store_events_trial.mat')    % heel strike, toe off indices (toe-in)    
load('store_events_baseline.mat') % heel strike, toe off indices (baseline)

% a. label ID cells
for i=1:1:12
    if i<10
        subID_baseline{i} = ['sub0' num2str(i,'%d') '_baseline.mat'];
        subID_trial{i} = ['sub0' num2str(i,'%d') '_toein.mat'];
    else
        subID_baseline{i} = ['sub' num2str(i,'%d') '_baseline.mat'];
        subID_trial{i} = ['sub' num2str(i,'%d') '_toein.mat'];
    end
end

% b. initialize parameters
% (12 subjects, up to 10 steps, within those cells store 101x3 in x, y, z)
store_GRF_b = cell([12, 10]); 
store_GRF_t = cell([12, 10]);
store_pelvis_b = cell([12, 10]); 
store_pelvis_t = cell([12, 10]);
store_COP_b = cell([12, 10]); 
store_COP_t = cell([12, 10]);
store_KJC_b = cell([12, 10]); 
store_KJC_t = cell([12, 10]);
store_FPA_b = cell([12, 10]); 
store_FPA_t = cell([12, 10]);

% filter vals (15 Hz 4th order zero lag butterworth low pass)
GRFFilterFreqHz = 15;
GRFFilterOrder=2;
grfSampFreq = 1000; % all subjects except for s08 sampled at 1 kHz
cutoffRatio = 0.05;

%% 1. load data from all subjects and set up parameter stores


for i = 1:1:12
    if i ==8
        grfSampFreq = 2000; 
    else
        grfSampFreq = 1000;
    end
    % a. load all raw data for the subject using ID name string
    s_baseline = load(subID_baseline{i});
    s_trial = load(subID_trial{i});
    events_baseline = store_events_baseline{i};
    events_trial = store_events_trial{i};
    
    % b. account for different sampling rates in subject 8 and between markers vs GRF/COP
    % (events_short is for markers, events_long is for GRF/COP)
    if i == 8
        events_baseline_short = round(store_events_baseline{i}/20);
        events_trial_short = round(store_events_trial{i}/20);
    else
        events_baseline_short = round(store_events_baseline{i}/16);
        events_trial_short = round(store_events_trial{i}/16);
    end
    
    % c. find FPA across entire trial for the subject (unsegmented)
    % (subjects 3, 4, 6 and 12 are RT foot)
    if i == 3 || i==4 || i==6 || i==12
        s_FPA_baseline = find_FPA(s_baseline.baseline.markers, 'RT'); % row by row
        s_FPA_trial = find_FPA(s_trial.tibia.markers, 'RT');
    else
        s_FPA_baseline = find_FPA(s_baseline.baseline.markers, 'LT');
        s_FPA_trial = find_FPA(s_trial.tibia.markers, 'LT');
    end
    
    % d. find pelvis centers across entire trial for subject(unsegmented)
    [s_pelvis_baseline, s_pelvis_trial] = find_PelvisCenter(i, s_baseline.baseline.markers, s_trial.tibia.markers);
    
    % e. save GRF, COP, KJC, shank
    if i == 3 || i==4 || i==6 || i==12
        s_GRF_baseline = s_baseline.baseline.GRF.FP_RTleg;
        s_GRF_trial = s_trial.tibia.GRF.FP_RTleg;
        s_COP_baseline = s_baseline.baseline.COP.FP_RTleg;
        s_COP_trial = s_trial.tibia.COP.FP_RTleg;
        
        % add thresholding
        thresholdValue = subjectHeightsWeights(i,2)*9.81*cutoffRatio;
        contactIndexb = s_GRF_baseline(3,:)>thresholdValue;
        contactIndext = s_GRF_trial(3,:)>thresholdValue;
        s_GRF_baseline(:,~contactIndexb) = 0;
        s_GRF_trial(:,~contactIndext) = 0;
        s_COP_baseline(:,~contactIndexb) = 0;
        s_COP_trial(:,~contactIndext) = 0;
        
        % add 15Hz LP filter (06.24.21)
        ratio = GRFFilterFreqHz/(grfSampFreq/2);
        [b,a] = butter(GRFFilterOrder, ratio);
        s_GRF_baseline = filtfilt(b, a, s_GRF_baseline')';
        s_GRF_trial = filtfilt(b, a, s_GRF_trial')';
        s_COP_baseline = filtfilt(b, a, s_COP_baseline')';
        s_COP_trial = filtfilt(b, a, s_COP_trial')';
        
    else
        s_GRF_baseline = s_baseline.baseline.GRF.FP_LTleg;
        s_GRF_trial = s_trial.tibia.GRF.FP_LTleg;
        s_COP_baseline = s_baseline.baseline.COP.FP_LTleg;
        s_COP_trial = s_trial.tibia.COP.FP_LTleg;
        
        % add thresholding
        thresholdValue = subjectHeightsWeights(i,2)*9.81*cutoffRatio;
        contactIndexb = s_GRF_baseline(3,:)>thresholdValue;
        contactIndext = s_GRF_trial(3,:)>thresholdValue;
        s_GRF_baseline(:,~contactIndexb) = 0;
        s_GRF_trial(:,~contactIndext) = 0;
        s_COP_baseline(:,~contactIndexb) = 0;
        s_COP_trial(:,~contactIndext) = 0;
        
        % add 15Hz LP filter (06.24.21)
        ratio = GRFFilterFreqHz/(grfSampFreq/2);
        [b,a] = butter(GRFFilterOrder, ratio);
        s_GRF_baseline = filtfilt(b, a, s_GRF_baseline')';
        s_GRF_trial = filtfilt(b, a, s_GRF_trial')';
        s_COP_baseline = filtfilt(b, a, s_COP_baseline')';
        s_COP_trial = filtfilt(b, a, s_COP_trial')';
    end 
    s_KJC_baseline = s_baseline.baseline.markers.KJC;
    s_KJC_trial = s_trial.tibia.markers.KJC;


     % f. loop through every step and save to cell
    step_iter = 1;
    for j = 1:1:10
        % FPA
        store_FPA_b{i, j} = mean(s_FPA_baseline(events_baseline_short(step_iter):events_baseline_short(step_iter+1)));
        store_FPA_t{i, j} = mean(s_FPA_trial(events_trial_short(step_iter):events_trial_short(step_iter+1)));
        
        % GRF
        store_GRF_b{i,j} = s_GRF_baseline(:, events_baseline(step_iter):events_baseline(step_iter+1));
        store_GRF_t{i,j} = s_GRF_trial(:, events_trial(step_iter):events_trial(step_iter+1));
        
          % change ML dir to be opposite for R leg
        if i == 3||i==4||i==6||i==12
            store_GRF_b{i,j}(2,:) = -store_GRF_b{i,j}(2,:);
            store_GRF_t{i,j}(2,:) = -store_GRF_t{i,j}(2,:);
        end
        
        
        % COP
        store_COP_b{i,j} = s_COP_baseline(:, events_baseline(step_iter):events_baseline(step_iter+1));
        store_COP_t{i,j} = s_COP_trial(:, events_trial(step_iter):events_trial(step_iter+1));
        
        % KJC
        store_KJC_b{i,j} = s_KJC_baseline(:, events_baseline_short(step_iter):events_baseline_short(step_iter+1));
        store_KJC_t{i,j} = s_KJC_trial(:, events_trial_short(step_iter):events_trial_short(step_iter+1));
        
        % pelvis
        store_pelvis_b{i,j} = s_pelvis_baseline{1,i}(:, events_baseline_short(step_iter):events_baseline_short(step_iter+1));
        store_pelvis_t{i,j} = s_pelvis_trial{1,i}(:, events_trial_short(step_iter):events_trial_short(step_iter+1));
        
        step_iter = step_iter+2;
    end
    
end

% g. interpolate to 100 (% stance):
store_GRF_b = interp100(store_GRF_b);
store_GRF_t = interp100(store_GRF_t);
store_COP_b = interp100(store_COP_b);
store_COP_t = interp100(store_COP_t);
store_KJC_b = interp100(store_KJC_b);
store_KJC_t = interp100(store_KJC_t);
store_pelvis_b = interp100(store_pelvis_b);
store_pelvis_t = interp100(store_pelvis_t);

% h. relate COP, KJC to pelvis (in x,y,z):
store_KJC_b_p = wrtPelvis(store_KJC_b, store_pelvis_b);
store_KJC_t_p = wrtPelvis(store_KJC_t, store_pelvis_t);
store_COP_b_p = wrtPelvis(store_COP_b, store_pelvis_b);
store_COP_t_p = wrtPelvis(store_COP_t, store_pelvis_t);

% i. compute average trajectories at baseline and trial for future validation:
[baseline_pelvis_x, baseline_pelvis_y, baseline_pelvis_z] = subjectAvgTrajectory(store_pelvis_t);
[baseline_GRF_x, baseline_GRF_y, baseline_GRF_z] = subjectAvgTrajectory(store_GRF_b);
[baseline_COP_x, baseline_COP_y, baseline_COP_z] = subjectAvgTrajectory(store_COP_b_p);
[baseline_KJC_x, baseline_KJC_y, baseline_KJC_z] = subjectAvgTrajectory(store_KJC_b_p);

[trial_pelvis_x, trial_pelvis_y, trial_pelvis_z] = subjectAvgTrajectory(store_pelvis_t);
[trial_GRF_x, trial_GRF_y, trial_GRF_z] = subjectAvgTrajectory(store_GRF_t);
[trial_COP_x, trial_COP_y, trial_COP_z] = subjectAvgTrajectory(store_COP_t_p);
[trial_KJC_x, trial_KJC_y, trial_KJC_z] = subjectAvgTrajectory(store_KJC_t_p);

% j. compute target FPA

allfpab = cell2mat(store_FPA_b); 
allfpat = cell2mat(store_FPA_t);
fpabmean = mean(allfpab,2);         % find subject's average baseline FPA
fpat = allfpat - fpabmean;          % relate subject's toe-in FPAs to baseline
target_FPA = round(mean(allfpab,2) - mean(allfpat,2));


%% 2. prepare data and export to Python for binning and FDA

% a. compute average trajectories wrt pelvis:

[baseline_KJC_x_p, baseline_KJC_y_p, baseline_KJC_z_p] = subjectAvgTrajectory(store_KJC_b_p);
[baseline_COP_x_p, baseline_COP_y_p, baseline_COP_z_p] = subjectAvgTrajectory(store_COP_b_p);

[trial_COP_x_p, trial_COP_y_p, trial_COP_z_p] = subjectAvgTrajectory(store_COP_t_p);
[trial_KJC_x_p, trial_KJC_y_p, trial_KJC_z_p] = subjectAvgTrajectory(store_KJC_t_p);

% b. append toe-in FPA 

for i = 1:1:12
    for j = 1:1:10
        store_KJC_t_p{i,j} = [ round(fpat(i,j)), store_KJC_t_p{i,j}(1,:); round(fpat(i,j)), store_KJC_t_p{i,j}(2,:); round(fpat(i,j)), store_KJC_t_p{i,j}(3,:)];
        store_COP_t_p{i,j} = [ round(fpat(i,j)), store_COP_t_p{i,j}(1,:); round(fpat(i,j)), store_COP_t_p{i,j}(2,:); round(fpat(i,j)), store_COP_t_p{i,j}(3,:)];
    end
end

store_KJC_t_p_diff = store_KJC_t_p;
store_COP_t_p_diff = store_COP_t_p;

% c. find difference between toe-in and baseline trajectories
for ii = 1:1:12
    for j = 1:1:10

        % subtract trajectory (AP and ML)
        store_KJC_t_p_diff{ii,j}(1, 2:end) = store_KJC_t_p{ii,j}(1, 2:end) - baseline_KJC_x_p(ii, :);
        store_KJC_t_p_diff{ii,j}(2, 2:end) = store_KJC_t_p{ii,j}(2, 2:end) - baseline_KJC_y_p(ii, :);
        store_KJC_t_p_diff{ii,j}(3, 2:end) = store_KJC_t_p{ii,j}(3, 2:end);
        store_COP_t_p_diff{ii,j}(1, 2:end) = store_COP_t_p{ii,j}(1, 2:end) - baseline_COP_x_p(ii, :);
        store_COP_t_p_diff{ii,j}(2, 2:end) = store_COP_t_p{ii,j}(2, 2:end) - baseline_COP_y_p(ii, :);
        store_COP_t_p_diff{ii,j}(3, 2:end) = store_COP_t_p{ii,j}(3, 2:end);
        
    end
end


% d. select out x,y params:
all_KJC_t_x = cellfun(@(v)v(1,:), store_KJC_t_p_diff, 'UniformOutput', false);
all_KJC_t_y = cellfun(@(v)v(2,:), store_KJC_t_p_diff, 'UniformOutput', false);

all_COP_t_x = cellfun(@(v)v(1,:), store_COP_t_p_diff, 'UniformOutput', false);
all_COP_t_y = cellfun(@(v)v(2,:), store_COP_t_p_diff, 'UniformOutput', false);


% e. save to destination folder
save(strcat(destFolder, "KJC_x.mat"), 'all_KJC_t_x');
save(strcat(destFolder, "KJC_y.mat"), 'all_KJC_t_y');
save(strcat(destFolder, "COP_x.mat"), 'all_COP_t_x');
save(strcat(destFolder, "COP_y.mat"), 'all_COP_t_y');

% now run Python script, to generate LOOCV for 12 subjects or overall fits

