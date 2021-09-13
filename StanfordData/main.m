% Synthetic KAM generation
% nrokh 08.05.2020
clear all; clc;

%% 1. load data for all subjects
load('subjectHeightsWeights.mat')
% a. label ID cells

load('store_events_trial.mat')          % load foot strike indices
load('store_events_baseline.mat')
for i=1:1:12
    if i<10
        subID_baseline{i} = ['sub0' num2str(i,'%d') '_baseline.mat'];
        subID_trial{i} = ['sub0' num2str(i,'%d') '_toein.mat'];
    else
        subID_baseline{i} = ['sub' num2str(i,'%d') '_baseline.mat'];
        subID_trial{i} = ['sub' num2str(i,'%d') '_toein.mat'];
    end
end

% b. initialize storage cells
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
store_shankCM_b = cell([12, 10]); 
store_shankCM_t = cell([12, 10]);
store_KAM_baseline_debug = cell([12,10]);
store_KAM_trial_debug = cell([12,10]);
store_KAM_baseline_debug2 = cell([12,10]);
store_KAM_trial_debug2 = cell([12,10]);
store_KAM_baseline_debug3 = cell([12,10]);
store_KAM_trial_debug3 = cell([12,10]);

    peaks_KAM = cell([12,1]);
    peaks_tKAM = cell([12,1]);
    
    GRFFilterFreqHz = 15;
    GRFFilterOrder=2;
    grfSampFreq = 1000;
    cutoffRatio = 0.05;

for i = 1:1:12
    if i ==8
        grfSampFreq = 2000; % idk if this is the right 
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
    s_shank_baseline = s_baseline.baseline.markers.shank_cM;
    s_shank_trial = s_trial.tibia.markers.shank_cM;
    
    %%%%% DEBUG 08.27.20
    
    % save KAM
    sCOPbx = interp1((linspace(1, length(s_COP_baseline), length(s_COP_baseline))), s_COP_baseline(1,:), linspace(1, length(s_COP_baseline), length(s_KJC_baseline)));
sCOPby = interp1((linspace(1, length(s_COP_baseline), length(s_COP_baseline))), s_COP_baseline(2,:), linspace(1, length(s_COP_baseline), length(s_KJC_baseline)));
sCOPbz = interp1((linspace(1, length(s_COP_baseline), length(s_COP_baseline))), s_COP_baseline(3,:), linspace(1, length(s_COP_baseline), length(s_KJC_baseline)));
sCOPb = [sCOPbx/1000; sCOPby/1000; sCOPbz/1000];

sCOPtx = interp1((linspace(1, length(s_COP_trial), length(s_COP_trial))), s_COP_trial(1,:), linspace(1, length(s_COP_trial), length(s_KJC_trial)));
sCOPty = interp1((linspace(1, length(s_COP_trial), length(s_COP_trial))), s_COP_trial(2,:), linspace(1, length(s_COP_trial), length(s_KJC_trial)));
sCOPtz = interp1((linspace(1, length(s_COP_trial), length(s_COP_trial))), s_COP_trial(3,:), linspace(1, length(s_COP_trial), length(s_KJC_trial)));
sCOPt = [sCOPtx/1000; sCOPty/1000; sCOPtz/1000];

sGRFbx = interp1((linspace(1, length(s_GRF_baseline), length(s_GRF_baseline))), s_GRF_baseline(1,:), linspace(1, length(s_GRF_baseline), length(s_KJC_baseline)));
sGRFby = interp1((linspace(1, length(s_GRF_baseline), length(s_GRF_baseline))), s_GRF_baseline(2,:), linspace(1, length(s_GRF_baseline), length(s_KJC_baseline)));
sGRFbz = interp1((linspace(1, length(s_GRF_baseline), length(s_GRF_baseline))), s_GRF_baseline(3,:), linspace(1, length(s_GRF_baseline), length(s_KJC_baseline)));
sGRFb = [sGRFbx; sGRFby; sGRFbz];

sGRFtx = interp1((linspace(1, length(s_GRF_trial), length(s_GRF_trial))), s_GRF_trial(1,:), linspace(1, length(s_GRF_trial), length(s_KJC_trial)));
sGRFty = interp1((linspace(1, length(s_GRF_trial), length(s_GRF_trial))), s_GRF_trial(2,:), linspace(1, length(s_GRF_trial), length(s_KJC_trial)));
sGRFtz = interp1((linspace(1, length(s_GRF_trial), length(s_GRF_trial))), s_GRF_trial(3,:), linspace(1, length(s_GRF_trial), length(s_KJC_trial)));
sGRFt = [sGRFtx; sGRFty; sGRFtz];

r = sCOPb  - s_KJC_baseline/1000;
F = sGRFb;
KM1 = cross(r,F);
KAM_baseline_debug = KM1(1,:)/((subjectHeightsWeights(i,2)*subjectHeightsWeights(i,3)/100))*9.81;
KAM_baseline_debug2 = s_baseline.baseline.kinetics.knee_adduction_moment  ;

KFM_baseline_debug = KM1(2,:)/((subjectHeightsWeights(i,2)*subjectHeightsWeights(i,3)/100))*9.81;
KFM_baseline_debug2 = s_baseline.baseline.kinetics.knee_flexion_moment  ;

r = sCOPt - s_KJC_trial/1000;
F = sGRFt;
KM2 = cross(r,F);
KAM_trial_debug = KM2(1,:)/((subjectHeightsWeights(i,2)*subjectHeightsWeights(i,3)/100))*9.81;
KAM_trial_debug2 = s_trial.tibia.kinetics.knee_adduction_moment  ;




% NR 05.27.2021 - building T_LAB_ACS:
KAMb_ACS = [];
KAMb_ACS2 = [];
for k = 1:1:length(s_baseline.baseline.markers.lfem)
    % to check: if you use the lateral marker instead of medial does it work
    % also: do you need to change coordinate system for L vs R
  
    yacs = s_baseline.baseline.markers.lfem(:,k) - s_baseline.baseline.markers.KJC(:,k);
    yacs = yacs/norm(yacs);
    xacs = cross( s_baseline.baseline.markers.KJC(:,k) -s_baseline.baseline.markers.AJC(:,k), yacs);
    xacs = xacs/norm(xacs);
    zacs = cross (xacs, yacs);
    
    bT_LAB_ACS = [[xacs;0],[yacs;0],[zacs;0],[s_baseline.baseline.markers.KJC(:,k);1]] ;
    
    % try multiplying KAM by this 
    KMb_ACS = bT_LAB_ACS(1:3,1:3)'*KM1(:,k);
    KAMb_ACS(k) = KMb_ACS(1)/((subjectHeightsWeights(i,2)*subjectHeightsWeights(i,3)/100))*9.81; 
    
    % NOTES AND TODO:
    %i think it's close, or at least better; incorporate this into trial
    %and also the computing KAM down below and see if it improves second pk
    
    % try multiplying COP, KJC, GRF by ACS and recomputing KAM
    r = bT_LAB_ACS(1:3,1:3)'*sCOPb(:,k)  - bT_LAB_ACS(1:3,1:3)'*(s_KJC_baseline(:,k)/1000);
    F = bT_LAB_ACS(1:3,1:3)'*sGRFb(:,k);
    KMb_ACS2 = cross(r,F);
    KAMb_ACS2(k) = KMb_ACS2(1)/((subjectHeightsWeights(i,2)*subjectHeightsWeights(i,3)/100))*9.81; 
    
    
    

    
    
    
    
end

figure;
hold on
plot(KAMb_ACS)
plot(KAMb_ACS2)


KAMt_ACS = [];
for k = 1:1:length(s_trial.tibia.markers.lfem)
    % to check: if you use the lateral marker instead of medial does it work
    % also: do you need to change coordinate system for L vs R
  
    yacs = s_trial.tibia.markers.lfem(:,k) - s_trial.tibia.markers.KJC(:,k);
    yacs = yacs/norm(yacs);
    xacs = cross( s_trial.tibia.markers.KJC(:,k) -s_trial.tibia.markers.AJC(:,k), yacs);
    xacs = xacs/norm(xacs);
    zacs = cross (xacs, yacs);
    
    bT_LAB_ACS = [[xacs;0],[yacs;0],[zacs;0],[s_trial.tibia.markers.KJC(:,k);1]] ;
    
    % try multiplying KAM by this 
    KMt_ACS = bT_LAB_ACS(1:3,1:3)'*KM2(:,k);
    KAMt_ACS(k) = KMt_ACS(1)/((subjectHeightsWeights(i,2)*subjectHeightsWeights(i,3)/100))*9.81; 
    
    % NOTES AND TODO:
    %i think it's close, or at least better; incorporate this into trial
    %and also the computing KAM down below and see if it improves second pk
    r = bT_LAB_ACS(1:3,1:3)'*sCOPt(:,k)  - bT_LAB_ACS(1:3,1:3)'*(s_KJC_trial(:,k)/1000);
    F = bT_LAB_ACS(1:3,1:3)'*sGRFt(:,k);
    KMt_ACS2 = cross(r,F);
    KAMt_ACS2(k) = KMt_ACS2(1)/((subjectHeightsWeights(i,2)*subjectHeightsWeights(i,3)/100))*9.81; 
    
    
end
figure;
hold on
plot(KAMt_ACS)
plot(KAMt_ACS2)





    
    % f. loop through every step and save to cell
    step_iter = 1;
    for j = 1:1:10
        % FPA
        store_FPA_b{i, j} = mean(s_FPA_baseline(events_baseline_short(step_iter):events_baseline_short(step_iter+1)));
        store_FPA_t{i, j} = mean(s_FPA_trial(events_trial_short(step_iter):events_trial_short(step_iter+1)));
        
        % GRF
        store_GRF_b{i,j} = s_GRF_baseline(:, events_baseline(step_iter):events_baseline(step_iter+1));
        store_GRF_t{i,j} = s_GRF_trial(:, events_trial(step_iter):events_trial(step_iter+1));
        
        store_KAM_baseline_debug{i,j} = KAM_baseline_debug(events_baseline_short(step_iter):events_baseline_short(step_iter+1));
        store_KAM_trial_debug{i,j} = KAM_trial_debug(events_trial_short(step_iter):events_trial_short(step_iter+1));
        store_KAM_baseline_debug2{i,j} = KAM_baseline_debug2(events_baseline_short(step_iter):events_baseline_short(step_iter+1));
        store_KAM_trial_debug2{i,j} = KAM_trial_debug2(events_trial_short(step_iter):events_trial_short(step_iter+1));
        
        % ACS dbug
        store_KAM_baseline_debug3{i,j} = KAMb_ACS( events_baseline_short(step_iter):events_baseline_short(step_iter+1));
        store_KAM_trial_debug3{i,j} = KAMt_ACS( events_trial_short(step_iter):events_trial_short(step_iter+1));
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
        
        % shank
        store_shankCM_b{i,j} = s_shank_baseline(:, events_baseline_short(step_iter):events_baseline_short(step_iter+1));
        store_shankCM_t{i,j} = s_shank_trial(:, events_trial_short(step_iter):events_trial_short(step_iter+1));
        
        % pelvis
        store_pelvis_b{i,j} = s_pelvis_baseline{1,i}(:, events_baseline_short(step_iter):events_baseline_short(step_iter+1));
        store_pelvis_t{i,j} = s_pelvis_trial{1,i}(:, events_trial_short(step_iter):events_trial_short(step_iter+1));
        
        step_iter = step_iter+2;
    end
    
end

%% NR 06.02.2021 trying to convert to TCS step-by-step
% COP, KJC, pelvis

% i don't think this will work the way you want it to
% store_COP_b2 = cell(12,10);
% store_COP_t2 = cell(12,10);
% store_GRF_b2 = cell(12,10);
% store_GRF_t2 = cell(12,10);
% for i=1:1:12
%     s_baseline = load(subID_baseline{i});
%     
%     for j = 1:1:10
%         % downsample COP to length of KJC
%         store_COP_b2{i,j}(1,:) = interp1(linspace(1, length(store_COP_b{i,j}(1,:)), length(store_COP_b{i,j}(1,:))), store_COP_b{i,j}(1,:), linspace(1, length(store_COP_b{i,j}(1,:)), length(store_KJC_b{i,j}(1,:))));
%         store_COP_b2{i,j}(2,:) = interp1(linspace(1, length(store_COP_b{i,j}(2,:)), length(store_COP_b{i,j}(2,:))), store_COP_b{i,j}(2,:), linspace(1, length(store_COP_b{i,j}(2,:)), length(store_KJC_b{i,j}(2,:))));
%         store_COP_b2{i,j}(3,:) = interp1(linspace(1, length(store_COP_b{i,j}(3,:)), length(store_COP_b{i,j}(3,:))), store_COP_b{i,j}(3,:), linspace(1, length(store_COP_b{i,j}(3,:)), length(store_KJC_b{i,j}(3,:))));
%         
%         store_GRF_b2{i,j}(1,:) = interp1(linspace(1, length(store_GRF_b{i,j}(1,:)), length(store_GRF_b{i,j}(1,:))), store_GRF_b{i,j}(1,:), linspace(1, length(store_GRF_b{i,j}(1,:)), length(store_KJC_b{i,j}(1,:))));
%         store_GRF_b2{i,j}(2,:) = interp1(linspace(1, length(store_GRF_b{i,j}(2,:)), length(store_GRF_b{i,j}(2,:))), store_GRF_b{i,j}(2,:), linspace(1, length(store_GRF_b{i,j}(2,:)), length(store_KJC_b{i,j}(2,:))));
%         store_GRF_b2{i,j}(3,:) = interp1(linspace(1, length(store_GRF_b{i,j}(3,:)), length(store_GRF_b{i,j}(3,:))), store_GRF_b{i,j}(3,:), linspace(1, length(store_GRF_b{i,j}(3,:)), length(store_KJC_b{i,j}(3,:))));
%        
%         % get events
%         if i == 8
%             eventsb = round(store_events_baseline{1,i}/20);
%         else
%             eventsb = round(store_events_baseline{1,i}/16);
%         end
%         
%         % create ACS at each point during that step
%         for k= eventsb(2*j-1):1:eventsb(2*j)
%             yacs = s_baseline.baseline.markers.lfem(:,k) - s_baseline.baseline.markers.KJC(:,k);
%             yacs = yacs/norm(yacs);
%             xacs = cross( s_baseline.baseline.markers.KJC(:,k) -s_baseline.baseline.markers.AJC(:,k), yacs);
%             xacs = xacs/norm(xacs);
%             zacs = cross (xacs, yacs);
%             
%             bT_LAB_ACS = [[xacs;0],[yacs;0],[zacs;0],[s_baseline.baseline.markers.KJC(:,k);1]] ;
%             
%             % multiply COP pelvis KJC 
%             store_COP_b2{i,j}(:,1+(k-eventsb(2*j-1))) = bT_LAB_ACS(1:3,1:3)'*store_COP_b2{i,j}(:,1+(k-eventsb(2*j-1)));
%             store_GRF_b2{i,j}(:,1+(k-eventsb(2*j-1))) = bT_LAB_ACS(1:3,1:3)'*store_GRF_b2{i,j}(:,1+(k-eventsb(2*j-1)));
%             store_KJC_b{i,j}(:,1+(k-eventsb(2*j-1))) = bT_LAB_ACS(1:3,1:3)'*store_KJC_b{i,j}(:,1+(k-eventsb(2*j-1)));
%             store_pelvis_b{i,j}(:,1+(k-eventsb(2*j-1))) = bT_LAB_ACS(1:3,1:3)'*store_pelvis_b{i,j}(:,1+(k-eventsb(2*j-1)));
%             
%             % make COPz all zeros again
%             store_COP_b2{i,j}(3,1+(k-eventsb(2*j-1))) = 0;
%         end
%         
%         
%     end
% end
%         
% for i=1:1:12
%     s_trial = load(subID_trial{i});
%     
%     for j = 1:1:10
%         % downsample COP to length of KJC
%         store_COP_t2{i,j}(1,:) = interp1(linspace(1, length(store_COP_t{i,j}(1,:)), length(store_COP_t{i,j}(1,:))), store_COP_t{i,j}(1,:), linspace(1, length(store_COP_t{i,j}(1,:)), length(store_KJC_t{i,j}(1,:))));
%         store_COP_t2{i,j}(2,:) = interp1(linspace(1, length(store_COP_t{i,j}(2,:)), length(store_COP_t{i,j}(2,:))), store_COP_t{i,j}(2,:), linspace(1, length(store_COP_t{i,j}(2,:)), length(store_KJC_t{i,j}(2,:))));
%         store_COP_t2{i,j}(3,:) = interp1(linspace(1, length(store_COP_t{i,j}(3,:)), length(store_COP_t{i,j}(3,:))), store_COP_t{i,j}(3,:), linspace(1, length(store_COP_t{i,j}(3,:)), length(store_KJC_t{i,j}(3,:))));
%         
%         store_GRF_t2{i,j}(1,:) = interp1(linspace(1, length(store_GRF_t{i,j}(1,:)), length(store_GRF_t{i,j}(1,:))), store_GRF_t{i,j}(1,:), linspace(1, length(store_GRF_t{i,j}(1,:)), length(store_KJC_t{i,j}(1,:))));
%         store_GRF_t2{i,j}(2,:) = interp1(linspace(1, length(store_GRF_t{i,j}(2,:)), length(store_GRF_t{i,j}(2,:))), store_GRF_t{i,j}(2,:), linspace(1, length(store_GRF_t{i,j}(2,:)), length(store_KJC_t{i,j}(2,:))));
%         store_GRF_t2{i,j}(3,:) = interp1(linspace(1, length(store_GRF_t{i,j}(3,:)), length(store_GRF_t{i,j}(3,:))), store_GRF_t{i,j}(3,:), linspace(1, length(store_GRF_t{i,j}(3,:)), length(store_KJC_t{i,j}(3,:))));
%        
%         % get events
%         if i == 8
%             eventst = round(store_events_trial{1,i}/20);
%         else
%             eventst = round(store_events_trial{1,i}/16);
%         end
%         
%         % create ACS at each point during that step
%         for k= eventst(2*j-1):1:eventst(2*j)
%             yacs = s_trial.tibia.markers.lfem(:,k) - s_trial.tibia.markers.KJC(:,k);
%             yacs = yacs/norm(yacs);
%             xacs = cross( s_trial.tibia.markers.KJC(:,k) -s_trial.tibia.markers.AJC(:,k), yacs);
%             xacs = xacs/norm(xacs);
%             zacs = cross (xacs, yacs);
%             
%             bT_LAB_ACS = [[xacs;0],[yacs;0],[zacs;0],[s_trial.tibia.markers.KJC(:,k);1]] ;
%             
%             % multiply GRF COP pelvis KJC 
%             store_COP_t2{i,j}(:,1+(k-eventst(2*j-1))) = bT_LAB_ACS(1:3,1:3)'*store_COP_t2{i,j}(:,1+(k-eventst(2*j-1)));
%             store_GRF_t2{i,j}(:,1+(k-eventst(2*j-1))) = bT_LAB_ACS(1:3,1:3)'*store_GRF_t2{i,j}(:,1+(k-eventst(2*j-1)));
%             store_KJC_t{i,j}(:,1+(k-eventst(2*j-1))) = bT_LAB_ACS(1:3,1:3)'*store_KJC_t{i,j}(:,1+(k-eventst(2*j-1)));
%             store_pelvis_t{i,j}(:,1+(k-eventst(2*j-1))) = bT_LAB_ACS(1:3,1:3)'*store_pelvis_t{i,j}(:,1+(k-eventst(2*j-1)));
%             
%             store_COP_t2{i,j}(3,1+(k-eventst(2*j-1))) = 0;
%         end
%         
%     end
% end
    
    
    

%%
% g. interpolate to 100:
store_GRF_b = interp100(store_GRF_b);
store_GRF_t = interp100(store_GRF_t);
store_COP_b = interp100(store_COP_b);
store_COP_t = interp100(store_COP_t);
% store_GRF_b = interp100(store_GRF_b2);
% store_GRF_t = interp100(store_GRF_t2);
% store_COP_b = interp100(store_COP_b2);
% store_COP_t = interp100(store_COP_t2);
store_KJC_b = interp100(store_KJC_b);
store_KJC_t = interp100(store_KJC_t);
store_shankCM_b = interp100(store_shankCM_b);
store_shankCM_t = interp100(store_shankCM_t);
store_pelvis_b = interp100(store_pelvis_b);
store_pelvis_t = interp100(store_pelvis_t);
%%
%%%% DEBUG INTEPR KAM
for i = 1:1:12
    for j = 1:1:10
        % break down into x,y,z
        inp_x = store_KAM_baseline_debug{i,j}(1,:);
        
        % resample to 100
        x = interp1(linspace(1, length(inp_x), length(inp_x)), inp_x, linspace(1, length(inp_x), 100));
        
        % save to resized cell
        store_KAM_baseline_debug{i,j} = x;
        
        
        inp_x = store_KAM_baseline_debug2{i,j}(1,:);
        
        % resample to 100
        x = interp1(linspace(1, length(inp_x), length(inp_x)), inp_x, linspace(1, length(inp_x), 100));
        
        % save to resized cell
        store_KAM_baseline_debug2{i,j} = x;
        
        
        inp_x = store_KAM_baseline_debug3{i,j}(1,:);
        
        % resample to 100
        x = interp1(linspace(1, length(inp_x), length(inp_x)), inp_x, linspace(1, length(inp_x), 100));
        
        % save to resized cell
        store_KAM_baseline_debug3{i,j} = x;
    end
    for j = 1:1:10
        % break down into x,y,z
        inp_x = store_KAM_trial_debug{i,j}(1,:);
        
        % resample to 100
        x = interp1(linspace(1, length(inp_x), length(inp_x)), inp_x, linspace(1, length(inp_x), 100));
        
        % save to resized cell
        store_KAM_trial_debug{i,j} = x;
                % break down into x,y,z
        inp_x = store_KAM_trial_debug2{i,j}(1,:);
        
        % resample to 100
        x = interp1(linspace(1, length(inp_x), length(inp_x)), inp_x, linspace(1, length(inp_x), 100));
        
        % save to resized cell
        store_KAM_trial_debug2{i,j} = x;
        
        
        
        inp_x = store_KAM_trial_debug3{i,j}(1,:);
        
        % resample to 100
        x = interp1(linspace(1, length(inp_x), length(inp_x)), inp_x, linspace(1, length(inp_x), 100));
        
        % save to resized cell
        store_KAM_trial_debug3{i,j} = x;
    end
end
%% 2. relate COP, KJC, shank to pelvis (in x, y, z)

store_KJC_b_p = wrtPelvis(store_KJC_b, store_pelvis_b);
store_KJC_t_p = wrtPelvis(store_KJC_t, store_pelvis_t);
store_COP_b_p = wrtPelvis(store_COP_b, store_pelvis_b);
store_COP_t_p = wrtPelvis(store_COP_t, store_pelvis_t);
store_shankCM_b_p = wrtPelvis(store_shankCM_b, store_pelvis_b);
store_shankCM_t_p = wrtPelvis(store_shankCM_t, store_pelvis_t);

% compute average trajectories at baseline and trial for future validation:
[baseline_pelvis_x, baseline_pelvis_y, baseline_pelvis_z] = subjectAvgTrajectory(store_pelvis_t);
[baseline_GRF_x, baseline_GRF_y, baseline_GRF_z] = subjectAvgTrajectory(store_GRF_b);
[baseline_COP_x, baseline_COP_y, baseline_COP_z] = subjectAvgTrajectory(store_COP_b_p);
[baseline_KJC_x, baseline_KJC_y, baseline_KJC_z] = subjectAvgTrajectory(store_KJC_b_p);
[baseline_shankCM_x, baseline_shankCM_y, baseline_shankCM_z] = subjectAvgTrajectory(store_shankCM_b_p);


[trial_pelvis_x, trial_pelvis_y, trial_pelvis_z] = subjectAvgTrajectory(store_pelvis_t);
[trial_GRF_x, trial_GRF_y, trial_GRF_z] = subjectAvgTrajectory(store_GRF_t);
[trial_COP_x, trial_COP_y, trial_COP_z] = subjectAvgTrajectory(store_COP_t_p);
[trial_KJC_x, trial_KJC_y, trial_KJC_z] = subjectAvgTrajectory(store_KJC_t_p);
[trial_shankCM_x, trial_shankCM_y, trial_shankCM_z] = subjectAvgTrajectory(store_shankCM_t_p);


%% 3. prepare format for linreg 

% a. compute target toe-in FPA
allfpab = cell2mat(store_FPA_b); 
allfpat = cell2mat(store_FPA_t);
fpabmean = mean(allfpab,2);         % find subject's average baseline FPA
fpat = allfpat - fpabmean;          % relate subject's toe-in FPAs to baseline
target_FPA = round(mean(allfpab,2) - mean(allfpat,2));

% % b. find subject average baselines
[baseline_KJC_x, baseline_KJC_y, baseline_KJC_z] = subjectAvgTrajectory(store_KJC_b_p);
[baseline_COP_x, baseline_COP_y, baseline_COP_z] = subjectAvgTrajectory(store_COP_b_p);
[baseline_shankCM_x, baseline_shankCM_y, baseline_shankCM_z] = subjectAvgTrajectory(store_shankCM_b_p);

% c. append toe-in FPA to toe-in COP, KJC, shank
for i = 1:1:12
    for j = 1:1:10
        store_shankCM_t_p{i,j} = [ round(fpat(i,j)), store_shankCM_t_p{i,j}(1,:); round(fpat(i,j)), store_shankCM_t_p{i,j}(2,:); round(fpat(i,j)), store_shankCM_t_p{i,j}(3,:)];
        store_KJC_t_p{i,j} = [ round(fpat(i,j)), store_KJC_t_p{i,j}(1,:); round(fpat(i,j)), store_KJC_t_p{i,j}(2,:); round(fpat(i,j)), store_KJC_t_p{i,j}(3,:)];
        store_COP_t_p{i,j} = [ round(fpat(i,j)), store_COP_t_p{i,j}(1,:); round(fpat(i,j)), store_COP_t_p{i,j}(2,:); round(fpat(i,j)), store_COP_t_p{i,j}(3,:)];
    end
end

store_KJC_t_p_diff= store_KJC_t_p;
store_COP_t_p_diff = store_COP_t_p;
store_shankCM_t_p_diff = store_shankCM_t_p;

% d. find difference between toe-in and baseline trajectories
for ii = 1:1:12
    for j = 1:1:10

        % subtract trajectory (AP and ML)
        store_KJC_t_p_diff{ii,j}(1, 2:end) = store_KJC_t_p{ii,j}(1, 2:end) - baseline_KJC_x(ii, :);
        store_KJC_t_p_diff{ii,j}(2, 2:end) = store_KJC_t_p{ii,j}(2, 2:end) - baseline_KJC_y(ii, :);
        store_KJC_t_p_diff{ii,j}(3, 2:end) = store_KJC_t_p{ii,j}(3, 2:end);
        store_COP_t_p_diff{ii,j}(1, 2:end) = store_COP_t_p{ii,j}(1, 2:end) - baseline_COP_x(ii, :);
        store_COP_t_p_diff{ii,j}(2, 2:end) = store_COP_t_p{ii,j}(2, 2:end) - baseline_COP_y(ii, :);
        store_COP_t_p_diff{ii,j}(3, 2:end) = store_COP_t_p{ii,j}(3, 2:end);
        store_shankCM_t_p_diff{ii,j}(1, 2:end) = store_shankCM_t_p{ii,j}(1, 2:end) - baseline_shankCM_x(ii, :);
        store_shankCM_t_p_diff{ii,j}(2, 2:end) = store_shankCM_t_p{ii,j}(2, 2:end) - baseline_shankCM_y(ii, :);
        store_shankCM_t_p_diff{ii,j}(3, 2:end) = store_shankCM_t_p{ii,j}(3, 2:end);

    end
end

% e. select out x, y, z params
all_shank_t_x = cellfun(@(v)v(1,:), store_shankCM_t_p_diff, 'UniformOutput', false);
all_shank_t_y = cellfun(@(v)v(2,:), store_shankCM_t_p_diff, 'UniformOutput', false);
all_shank_t_z = cellfun(@(v)v(3,:), store_shankCM_t_p_diff, 'UniformOutput', false);
all_KJC_t_x = cellfun(@(v)v(1,:), store_KJC_t_p_diff, 'UniformOutput', false);
all_KJC_t_y = cellfun(@(v)v(2,:), store_KJC_t_p_diff, 'UniformOutput', false);
all_KJC_t_z = cellfun(@(v)v(3,:), store_KJC_t_p_diff, 'UniformOutput', false);
all_COP_t_x = cellfun(@(v)v(1,:), store_COP_t_p_diff, 'UniformOutput', false);
all_COP_t_y = cellfun(@(v)v(2,:), store_COP_t_p_diff, 'UniformOutput', false);
all_COP_t_z = cellfun(@(v)v(3,:), store_COP_t_p_diff, 'UniformOutput', false);

%% LOOCV
for leftout = 1:1:12
    disp(['Beginning subject: ', num2str(leftout)]);
    % a. save as test and training variables
    Tshank_t_x = all_shank_t_x;
    Tshank_t_y = all_shank_t_y;
    Tshank_t_z = all_shank_t_z;
    TKJC_t_x = all_KJC_t_x;
    TKJC_t_y = all_KJC_t_y;
    TKJC_t_z = all_KJC_t_z;
    TCOP_t_x = all_COP_t_x;
    TCOP_t_y = all_COP_t_y;
    TCOP_t_z = all_COP_t_z;
    
    % b. save testing var
    TESTshank_t_x = Tshank_t_x(leftout,:);
    TESTshank_t_y = Tshank_t_y(leftout,:);
    TESTshank_t_z = Tshank_t_z(leftout,:);
    TESTKJC_t_x = TKJC_t_x(leftout,:);
    TESTKJC_t_y = TKJC_t_y(leftout,:);
    TESTKJC_t_z = TKJC_t_z(leftout,:);
    TESTCOP_t_x = TCOP_t_x(leftout,:);
    TESTCOP_t_y = TCOP_t_y(leftout,:);
    TESTCOP_t_z = TCOP_t_z(leftout,:);
    
    % c. remove testing from training
    % NR 08.25.2020 commented this to save 5deg heuristics:
    Tshank_t_x(leftout,:) = [];
    Tshank_t_y(leftout,:) = [];
    Tshank_t_z(leftout,:) = [];
    TKJC_t_x(leftout,:) = [];
    TKJC_t_y(leftout,:) = [];
    TKJC_t_z(leftout,:) = [];
    TCOP_t_x(leftout,:) = [];
    TCOP_t_y(leftout,:) = [];
    TCOP_t_z(leftout,:) = [];
    
    % d. sort training data as ascending
    Tshank_t_x = sortrows(vertcat(Tshank_t_x{:}),1);
    Tshank_t_y = sortrows(vertcat(Tshank_t_y{:}),1);
    Tshank_t_z = sortrows(vertcat(Tshank_t_z{:}),1);
    TKJC_t_x = sortrows(vertcat(TKJC_t_x{:}),1);
    TKJC_t_y = sortrows(vertcat(TKJC_t_y{:}),1);
    TKJC_t_z = sortrows(vertcat(TKJC_t_z{:}),1);
    TCOP_t_x = sortrows(vertcat(TCOP_t_x{:}),1);
    TCOP_t_y = sortrows(vertcat(TCOP_t_y{:}),1);
    TCOP_t_z = sortrows(vertcat(TCOP_t_z{:}),1);
    
    % e. remove any rows greater than 0 or less than -14
    Tshank_t_x = Tshank_t_x(Tshank_t_x(:,1)<0,:);
    Tshank_t_y = Tshank_t_y(Tshank_t_y(:,1)<0,:);
    Tshank_t_z = Tshank_t_z(Tshank_t_z(:,1)<0,:);
    Tshank_t_x = Tshank_t_x(Tshank_t_x(:,1)>-14,:);
    Tshank_t_y = Tshank_t_y(Tshank_t_y(:,1)>-14,:);
    Tshank_t_z = Tshank_t_z(Tshank_t_z(:,1)>-14,:);
    TKJC_t_x = TKJC_t_x(TKJC_t_x(:,1)<0,:);
    TKJC_t_y = TKJC_t_y(TKJC_t_y(:,1)<0,:);
    TKJC_t_z = TKJC_t_z(TKJC_t_z(:,1)<0,:);
    TKJC_t_x = TKJC_t_x(TKJC_t_x(:,1)>-14,:);
    TKJC_t_y = TKJC_t_y(TKJC_t_y(:,1)>-14,:);
    TKJC_t_z = TKJC_t_z(TKJC_t_z(:,1)>-14,:);
    TCOP_t_x = TCOP_t_x(TCOP_t_x(:,1)<0,:);
    TCOP_t_y = TCOP_t_y(TCOP_t_y(:,1)<0,:);
    TCOP_t_z = TCOP_t_z(TCOP_t_z(:,1)<0,:);
    TCOP_t_x = TCOP_t_x(TCOP_t_x(:,1)>-14,:);
    TCOP_t_y = TCOP_t_y(TCOP_t_y(:,1)>-14,:);
    TCOP_t_z = TCOP_t_z(TCOP_t_z(:,1)>-14,:);
    
    disp(['Number of outlier rows removed: ', num2str(110 - length(TCOP_t_z))]);
    
    % f. group into bins
    store_shankx_at_FPA = {};
    store_shanky_at_FPA = {};
    store_shankz_at_FPA = {};
    store_KJCx_at_FPA = {};
    store_KJCy_at_FPA = {};
    store_KJCz_at_FPA = {};
    store_COPx_at_FPA = {};
    store_COPy_at_FPA = {};
    store_COPz_at_FPA = {};
    for ii = max(Tshank_t_x(:,1)) : -1: -10  %min(Tshank_t_x(:,1))
        iter = 1;
        for j = 1:1:size(Tshank_t_x,1)
            if Tshank_t_x(j, 1) == ii
                store_shankx_at_FPA{-ii, iter} = Tshank_t_x(j,:);
                store_shanky_at_FPA{-ii, iter} = Tshank_t_y(j,:);
                store_shankz_at_FPA{-ii, iter} = Tshank_t_z(j,:);
                store_KJCx_at_FPA{-ii, iter} = TKJC_t_x(j,:);
                store_KJCy_at_FPA{-ii, iter} = TKJC_t_y(j,:);
                store_KJCz_at_FPA{-ii, iter} = TKJC_t_z(j,:);
                store_COPx_at_FPA{-ii, iter} = TCOP_t_x(j,:);
                store_COPy_at_FPA{-ii, iter} = TCOP_t_y(j,:);
                store_COPz_at_FPA{-ii, iter} = TCOP_t_z(j,:);
                iter = iter+1;
            end
        end
    end
    % i. find average at each FPA
    % (capped at 10 deg, uncomment if you want full range)
    shape = 10; %-min(Tshank_t_x(:,1));
    
    y_shank = zeros(shape,100);
    x_shank = zeros(shape,100);
    z_shank = zeros(shape,100);
    y_KJC = zeros(shape,100);
    x_KJC = zeros(shape,100);
    z_KJC = zeros(shape,100);
    y_COP = zeros(shape,100);
    x_COP = zeros(shape,100);
    z_COP = zeros(shape,100);
    
    for ii = 1:1:10%-min(Tshank_t_x(:,1))
        arranged = cell2mat(store_shanky_at_FPA(ii,:));
        num = length(arranged)/101;
        all_shanky = reshape(arranged, [101, num])';
        
        arranged = cell2mat(store_shankx_at_FPA(ii,:));
        num = length(arranged)/101;
        all_shankx = reshape(arranged, [101, num])';
        
        arranged = cell2mat(store_shankz_at_FPA(ii,:));
        num = length(arranged)/101;
        all_shankz = reshape(arranged, [101, num])';
        
        arranged = cell2mat(store_KJCx_at_FPA(ii,:));
        num = length(arranged)/101;
        all_KJCx = reshape(arranged, [101, num])';
        
        arranged = cell2mat(store_KJCy_at_FPA(ii,:));
        num = length(arranged)/101;
        all_KJCy = reshape(arranged, [101, num])';
        
        arranged = cell2mat(store_KJCz_at_FPA(ii,:));
        num = length(arranged)/101;
        all_KJCz = reshape(arranged, [101, num])';
        
        arranged = cell2mat(store_COPx_at_FPA(ii,:));
        num = length(arranged)/101;
        all_COPx = reshape(arranged, [101, num])';
        
        arranged = cell2mat(store_COPy_at_FPA(ii,:));
        num = length(arranged)/101;
        all_COPy = reshape(arranged, [101, num])';
        
        arranged = cell2mat(store_COPz_at_FPA(ii,:));
        num = length(arranged)/101;
        all_COPz = reshape(arranged, [101, num])';
        
        % j. find and save value at each point over 100
        for j = 1:1:100
            y_shank(ii, j) = mean(all_shanky(:, j+1));
            x_shank(ii, j) = mean(all_shankx(:, j+1));
            z_shank(ii, j) = mean(all_shankz(:, j+1));
            y_KJC(ii, j) = mean(all_KJCy(:, j+1));
            x_KJC(ii, j) = mean(all_KJCx(:, j+1));
            z_KJC(ii, j) = mean(all_KJCz(:, j+1));
            y_COP(ii, j) = mean(all_COPy(:, j+1));
            x_COP(ii, j) = mean(all_COPx(:, j+1));
            z_COP(ii, j) = mean(all_COPz(:, j+1));
        end
    end
    
    % k. fit a regression line at each bin
    y_pred_shank = zeros(shape,100);
    x_pred_shank = zeros(shape,100);
    z_pred_shank = zeros(shape,100);
    y_pred_KJC = zeros(shape,100);
    x_pred_KJC = zeros(shape,100);
    z_pred_KJC = zeros(shape,100);
    y_pred_COP = zeros(shape,100);
    x_pred_COP = zeros(shape,100);
    z_pred_COP = zeros(shape,100);
    
    for ii = 1:1:100
        mdl = fitlm(linspace(1,shape,shape), y_shank(:,ii));
        y_pred_shank(:, ii) = mdl.Fitted;
        %y_pred_shank_R2(i,ii) = mdl.Rsquared.Ordinary;
        
        mdl = fitlm(linspace(1,shape,shape), x_shank(:,ii));
        x_pred_shank(:, ii) = mdl.Fitted;
        %x_pred_shank_R2(i,ii) = mdl.Rsquared.Ordinary;
        
        mdl = fitlm(linspace(1,shape,shape), z_shank(:,ii));
        z_pred_shank(:,ii) = mdl.Fitted;
        %z_pred_shank_R2(i,ii) = mdl.Rsquared.Ordinary;
        
        mdl = fitlm(linspace(1,shape,shape), y_KJC(:,ii));
        y_pred_KJC(:, ii) = mdl.Fitted;
        y_pred_KJC_R2(leftout,ii) = mdl.Rsquared.Ordinary;
        
        mdl = fitlm(linspace(1,shape,shape), x_KJC(:,ii));
        x_pred_KJC(:, ii) = mdl.Fitted;
        x_pred_KJC_R2(leftout,ii) = mdl.Rsquared.Ordinary;
        
        mdl = fitlm(linspace(1,shape,shape), z_KJC(:,ii));
        z_pred_KJC(:, ii) = mdl.Fitted;
        z_pred_KJC_R2(leftout,ii) = mdl.Rsquared.Ordinary;
        
        mdl = fitlm(linspace(1,shape,shape), y_COP(:,ii));
        y_pred_COP(:, ii) = mdl.Fitted;
        y_pred_COP_R2(leftout,ii) = mdl.Rsquared.Ordinary;
        
        mdl = fitlm(linspace(1,shape,shape), x_COP(:,ii));
        x_pred_COP(:, ii) = mdl.Fitted;
        x_pred_COP_R2(leftout,ii) = mdl.Rsquared.Ordinary;
        
        mdl = fitlm(linspace(1,shape,shape), z_COP(:,ii));
        z_pred_COP(:, ii) = mdl.Fitted;
        z_pred_COP_R2(leftout,ii) = mdl.Rsquared.Ordinary;
        
    end

    disp(['Completed training subject: ', num2str(leftout)]);
    
    % l. generate predicted trajectories for leftout subject
    fpa = target_FPA(leftout);       
   
    
    if fpa > shape % if target FPA is outside of training range then set as max
        fpa = shape;
        disp('Target FPA outside of training range, setting to max FPA');
    end
    
    % pred x and y-component use baseline for  z
    pred_kjc = [x_pred_KJC(fpa,:); y_pred_KJC(fpa,:); baseline_KJC_z(leftout,:)];
    pred_cop = [x_pred_COP(fpa,:); y_pred_COP(fpa,:); baseline_COP_z(leftout,:)];
    pred_shank = [x_pred_shank(fpa,:); y_pred_shank(fpa,:); baseline_shankCM_z(leftout,:)];
    pred_grf = [baseline_GRF_x(leftout,:); baseline_GRF_y(leftout,:); baseline_GRF_z(leftout,:)];
    pred_pelvis = [baseline_pelvis_x(leftout,:); baseline_pelvis_y(leftout,:); baseline_pelvis_z(leftout,:)];
    
    % add back baseline_XXXX_XXX
    pred_kjc = [pred_kjc(1,:) + baseline_KJC_x(leftout,:); pred_kjc(2,:) + baseline_KJC_y(leftout,:); pred_kjc(3,:)];
    pred_cop = [pred_cop(1,:) + baseline_COP_x(leftout,:); pred_cop(2,:) + baseline_COP_y(leftout,:); pred_cop(3,:)];
    pred_shank = [pred_shank(1,:) + baseline_shankCM_x(leftout,:); pred_shank(2,:) + baseline_shankCM_y(leftout,:); pred_shank(3,:)];
    
    % relate all back to global by adding back pelvis 
    pred_kjc_g = pred_kjc + pred_pelvis;
    pred_cop_g = pred_cop + pred_pelvis;
    pred_shank_g = pred_shank + pred_pelvis;
    
    % relate kjc and cop to shank
    pred_kjc_l = pred_kjc_g - pred_shank_g;
    pred_cop_l = [pred_cop_g(1,:) - pred_shank_g(1,:); pred_cop_g(2,:) - pred_shank_g(2,:); zeros(1,100)];
    
    %%% NROKH DEBUG 08.27.2020
    figure()
    hold on
    r = pred_cop_g/1000 - pred_kjc_g/1000;
    KM = cross(r, pred_grf);
    %plot(-KM(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)*9.8)*100)
    plot(  (-KM(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100*9.8))*100)
    
    
    
    
    
    
%     %
%      figure()
%    % subplot(4,3,[1:3])
%     KM = cross(pred_cop_l - pred_grf, pred_kjc_l - pred_grf);
%     plot(KM(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)));
%     hold on
    
    
    
    
    %%%% COMPARE TO REAL TOE-IN
    real_kjc = [trial_KJC_x(leftout,:); trial_KJC_y(leftout,:); trial_KJC_z(leftout,:)];
    real_cop = [trial_COP_x(leftout,:); trial_COP_y(leftout,:); trial_COP_z(leftout,:)];
    real_shank = [trial_shankCM_x(leftout,:); trial_shankCM_y(leftout,:); trial_shankCM_z(leftout,:)];
    real_grf = [trial_GRF_x(leftout,:); trial_GRF_y(leftout,:); trial_GRF_z(leftout,:)];
    real_pelvis = [trial_pelvis_x(leftout,:); trial_pelvis_y(leftout,:); trial_pelvis_z(leftout,:)];
    real_kjc_g = real_kjc + real_pelvis;
    real_cop_g = real_cop + real_pelvis;
    real_shank_g = real_shank + real_pelvis;
    real_kjc_l = real_kjc_g - real_shank_g;
    real_cop_l = [real_cop_g(1,:) - real_shank_g(1,:); real_cop_g(2,:) - real_shank_g(2,:); zeros(1,100)];
    
%     KM2 = cross(real_cop_l - real_grf, real_kjc_l - real_grf);
%     plot(KM2(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)));
    r = real_cop_g/1000 - real_kjc_g/1000;
    KM2 = cross(r, real_grf);
    %plot(-KM2(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100)*10)
    plot(  (-KM2(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100*9.8))*100)
    
    %%%% COMPARE TO REAL BASELINE
    base_kjc = [baseline_KJC_x(leftout,:); baseline_KJC_y(leftout,:); baseline_KJC_z(leftout,:)];
    base_cop = [baseline_COP_x(leftout,:); baseline_COP_y(leftout,:); baseline_COP_z(leftout,:)];
    base_shank = [baseline_shankCM_x(leftout,:); baseline_shankCM_y(leftout,:); baseline_shankCM_z(leftout,:)];
    base_grf = [baseline_GRF_x(leftout,:); baseline_GRF_y(leftout,:); baseline_GRF_z(leftout,:)];
    base_pelvis = [baseline_pelvis_x(leftout,:); baseline_pelvis_y(leftout,:); baseline_pelvis_z(leftout,:)];
    base_kjc_g = base_kjc + base_pelvis;
    base_cop_g = base_cop + base_pelvis;
    base_shank_g = base_shank + base_pelvis;
    base_kjc_l = base_kjc_g - base_shank_g;
    base_cop_l = [base_cop_g(1,:) - base_shank_g(1,:); base_cop_g(2,:) - base_shank_g(2,:); zeros(1,100)];
    
        r = base_cop_g/1000 - base_kjc_g/1000;
    KM3 = cross(r, base_grf);
    %plot(-KM3(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100)*10)
    plot(  (-KM3(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100*9.8))*100)
%     KM3 = cross(base_cop_l - base_grf, base_kjc_l - base_grf);
%     plot(KM3(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)));
    legend('predicted toe-in KAM', 'real toe-in KAM', 'baseline KAM');
    xlabel('% stance')
    ylabel('KAM [%BW*H]')
    
    figure()
    % compare predicted params
    subplot(4,3,4)
    plot(real_kjc(1,:), 'r'); hold on; plot(real_kjc(2,:), 'k'); plot(real_kjc(3,:),'b'); plot(pred_kjc(1,:), 'r--');plot(pred_kjc(2,:), 'k--');plot(pred_kjc(3,:), 'b--');
    title('KJC real vs pred--')
    xlabel('% stance')
    ylabel('position [mm]')
    
    subplot(4,3,5)
    plot(real_cop(1,:), 'r'); hold on; plot(real_cop(2,:), 'k'); plot(real_cop(3,:),'b'); plot(pred_cop(1,:), 'r--');plot(pred_cop(2,:), 'k--');plot(pred_cop(3,:), 'b--');
    title('COP real vs pred--')
    
    subplot(4,3,6)
    plot(real_shank(1,:), 'r'); hold on; plot(real_shank(2,:), 'k'); plot(real_shank(3,:),'b'); plot(pred_shank(1,:), 'r--');plot(pred_shank(2,:), 'k--');plot(pred_shank(3,:), 'b--');
    title('shank real vs pred--')
    
    % predicted global: after adding back pelvis
    subplot(4,3,7)
    plot(real_kjc_g(1,:), 'r'); hold on; plot(real_kjc_g(2,:), 'k'); plot(real_kjc_g(3,:),'b'); plot(pred_kjc_g(1,:), 'r--');plot(pred_kjc_g(2,:), 'k--');plot(pred_kjc_g(3,:), 'b--');
    title('global KJC real vs pred--')
    
    subplot(4,3,8)
    plot(real_cop_g(1,:), 'r'); hold on; plot(real_cop_g(2,:), 'k'); plot(real_cop_g(3,:),'b'); plot(pred_cop_g(1,:), 'r--');plot(pred_cop_g(2,:), 'k--');plot(pred_cop_g(3,:), 'b--');
    title('global COP real vs pred--')
    
    subplot(4,3,9)
    plot(real_shank_g(1,:), 'r'); hold on; plot(real_shank_g(2,:), 'k'); plot(real_shank_g(3,:),'b'); plot(pred_shank_g(1,:), 'r--');plot(pred_shank_g(2,:), 'k--');plot(pred_shank_g(3,:), 'b--');
    title('global shank real vs pred--')
    
    % predicted local: after subtracting shank
    subplot(4,3,10)
    plot(real_kjc_l(1,:), 'r'); hold on; plot(real_kjc_l(2,:), 'k'); plot(real_kjc_l(3,:),'b'); plot(pred_kjc_l(1,:), 'r--');plot(pred_kjc_l(2,:), 'k--');plot(pred_kjc_l(3,:), 'b--');
    title('local KJC real vs pred--')
    
    subplot(4,3,11)
    plot(real_cop_l(1,:), 'r'); hold on; plot(real_cop_l(2,:), 'k'); plot(real_cop_l(3,:),'b'); plot(pred_cop_l(1,:), 'r--');plot(pred_cop_l(2,:), 'k--');plot(pred_cop_l(3,:), 'b--');
    title('local COP real vs pred--')
    
    subplot(4,3,12)
    plot(real_grf(1,:), 'r'); hold on; plot(real_grf(2,:), 'k'); plot(real_grf(3,:), 'b'); plot(pred_grf(1,:), 'r--'); plot(pred_grf(2,:), 'k--'); plot(pred_grf(3,:), 'b--');
    title('GRF real vs pred--')
    xlabel('% stance')
    ylabel('GRF [N]')
    
    sgtitle(['S', num2str(leftout), ' diagnostics']);
    legend('x', 'y', 'z', 'pred x', 'pred y', 'pred z');
    
    % compute error
    RMSE_KJC_x(leftout) = sqrt(mean( (pred_kjc(1,:) - real_kjc(1,:)).^2));
    RMSE_KJC_y(leftout) = sqrt(mean( (pred_kjc(2,:) - real_kjc(2,:)).^2));
    RMSE_KJC_z(leftout) = sqrt(mean( (pred_kjc(3,:) - real_kjc(3,:)).^2));
    
    RMSE_COP_x(leftout) = sqrt(mean( (pred_cop(1,:) - real_cop(1,:)).^2));
    RMSE_COP_y(leftout) = sqrt(mean( (pred_cop(2,:) - real_cop(2,:)).^2));
    RMSE_COP_z(leftout) = sqrt(mean( (pred_cop(3,:) - real_cop(3,:)).^2));
    
    RMSE_shank_x(leftout) = sqrt(mean( (pred_shank(1,:) - real_shank(1,:)).^2));
    RMSE_shank_y(leftout) = sqrt(mean( (pred_shank(2,:) - real_shank(2,:)).^2));
    RMSE_shank_z(leftout) = sqrt(mean( (pred_shank(3,:) - real_shank(3,:)).^2));
    
    % comparative RMSE (what's the difference between real and baseline toe-in)
    cRMSE_KJC_x(leftout) = sqrt(mean( (real_kjc_g(1,:) - base_kjc_g(1,:)).^2));
    cRMSE_KJC_y(leftout) = sqrt(mean( (real_kjc_g(2,:) - base_kjc_g(2,:)).^2));
    cRMSE_KJC_z(leftout) = sqrt(mean( (real_kjc_g(3,:) - base_kjc_g(3,:)).^2));
    
    cRMSE_COP_x(leftout) = sqrt(mean( (real_cop(1,:) - base_cop(1,:)).^2));
    cRMSE_COP_y(leftout) = sqrt(mean( (real_cop(2,:) - base_cop(2,:)).^2));
    cRMSE_COP_z(leftout) = sqrt(mean( (real_cop(3,:) - base_cop(3,:)).^2));
    
    cRMSE_shank_x(leftout) = sqrt(mean( (real_shank(1,:) - base_shank(1,:)).^2));
    cRMSE_shank_y(leftout) = sqrt(mean( (real_shank(2,:) - base_shank(2,:)).^2));
    cRMSE_shank_z(leftout) = sqrt(mean( (real_shank(3,:) - base_shank(3,:)).^2));
    
    % compute KAM RMSE
    KAMpred(leftout,:) = (KM(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100)*10);
    KAMreal(leftout,:) = (KM2(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100)*10);
    KAMbaseline(leftout,:) = (KM3(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100)*10);
    RMSE_KAM(leftout) = sqrt (mean( (KAMpred(leftout,:) - KAMreal(leftout,:)).^2));
    RMSE_KAM_compare(leftout) = sqrt (mean( (KAMbaseline(leftout,:) - KAMreal(leftout,:)).^2));
    

    

end

%% determine if first peak decreased for each subject
psor1 = [];
lsor1 = [];
psor2 = [];
lsor2 = [];
peaks_KAM = cell(12,1);
for i = 1:1:12
    % compute peak KAM decreases
    % peak detect
   
    x = linspace(1,100,100);
    % find peaks
    %[psor1, lsor1] = findpeaks(-KAMbaseline(i,:), x, 'SortStr', 'descend');%, 'MinPeakDistance', 30, 'MinPeakProminence', 2);
    [psor1, lsor1] = findpeaks(-KAMpred(i,:), x, 'SortStr', 'descend');%, 'MinPeakDistance', 30, 'MinPeakProminence', 2);
    [psor2, lsor2] = findpeaks(-KAMreal(i,:), x, 'SortStr', 'descend');%, 'MinPeakDistance', 30);
    
    % save just first two peaks
    peaks_KAM{i} = [psor1(1), lsor1(1)];%; psor1(2), lsor1(2)];
    peaks_tKAM{i} = [psor2(1), lsor2(1); psor2(2), lsor2(2)];
end

firstPeakDecrease = zeros(12,1);
secondPeakIncrease = zeros(12,1);
for i = 1:1:12
    %disp(['              SUBJECT: ', num2str(i)])
    [a, firstPeakloc] = min(peaks_KAM{i,1}(:,2)); 
    firstPeakKAM = peaks_KAM{i,1}(firstPeakloc,:); % first peak value
    [b, secondPeakloc] = max(peaks_KAM{i,1}(:,2)); 
    secondPeakKAM = peaks_KAM{i,1}(secondPeakloc,:); % second peak value
    
    [c, firstPeakloc] = min(peaks_tKAM{i,1}(:,2)); 
    firstPeakSKAM = peaks_tKAM{i,1}(firstPeakloc,:); % first peak value
    [d, secondPeakloc] = max(peaks_tKAM{i,1}(:,2)); 
    secondPeakSKAM = peaks_tKAM{i,1}(secondPeakloc,:); % second peak value
    
    
    % compare first peak vals
    if firstPeakSKAM(1) < firstPeakKAM(1)
        disp(['first peak decreased by ', num2str(firstPeakKAM(1) - firstPeakSKAM(1))])
        
        firstPeakDecrease(i) = firstPeakKAM(1) - firstPeakSKAM(1);
        
    elseif firstPeakSKAM(1) > firstPeakKAM(1)
        disp(['first peak increased by ', num2str(firstPeakSKAM(1) - firstPeakKAM(1))])
        firstPeakDecrease(i) = firstPeakKAM(1) - firstPeakSKAM(1);
    else
        disp('Equal????')
    end
    
    % compare second peak vals
    if secondPeakSKAM(1) < secondPeakKAM(1)
        disp(['second peak decreased by ', num2str(secondPeakKAM(1) - secondPeakSKAM(1))]);
    elseif secondPeakSKAM(1) > secondPeakKAM(1)
        disp(['second peak increased by ', num2str(secondPeakSKAM(1) - secondPeakKAM(1))]);
        secondPeakIncrease(i) = secondPeakSKAM(1) - secondPeakKAM(1);
    else
        disp('Equal????')
    end

    % check if second peak increase is greater than first peak decrease
    if secondPeakIncrease > firstPeakDecrease
        disp('            WARNING second peak increased more than first peak') % both s32 and 102 aren't actually a case for this
    end
end

%% recalculating first peak error
for i = 1:1:12
    figure();
    hold on
    [pks1, locs1 ] = findpeaks(-KAMpred(i,:), 'MinPeakDistance', 30, 'MinPeakProminence', 0);
    [pks2, locs2 ] = findpeaks(-KAMreal(i,:), 'MinPeakDistance', 30, 'MinPeakProminence', 0);
    [pks3, locs3 ] = findpeaks(KAM_fda(i,:), 'MinPeakDistance', 30, 'MinPeakProminence', 0);
    findpeaks(-KAMpred(i,:), 'MinPeakDistance', 30, 'MinPeakProminence', 0);
    findpeaks(-KAMreal(i,:), 'MinPeakDistance', 30, 'MinPeakProminence', 0);
    findpeaks(KAM_fda(i,:), 'MinPeakDistance', 30, 'MinPeakProminence', 0);
    text(locs1-1, pks1-1, num2str(pks1))
    text(locs2+1, pks2+0.5, num2str(pks2))
    text(locs3+1, pks3-3, num2str(pks3))
end
%% RESULTS:
% RMSE COP x: 13.8774mm +- 5.4333 (vs 14.7552 mm
% RMSE COP y: 8.1251 mm +- 5.5823(vs 7.8800 mm
% RMSE KJC x: 10.1361 mm +_ 5.1334(vs 8.6452 mm
% RMSE KJC y: 4.9646 mm +_ 2.6756 (vs 4.7828
% RMSE shank x: 9.5051 mm +_ 4.5808(vs 7.9355
% RMSE shank y: 5.8343 mm +_ 4.2048(vs 5.5354
% RMSE KAM: 0.1746 +- 0.0912
% RMSE FP KAM: 0.1732 +-0.1663



figure()
hold on


p1 = plot(-mean(KAMpred), 'LineWidth', 2);
p2 = plot(-mean(KAMreal), '--', 'LineWidth', 2);

p3 = plot(-mean(KAMbaseline), '-.', 'LineWidth', 2);

p1.Color = [0 0 0];
p2.Color = [239,138,98]/256;
p3.Color = [103,169,207]/256;


box off
set(gca,'TickLength',[0 0])
set(gca,'FontSize',17)

grid on
ax = gca;
ax.GridColor = [0.01, 0.01, 0.01];  % [R, G, B]
ax.GridAlpha = 0.075;
ylim([-1, 3])


%% KAM VISUALIZATION 04.29.2021 (wow it's been a minute since I was here last)
figure()
hold on

c10 = -(mean(KAMpred) + std(KAMpred));
c11 = -(mean(KAMpred) - std(KAMpred));
x = 1:length(c10);
x2 = [x, fliplr(x)];
inbetween = [c10, fliplr(c11)];
color = [1 0 0];
a2 = fill(x2, inbetween, color);
a2.FaceAlpha = 0.5;
a2.EdgeAlpha = 0;


c20 = -(mean(KAMreal) + std(KAMreal));
c21 = -(mean(KAMreal) - std(KAMreal));
inbetween = [c20, fliplr(c21)];
color = [0 1 0];
a2 = fill(x2, inbetween, color);
a2.FaceAlpha = 0.5;
a2.EdgeAlpha = 0;

c30 = -(mean(KAMbaseline) + std(KAMbaseline));
c31 = -(mean(KAMbaseline) - std(KAMbaseline));
inbetween = [c30, fliplr(c31)];
color = [0 0 1];
a2 = fill(x2, inbetween, color);
a2.FaceAlpha = 0.5;
a2.EdgeAlpha = 0;


p3 = plot(-mean(KAMbaseline));

p2 = plot(-mean(KAMreal));
p1 = plot(-mean(KAMpred));

p1.Color = [0 0 0];
p2.Color = [239,138,98]/256;
p3.Color = [103,169,207]/256;


%%

plot(x_pred_COP(:,15), 'LineWidth', 2); hold on; plot(x_COP(:,15),'o', 'LineWidth', 2); box off
set(gca,'TickLength',[0 0])
set(gca,'FontSize',17); xlim([0, 11])

%%
figure()
hold on
i = 11
curve1 = -(KAMpred(i,:) + RMSE_KAM(i));
curve2 = -(KAMpred(i,:) - RMSE_KAM(i));
x = 1:length(curve1);
x2 = [x, fliplr(x)];
inbetween = [curve1, fliplr(curve2)];
color = [0, 0.1, 0.1];
a = fill(x2, inbetween, color);
a.FaceAlpha = 0.1;
a.EdgeAlpha = 0;
p1 = plot(-KAMpred(i,:));
p2 = plot(-KAMreal(i,:), '--', 'LineWidth', 2);
p3 = plot(-KAMbaseline(i,:), '-.', 'LineWidth', 2);

p1.Color = [0 0 0];
p2.Color = [239,138,98]/256;
p3.Color = [103,169,207]/256;


box off
set(gca,'TickLength',[0 0])
set(gca,'FontSize',17)

grid on
ax = gca;
ax.GridColor = [0.01, 0.01, 0.01];  % [R, G, B]
ax.GridAlpha = 0.075;

ylim([-1, 5])


%% 05.25.2021 testing FDA prediction on s01
addpath 'D:\sKAM Classifier (2021)\FDA'
load('KJCy_FDA.mat')
load('KJCx_FDA.mat')
load('COPy_FDA.mat')
load('COPx_FDA.mat')




for leftout=1:1:12
    fpa = target_FPA(leftout);
    % linreg pred:
    pred_kjc = [x_pred_KJC(fpa,:); y_pred_KJC(fpa,:); baseline_KJC_z(leftout,:)];
    pred_cop = [x_pred_COP(fpa,:); y_pred_COP(fpa,:); baseline_COP_z(leftout,:)];
    pred_shank = [x_pred_shank(fpa,:); y_pred_shank(fpa,:); baseline_shankCM_z(leftout,:)];
    pred_grf = [baseline_GRF_x(leftout,:); baseline_GRF_y(leftout,:); baseline_GRF_z(leftout,:)];
    pred_pelvis = [baseline_pelvis_x(leftout,:); baseline_pelvis_y(leftout,:); baseline_pelvis_z(leftout,:)];
    
    % add back baseline_XXXX_XXX
    pred_kjc = [pred_kjc(1,:) + baseline_KJC_x(leftout,:); pred_kjc(2,:) + baseline_KJC_y(leftout,:); pred_kjc(3,:)];
    pred_cop = [pred_cop(1,:) + baseline_COP_x(leftout,:); pred_cop(2,:) + baseline_COP_y(leftout,:); pred_cop(3,:)];
    pred_shank = [pred_shank(1,:) + baseline_shankCM_x(leftout,:); pred_shank(2,:) + baseline_shankCM_y(leftout,:); pred_shank(3,:)];
    
    % relate all back to global by adding back pelvis
    pred_kjc_g = pred_kjc + pred_pelvis;
    pred_cop_g = pred_cop + pred_pelvis;
    pred_shank_g = pred_shank + pred_pelvis;
    
    
    %%% NROKH DEBUG 08.27.2020
    figure()
    hold on
    r = pred_cop_g/1000 - pred_kjc_g/1000;
    KM0 = cross(r, pred_grf);
    plot(-KM0(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100)*10)
    KAM_lin(leftout,:) = (-KM0(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100)*10);
    
    % FDA pred:
    pred_kjc = [KJCx_offset(leftout)' + baseline_KJC_x(leftout,:); KJCy_offset(leftout)' + baseline_KJC_y(leftout,:); pred_kjc(3,:)];
    pred_cop = [COPx_offset(leftout)' + baseline_COP_x(leftout,:); COPy_offset(leftout)' + baseline_COP_y(leftout,:); pred_cop(3,:)];
    % relate all back to global by adding back pelvis
    pred_kjc_g = pred_kjc + pred_pelvis;
    pred_cop_g = pred_cop + pred_pelvis;
    pred_shank_g = pred_shank + pred_pelvis;
    
    r = pred_cop_g/1000 - pred_kjc_g/1000;
    KM = cross(r, pred_grf);
    
    plot(  (-KM(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100*9.8))*100)
    KAM_fda(leftout,:) =  (-KM(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100*9.8))*100;
    
    
    
    
    %%%% COMPARE TO REAL TOE-IN
    real_kjc = [trial_KJC_x(leftout,:); trial_KJC_y(leftout,:); trial_KJC_z(leftout,:)];
    real_cop = [trial_COP_x(leftout,:); trial_COP_y(leftout,:); trial_COP_z(leftout,:)];
    real_shank = [trial_shankCM_x(leftout,:); trial_shankCM_y(leftout,:); trial_shankCM_z(leftout,:)];
    real_grf = [trial_GRF_x(leftout,:); trial_GRF_y(leftout,:); trial_GRF_z(leftout,:)];
    real_pelvis = [trial_pelvis_x(leftout,:); trial_pelvis_y(leftout,:); trial_pelvis_z(leftout,:)];
    real_kjc_g = real_kjc + real_pelvis;
    real_cop_g = real_cop + real_pelvis;
    real_shank_g = real_shank + real_pelvis;
    real_kjc_l = real_kjc_g - real_shank_g;
    real_cop_l = [real_cop_g(1,:) - real_shank_g(1,:); real_cop_g(2,:) - real_shank_g(2,:); zeros(1,100)];

    r = real_cop_g/1000 - real_kjc_g/1000;
    KM2 = cross(r, real_grf);

    plot(  (-KM2(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100*9.8))*100)
    KAM_toein(leftout,:) = (-KM2(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100*9.8))*100;
    
    %%%% COMPARE TO REAL BASELINE
    base_kjc = [baseline_KJC_x(leftout,:); baseline_KJC_y(leftout,:); baseline_KJC_z(leftout,:)];
    base_cop = [baseline_COP_x(leftout,:); baseline_COP_y(leftout,:); baseline_COP_z(leftout,:)];
    base_shank = [baseline_shankCM_x(leftout,:); baseline_shankCM_y(leftout,:); baseline_shankCM_z(leftout,:)];
    base_grf = [baseline_GRF_x(leftout,:); baseline_GRF_y(leftout,:); baseline_GRF_z(leftout,:)];
    base_pelvis = [baseline_pelvis_x(leftout,:); baseline_pelvis_y(leftout,:); baseline_pelvis_z(leftout,:)];
    base_kjc_g = base_kjc + base_pelvis;
    base_cop_g = base_cop + base_pelvis;
    base_shank_g = base_shank + base_pelvis;
    base_kjc_l = base_kjc_g - base_shank_g;
    base_cop_l = [base_cop_g(1,:) - base_shank_g(1,:); base_cop_g(2,:) - base_shank_g(2,:); zeros(1,100)];
    
    r = base_cop_g/1000 - base_kjc_g/1000;
    KM3 = cross(r, base_grf);
    plot((-KM3(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100*9.8))*100)
    KAM_base(leftout,:) = (-KM3(1,:)/(subjectHeightsWeights(leftout,2)*subjectHeightsWeights(leftout,3)/100*9.8))*100;

    legend('linreg predicted toe-in KAM', 'FDA predicted toe-in KAM', 'real toe-in KAM', 'baseline KAM');
    xlabel('% stance')
    ylabel('KAM [%BW*H]')
    
   
    
end


figure()
hold on
plot(mean(KAM_toein), 'LineWidth', 2)
plot(mean(KAM_base), 'LineWidth', 2)
plot(mean(KAM_lin), 'LineWidth', 2)
plot(mean(KAM_fda), 'LineWidth', 2)
legend('toe in', 'baseline', 'linreg', 'fda')



%% LOOCV for spline regression
 
target_FPA(6)=10;
for i =1:1:12
    fpa = target_FPA(i);
    
    load("KJCy_FDA" + i)
    load("KJCx_FDA" + i)
    load("COPy_FDA" + i)
    load("COPx_FDA" + i)
    
    y_FDA_KJC = zeros(shape,100);
    x_FDA_KJC = zeros(shape,100);
    
    y_FDA_COP = zeros(shape,100);
    x_FDA_COP = zeros(shape,100);

    
    
    % testing if i can remove regressions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_FDA_KJC = KJCx_offset;
    y_FDA_KJC = KJCy_offset;
    x_FDA_COP = COPx_offset;
    y_FDA_COP = COPy_offset;
    
    %%%%%%%%%%%%COMMENT AWAY ^^^ %%%%%%%%%%%%%%%%%%%%%%%%
    
    pred_kjc = [x_pred_KJC(fpa,:); y_pred_KJC(fpa,:); baseline_KJC_z(i,:)];
    pred_cop = [x_pred_COP(fpa,:); y_pred_COP(fpa,:); baseline_COP_z(i,:)];
    pred_grf = [baseline_GRF_x(i,:); baseline_GRF_y(i,:); baseline_GRF_z(i,:)];
    pred_pelvis = [baseline_pelvis_x(i,:); baseline_pelvis_y(i,:); baseline_pelvis_z(i,:)];
    
    % add back baseline_XXXX_XXX
    pred_kjc = [pred_kjc(1,:) + baseline_KJC_x(i,:); pred_kjc(2,:) + baseline_KJC_y(i,:); pred_kjc(3,:)];
    pred_cop = [pred_cop(1,:) + baseline_COP_x(i,:); pred_cop(2,:) + baseline_COP_y(i,:); pred_cop(3,:)];
    pred_shank = [pred_shank(1,:) + baseline_shankCM_x(i,:); pred_shank(2,:) + baseline_shankCM_y(i,:); pred_shank(3,:)];
    
    % relate all back to global by adding back pelvis
    pred_kjc_g = pred_kjc + pred_pelvis;
    pred_cop_g = pred_cop + pred_pelvis;
    
    
    % FDA pred:
    pred_kjc = [x_FDA_KJC(fpa,:) + baseline_KJC_x(i,:); y_FDA_KJC(fpa,:) + baseline_KJC_y(i,:); pred_kjc(3,:)];
    pred_cop = [x_FDA_COP(fpa,:) + baseline_COP_x(i,:); y_FDA_COP(fpa,:) + baseline_COP_y(i,:); pred_cop(3,:)];
    % relate all back to global by adding back pelvis
    pred_kjc_g = pred_kjc + pred_pelvis;
    pred_cop_g = pred_cop + pred_pelvis;
    pred_shank_g = pred_shank + pred_pelvis;
    
    r = pred_cop_g/1000 - pred_kjc_g/1000;
    KM = cross(r, pred_grf);

    KAM_fda(i,:) = (-KM(1,:)/(subjectHeightsWeights(i,2)*subjectHeightsWeights(i,3)/100*9.8))*100;
    
    % 06.24.21 calculating rmse
    real_kjc_g = [trial_KJC_x(i,:); trial_KJC_y(i,:); trial_KJC_z(i,:)] ;
    real_cop_g = [trial_COP_x(i,:); trial_COP_y(i,:); trial_COP_z(i,:)] ;
    
    rmseKJCap(i,1) = sqrt(mean((real_kjc_g(1,[10:95])-pred_kjc(1,[10:95])).^2));
    rmseKJCml(i,1) = sqrt(mean((real_kjc_g(2,[10:95])-pred_kjc(2,[10:95])).^2));
    rmseCOPap(i,1) = sqrt(mean((real_cop_g(1,[10:95])-pred_cop(1,[10:95])).^2));
    rmseCOPml(i,1) = sqrt(mean((real_cop_g(2,[10:95])-pred_cop(2,[10:95])).^2));
    
    
end




figure;
sgtitle('Shull LOOCV accuracy (no regression fit)')
subplot(1,4,1)
boxplot(rmseKJCap);
hold on
plot(ones(12,1), rmseKJCap,  '*')
title(["RMSE KJC AP", "mean = " + mean(rmseKJCap), "std = " + std(rmseKJCap)])
subplot(1,4,2)
boxplot(rmseKJCml);
hold on
plot(ones(12,1), rmseKJCml,  '*')
title(["RMSE KJC ML", "mean = " + mean(rmseKJCml), "std = " + std(rmseKJCml)])
subplot(1,4,3)
boxplot(rmseCOPap);
hold on
plot(ones(12,1), rmseCOPap,  '*')
title(["RMSE COP AP", "mean = " + mean(rmseCOPap), "std = " + std(rmseCOPap)])
subplot(1,4,4)
boxplot(rmseCOPml);
hold on
plot(ones(12,1), rmseCOPml,  '*')
title(["RMSE COP ML", "mean = " + mean(rmseCOPml), "std = " + std(rmseCOPml)])


% calc entire traj RMSE
for i = 1:1:12
    traj_rmse(i) = sqrt(mean((KAM_toein(i,:)-KAM_fda(i,:)).^2));
end
mean(traj_rmse)
std(traj_rmse)

%% Fig 4 with 95%CI
figure()
subplot(1,2,1)
hold on

sem1 = std(KAM_fda)/sqrt(length(KAM_fda));
ts1 = tinv([0.025  0.975],length(KAM_fda)-1);
ci1 = mean(KAM_fda) + ts1'.*sem1;
x = 1:length(ci1);
x2 = [x, fliplr(x)];
inbetween = [ci1(1,:), fliplr(ci1(2,:))];
color = [0.64,0.08,0.18];
a2 = fill(x2, inbetween, color);
a2.FaceAlpha = 0.3;
a2.EdgeAlpha = 0;

sem2 = std(KAM_toein)/sqrt(length(KAM_toein));
ci2 = mean(KAM_toein) + ts1'.*sem2;
inbetween = [ci2(1,:), fliplr(ci2(2,:))];
color = [0.47,0.72,0.77];
a2 = fill(x2, inbetween, color);
a2.FaceAlpha = 0.3;
a2.EdgeAlpha = 0;

sem3 = std(KAM_base)/sqrt(length(KAM_base));
ci3 = mean(KAM_base) + ts1'.*sem3;
inbetween = [ci3(1,:), fliplr(ci3(2,:))];
color = [0.80,0.80,0.80];
a2 = fill(x2, inbetween, color);
a2.FaceAlpha = 0.3;
a2.EdgeAlpha = 0;

p1 = plot(mean(KAM_base), 'LineWidth', 3);
p2 = plot(mean(KAM_toein), 'LineWidth', 3);
p3 = plot(mean(KAM_fda), 'LineWidth', 3);

p2.Color = [0.47,0.72,0.77];
p1.Color = [0 0 0];
p3.Color = [0.64,0.08,0.18];

ylim([-2,4])

subplot(1,2,2)
hold on
xlim([0,4])

plot([0.75,1.25],ones(1,2)*mean(fp_Shull(:,1)),'k:', 'LineWidth', 2)
plot([1.75,2.25],ones(1,2)*mean(fp_Shull(:,2)),'k:', 'LineWidth', 2)
plot([2.75,3.25],ones(1,2)*mean(fp_Shull(:,3)),'k:', 'LineWidth', 2)

sem1 = std(fp_Shull(:,1))/sqrt(length(fp_Shull(:,1)));
ts1 = tinv([0.025  0.975],length(fp_Shull(:,1))-1);
ci1 = mean(fp_Shull(:,1)) + ts1*sem1;
plot([1,1],[ci1(1),ci1(2)],'k', 'LineWidth', 2)
plot([0.875,1.125],[ci1(1),ci1(1)],'k','LineWidth', 2)
plot([0.875,1.125],[ci1(2),ci1(2)],'k','LineWidth', 2)

sem2 = std(fp_Shull(:,2))/sqrt(length(fp_Shull(:,2)));
ci2 = mean(fp_Shull(:,2)) + ts1*sem2;
plot([2,2],[ci2(1),ci2(2)],'k', 'LineWidth', 2)
plot([1.875,2.125],[ci2(1),ci2(1)],'k','LineWidth', 2)
plot([1.875,2.125],[ci2(2),ci2(2)],'k','LineWidth', 2)

sem3 = std(fp_Shull(:,3))/sqrt(length(fp_Shull(:,3)));
ci3 = mean(fp_Shull(:,3)) + ts1*sem3;
plot([3,3],[ci3(1),ci3(2)],'k', 'LineWidth', 2)
plot([2.875,3.125],[ci3(1),ci3(1)],'k','LineWidth', 2)
plot([2.875,3.125],[ci3(2),ci3(2)],'k','LineWidth', 2)

a1 = scatter(ones(12,1),fp_Shull(:,1));
a1.MarkerFaceAlpha = 0.8;
a1.MarkerFaceColor = [0 0 0];
a1.MarkerEdgeAlpha = 0;
a2 = scatter(2*ones(12,1),fp_Shull(:,2));
a2.MarkerFaceAlpha = 0.8;
a2.MarkerFaceColor = [0.47,0.72,0.77];
a2.MarkerEdgeAlpha = 0;
a3 = scatter(3*ones(12,1),fp_Shull(:,3));
a3.MarkerFaceAlpha = 0.8;
a3.MarkerFaceColor = [0.64,0.08,0.18];
a3.MarkerEdgeAlpha = 0;



%% plot real kjc vs pred for a sample sub
figure;
subplot(2,2,1)
hold on
plot(real_kjc_g(1,:))
plot(pred_kjc_g(1,:))
title('KJC AP')
subplot(2,2,2)
hold on
plot(real_kjc_g(2,:))
plot(pred_kjc_g(2,:))
title('KJC ML')
subplot(2,2,3)
hold on
plot(real_cop_g(1,:))
plot(pred_cop_g(1,:))
title('COP AP')
subplot(2,2,4)
hold on
plot(real_cop_g(2,:))
plot(pred_cop_g(2,:))
title('COP ML')
%%
figure()
hold on

plot(mean(KAM_base), 'LineWidth', 2)
plot(mean(KAM_toein), 'LineWidth', 2)
%plot(mean(KAM_lin), 'LineWidth', 2)
plot(mean(KAM_fda), 'LineWidth', 2)
legend('baseline','toe in',   'fda')


% with std
figure()
hold on

c10 = (mean(KAM_toein) + std(KAM_toein));
c11 = (mean(KAM_toein) - std(KAM_toein));
x = 1:length(c10);
x2 = [x, fliplr(x)];
inbetween = [c10, fliplr(c11)];
color = [1 0 0];
a2 = fill(x2, inbetween, color);
a2.FaceAlpha = 0.1;
a2.EdgeAlpha = 0;


c20 = (mean(KAM_base) + std(KAM_base));
c21 = (mean(KAM_base) - std(KAM_base));
inbetween = [c20, fliplr(c21)];
color = [0 1 0];
a2 = fill(x2, inbetween, color);
a2.FaceAlpha = 0.1;
a2.EdgeAlpha = 0;

c30 = (mean(KAM_fda) + std(KAM_fda));
c31 = (mean(KAM_fda) - std(KAM_fda));
inbetween = [c30, fliplr(c31)];
color = [0 0 1];
a2 = fill(x2, inbetween, color);
a2.FaceAlpha = 0.1;
a2.EdgeAlpha = 0;

p1 = plot(mean(KAM_toein), 'LineWidth', 3);
p2 = plot(mean(KAM_base), 'LineWidth', 3);
p3 = plot(mean(KAM_fda), 'LineWidth', 3);

p1.Color = [1 0 0];
p2.Color = [0 1 0];
p3.Color = [0 0 1];





