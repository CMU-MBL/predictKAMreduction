%% generating synthetic kam from calgary data
% nrokh 04.22.2021


%% 1. load data
clear all; 
disp('Loading Calgary walking file...')
load('eniGrfWalking.mat');
load('eniGrfMeta.mat')

%% load LR offsets
load('D:\Shull_FPA\FINAL\RESULTS_11272020\Shull_offsets.mat') % ap/ml_pred_COP/KJC

% load FDA offsets
load('D:\sKAM Classifier (2021)\FDA\offsets\y_FDA_KJC.mat')
load('D:\sKAM Classifier (2021)\FDA\offsets\x_FDA_KJC.mat')
load('D:\sKAM Classifier (2021)\FDA\offsets\y_FDA_COP.mat')
load('D:\sKAM Classifier (2021)\FDA\offsets\x_FDA_COP.mat')

% using raw splines instead of regressing on splines
% load('D:\sKAM Classifier (2021)\FDA\KJCx_FDA.mat')
% load('D:\sKAM Classifier (2021)\FDA\KJCy_FDA.mat')
% load('D:\sKAM Classifier (2021)\FDA\COPx_FDA.mat')
% load('D:\sKAM Classifier (2021)\FDA\COPy_FDA.mat')
x_FDA_KJC = KJCx_offset;
y_FDA_KJC = KJCy_offset;
x_FDA_COP = COPx_offset;
y_FDA_COP = COPy_offset;

load('store_alignment.mat')
disp('Files loaded')

%% 2. remove subjects with poor calibration accuracy
disp('Loading calibration file')
load('C:\Users\nrokh\Box\CMU_MBL\Data\Calgary\GRF\grf_calibrations_20210309.mat')

% a. remove empty row at top of ...
store_alignment(1,:) = [];
walkRunList(1,:) = [];

events(1,:) = [];
grf_data_meta(1,:) = [];
grfdata(1,:) = [];
jointdata(1,:) = [];
joints(1,:) = [];
markerdata(1,:) = [];
neutral(1,:) = [];
velocities(1,:) = [];

% b. find subject indices to remove
q1 = quantile(cell2mat(rms_cols3to5),[0,0.75])+iqr(cell2mat(rms_cols3to5));
excl35 = find(cell2mat(rms_cols3to5)>q1(2));

% c. remove those subjects
store_alignment(excl35,:) = [];
walkRunList(excl35,:) = [];
events(excl35,:) = [];
grf_data_meta(excl35,:) = [];
grfdata(excl35,:) = [];
jointdata(excl35,:) = [];
joints(excl35,:) = [];
markerdata(excl35,:) = [];
neutral(excl35,:) = [];
velocities(excl35,:) = [];
rms_cols3to5(excl35) = [];
rms_cols6to8(excl35) = [];
D_3to5(excl35,:) = [];
D_6to8(excl35,:) = [];
R_3to5(excl35,:) = [];
R_6to8(excl35,:) = [];

q2 = quantile(cell2mat(rms_cols6to8),[0,0.75])+iqr(cell2mat(rms_cols6to8));
excl68 = find(cell2mat(rms_cols6to8)>q2(2));

store_alignment(excl68,:) = [];
walkRunList(excl68,:) = [];
events(excl68,:) = [];
grf_data_meta(excl68,:) = [];
grfdata(excl68,:) = [];
jointdata(excl68,:) = [];
joints(excl68,:) = [];
markerdata(excl68,:) = [];
neutral(excl68,:) = [];
velocities(excl68,:) = [];
rms_cols3to5(excl68) = [];
rms_cols6to8(excl68) = [];
D_3to5(excl68,:) = [];
D_6to8(excl68,:) = [];
R_3to5(excl68,:) = [];
R_6to8(excl68,:) = [];

nan35 = find(isnan(cell2mat(rms_cols3to5)));

store_alignment(nan35,:) = [];
walkRunList(nan35,:) = [];
events(nan35,:) = [];
grf_data_meta(nan35,:) = [];
grfdata(nan35,:) = [];
jointdata(nan35,:) = [];
joints(nan35,:) = [];
markerdata(nan35,:) = [];
neutral(nan35,:) = [];
velocities(nan35,:) = [];
rms_cols3to5(nan35) = [];
rms_cols6to8(nan35) = [];
D_3to5(nan35,:) = [];
D_6to8(nan35,:) = [];
R_3to5(nan35,:) = [];
R_6to8(nan35,:) = [];

nan68 = find(isnan(cell2mat(rms_cols6to8)));

store_alignment(nan68,:) = [];
walkRunList(nan68,:) = [];
events(nan68,:) = [];
grf_data_meta(nan68,:) = [];
grfdata(nan68,:) = [];
jointdata(nan68,:) = [];
joints(nan68,:) = [];
markerdata(nan68,:) = [];
neutral(nan68,:) = [];
velocities(nan68,:) = [];
rms_cols3to5(nan68) = [];
rms_cols6to8(nan68) = [];
D_3to5(nan68,:) = [];
D_6to8(nan68,:) = [];
R_3to5(nan68,:) = [];
R_6to8(nan68,:) = [];
% 151 subjects remaining

%% 3. check size of GRF
it = 1;
for i = 1:1:151
    if length(grfdata{i, 1}.R_C_mm) ~= 60000
        excl_short(it) = i;
        it=it+1;
    end  
end
store_alignment(excl_short,:) = [];
walkRunList(excl_short,:) = [];
events(excl_short,:) = [];
grf_data_meta(excl_short,:) = [];
grfdata(excl_short,:) = [];
jointdata(excl_short,:) = [];
joints(excl_short,:) = [];
markerdata(excl_short,:) = [];
neutral(excl_short,:) = [];
velocities(excl_short,:) = [];
rms_cols3to5(excl_short) = [];
rms_cols6to8(excl_short) = [];
D_3to5(excl_short,:) = [];
D_6to8(excl_short,:) = [];
R_3to5(excl_short,:) = [];
R_6to8(excl_short,:) = [];

% 138 remaining
%% 4. loop through each subject (the fun bit)
bad_ids = [3,12,17,18,35,36,40,49,67,101,103,104,109,116,119,120,122,127,131,134,135];
good_ids =[1,2,4,5,6,7,10,14,21,25,42,43,44,47,48,51,53,64,72,76,83,85,98,107,108,113,117];

store_KJC = cell(length(events),1);
store_KAMb = cell(length(events),56);
store_KAMt = cell(length(events),56);
store_FPA = cell(length(events),56);
store_vGRF = cell(length(events),56);
store_meanKAMb = cell(length(events),1);
store_meanSKAM_LR = cell(length(events),10);
store_meanSKAM_FDA = cell(length(events),10);
store_meanVGRF = cell(length(events),1);
store_meanFPA = cell(length(events),1);

for subject =1:1:length(events)%:1:40
    % a. get dims
    BW_kg = grf_data_meta{subject,3};
    h_cm = grf_data_meta{subject,2};

    % b. get KJC
        %i. find static ACS
        
        % reorder markers
    L_med_knee = [joints{subject,1}.L_med_knee(3), joints{subject,1}.L_med_knee(1), joints{subject,1}.L_med_knee(2)];
    L_lat_knee = [joints{subject,1}.L_lat_knee(3), joints{subject,1}.L_lat_knee(1), joints{subject,1}.L_lat_knee(2)];
    L_hip = [joints{subject,1}.L_hip(3), joints{subject,1}.L_hip(1), joints{subject,1}.L_hip(2)];
    L_thigh_4 = [neutral{subject,1}.L_thigh_4(3), neutral{subject,1}.L_thigh_4(1), neutral{subject,1}.L_thigh_4(2)];
    L_thigh_3 = [neutral{subject,1}.L_thigh_3(3), neutral{subject,1}.L_thigh_3(1), neutral{subject,1}.L_thigh_3(2)];
    L_thigh_2 = [neutral{subject,1}.L_thigh_2(3), neutral{subject,1}.L_thigh_2(1), neutral{subject,1}.L_thigh_2(2)];
    L_thigh_1 = [neutral{subject,1}.L_thigh_1(3), neutral{subject,1}.L_thigh_1(1), neutral{subject,1}.L_thigh_1(2)];
    
        % find T_lab_femur
    KJC_static = (L_med_knee + L_lat_knee)/2;
    ACS_o = KJC_static;
    ACS_x = L_lat_knee - L_med_knee;
    ACS_x = ACS_x/norm(ACS_x); 
    temp = L_hip - L_lat_knee;
    temp = temp/norm(temp);
    ACS_z = cross(ACS_x, temp);
    ACS_z = ACS_z/norm(ACS_z);
    ACS_y = cross(ACS_z, ACS_x);
    T_lab_femur = [[ACS_x'; 0], [ACS_y'; 0], [ACS_z'; 0], [ACS_o'; 1]];
    
    
        if ismember(subject,bad_ids)
        % use different coordinate system for tracking
            TCS_o = L_thigh_3;
            TCS_x = L_thigh_1 - L_thigh_2;
            TCS_x = TCS_x/norm(TCS_x);
            temp = L_thigh_4 - L_thigh_1;
            temp = temp/norm(temp);
            TCS_z = cross(TCS_x, temp);
            TCS_z = temp/norm(TCS_z);
            TCS_y = cross(TCS_z, TCS_x);
        
        elseif ismember(subject,good_ids)
        
        % ii. find static TCs
            TCS_o = L_thigh_4;
            TCS_x = L_thigh_1 - L_thigh_2;
            TCS_x = TCS_x/norm(TCS_x);
            temp = L_thigh_3 - L_thigh_1;
            temp = temp/norm(temp);
            TCS_z = cross(TCS_x, temp);
            TCS_z = temp/norm(TCS_z);
            TCS_y = cross(TCS_z, TCS_x);
            
        else
            TCS_o = L_thigh_1;
            TCS_x = L_thigh_1 - L_thigh_2;
            TCS_x = TCS_x/norm(TCS_x);
            temp = L_thigh_3 - L_thigh_4;
            temp = temp/norm(temp);
            TCS_z = cross(TCS_x, temp);
            TCS_z = temp/norm(TCS_z);
            TCS_y = cross(TCS_z, TCS_x);
        end
    
        % assemble transformation matrix for the TCS:
    T_lab_thigh = [[TCS_x'; 0], [TCS_y'; 0], [TCS_z'; 0], [TCS_o'; 1]];
    
        % relate TCS to ACS
    T_TCS_ACS = inv(T_lab_thigh)*T_lab_femur;   
    
        % iii. loop through each timestep, find dynamic TCS
        for j = 1:1:12000
            
            % reorder markers:
            % dynamic
            L_thigh_4 = [markerdata{subject,1}.L_thigh_4(j,3), markerdata{subject,1}.L_thigh_4(j,1), markerdata{subject,1}.L_thigh_4(j,2)];
            L_thigh_3 = [markerdata{subject,1}.L_thigh_3(j,3), markerdata{subject,1}.L_thigh_3(j,1), markerdata{subject,1}.L_thigh_3(j,2)];
            L_thigh_2 = [markerdata{subject,1}.L_thigh_2(j,3), markerdata{subject,1}.L_thigh_2(j,1), markerdata{subject,1}.L_thigh_2(j,2)];
            L_thigh_1 = [markerdata{subject,1}.L_thigh_1(j,3), markerdata{subject,1}.L_thigh_1(j,1), markerdata{subject,1}.L_thigh_1(j,2)];
            
            if ismember(subject,bad_ids)
                dTCS_o = L_thigh_3;
                dTCS_x = L_thigh_1 - L_thigh_2;
                dTCS_x = dTCS_x/norm(dTCS_x);
                temp = L_thigh_4 - L_thigh_1;
                temp = temp/norm(temp);
                dTCS_z = cross(dTCS_x, temp);
                dTCS_z = temp/norm(dTCS_z);
                dTCS_y = cross(dTCS_z, dTCS_x);
                dT_lab_thigh = [[dTCS_x'; 0], [dTCS_y'; 0], [dTCS_z'; 0], [dTCS_o'; 1]];
            elseif ismember(subject,good_ids)
                dTCS_o = L_thigh_4;
                dTCS_x = L_thigh_1 - L_thigh_2;
                dTCS_x = dTCS_x/norm(dTCS_x);
                temp = L_thigh_3 - L_thigh_1;
                temp = temp/norm(temp);
                dTCS_z = cross(dTCS_x, temp);
                dTCS_z = temp/norm(dTCS_z);
                dTCS_y = cross(dTCS_z, dTCS_x);
                dT_lab_thigh = [[dTCS_x'; 0], [dTCS_y'; 0], [dTCS_z'; 0], [dTCS_o'; 1]];
            else
                dTCS_o = L_thigh_1;
                dTCS_x = L_thigh_1 - L_thigh_2;
                dTCS_x = dTCS_x/norm(dTCS_x);
                temp = L_thigh_3 - L_thigh_4;
                temp = temp/norm(temp);
                dTCS_z = cross(dTCS_x, temp);
                dTCS_z = temp/norm(dTCS_z);
                dTCS_y = cross(dTCS_z, dTCS_x);
                dT_lab_thigh = [[dTCS_x'; 0], [dTCS_y'; 0], [dTCS_z'; 0], [dTCS_o'; 1]];
            end
            % iv. compute dynamic ACS
            T_lab_femur = dT_lab_thigh * T_TCS_ACS;
            
            % v. save origin of dACS as KJC
            KJCr(j,:) = T_lab_femur(1:3,4)';
        end
    store_KJC{subject} = KJCr;
    
    
    % c. change GRF to be consistent with markers, compute COP offset
   L_cop = (R_6to8{subject}*grfdata{subject}.L_C_mm' + repmat(D_6to8{subject},1,size(grfdata{subject}.L_C_mm,1)))'; 
   %L_cop = [L_cop(:,1) L_cop(:,3) -L_cop(:,2)]; %ML,V,AP
   L_cop = [-L_cop(:,2) L_cop(:,1) L_cop(:,3)]; %AP,ML,V
   
   L_grf = [-grfdata{subject}.L_F_N(:,2) -grfdata{subject}.L_F_N(:,1)  -grfdata{subject}.L_F_N(:,3)]; 
    

    % d. loop through each step   
    stepct = 1;
    KAMbrunsum = 0;
    KAMtrunsum = 0;
    VGRFrunsum = 0;
    FPArunsum = 0;
    LR_KAM = cell(10,length(events{subject,1})-2);
    FDA_KAM = cell(10,length(events{subject,1})-2);
    for step =2:1:length(events{subject,1})-1 % skip first and last steps
        
        % interpolate GRF
        grf_step = L_grf(events{subject,1}(step,1)*5:events{subject,1}(step,2)*5,:);
        grf_step1 = interp1(linspace(1,length(grf_step(:,1)),length(grf_step(:,1))),grf_step(:,1),linspace(1,length(grf_step(:,1)),100));
        grf_step2 = interp1(linspace(1,length(grf_step(:,2)),length(grf_step(:,2))),grf_step(:,2),linspace(1,length(grf_step(:,2)),100));
        grf_step3 = interp1(linspace(1,length(grf_step(:,3)),length(grf_step(:,3))),grf_step(:,3),linspace(1,length(grf_step(:,3)),100));
        grf_step = [grf_step1; grf_step2; grf_step3];
        
        % interpolate COP
        cop_step = L_cop(events{subject,1}(step,1)*5:events{subject,1}(step,2)*5,:);
        cop_step1 = interp1(linspace(1,length(cop_step(:,1)),length(cop_step(:,1))),cop_step(:,1),linspace(1,length(cop_step(:,1)),100));
        cop_step2 = interp1(linspace(1,length(cop_step(:,2)),length(cop_step(:,2))),cop_step(:,2),linspace(1,length(cop_step(:,2)),100));
        cop_step3 = interp1(linspace(1,length(cop_step(:,3)),length(cop_step(:,3))),cop_step(:,3),linspace(1,length(cop_step(:,3)),100));
        cop_step = [cop_step1; cop_step2; cop_step3];
        
        % interpolate KJC
        kjc_step = KJCr(events{subject,1}(step,1):events{subject,1}(step,2),:);
        kjc_step1 = interp1(linspace(1,length(kjc_step(:,1)),length(kjc_step(:,1))),kjc_step(:,1),linspace(1,length(kjc_step(:,1)),100));
        kjc_step2 = interp1(linspace(1,length(kjc_step(:,2)),length(kjc_step(:,2))),kjc_step(:,2),linspace(1,length(kjc_step(:,2)),100));
        kjc_step3 = interp1(linspace(1,length(kjc_step(:,3)),length(kjc_step(:,3))),kjc_step(:,3),linspace(1,length(kjc_step(:,3)),100));
        kjc_step = [kjc_step1; kjc_step2; kjc_step3];
        
        % compute KAM
        r = cop_step/1000 - kjc_step/1000;
        KM = cross(grf_step-cop_step/1000, r);
        KAMb = (KM(1,:)/((BW_kg*h_cm/10)))*100;
        
        % compute sKAM
        if max(KAMb)>50
            ... # skip this step because it was no good
        else
            stepct = stepct+1;
            
            % compute LR and FDA sKAM at all FPAs
            
            for fpa = 1:1:10
                % LR
                % add offsets 
                KJCnew_x = kjc_step1+(ap_pred_KJC(fpa,:)*1000);
                KJCnew_y = kjc_step2+(ml_pred_KJC(fpa,:)*1000);
                KJCnew_z = kjc_step3;
                COPnew_x = cop_step1+(ap_pred_COP(fpa,:)*1000);
                COPnew_y = cop_step2+(ml_pred_COP(fpa,:)*1000);
                COPnew_z = cop_step3;
                
                
                % calc
                r = [COPnew_x;COPnew_y;COPnew_z]/1000-[KJCnew_x; KJCnew_y;KJCnew_z]/1000;
                KM = cross(grf_step-[COPnew_x;COPnew_y;COPnew_z]/1000,r);
                KAM = (KM(1,:)/(BW_kg*(h_cm/100)*9.8))*100;
                   
                
                % store that step (average later)
                LR_KAM{fpa,step-1} = KAM;
                
                % FDA
                % add offsets
                KJCnew_x = kjc_step1+x_FDA_KJC(fpa,:);
                KJCnew_y = kjc_step2+y_FDA_KJC(fpa,:);
                KJCnew_z = kjc_step3;
                COPnew_x = cop_step1+x_FDA_COP(fpa,:);
                COPnew_y = cop_step2+y_FDA_COP(fpa,:);
                COPnew_z = cop_step3;
                
                
                % calc
                r = [COPnew_x;COPnew_y;COPnew_z]/1000-[KJCnew_x; KJCnew_y;KJCnew_z]/1000;
                KM = cross(grf_step-[COPnew_x;COPnew_y;COPnew_z]/1000,r);
                KAM = (KM(1,:)/(BW_kg*(h_cm/100)*9.8))*100;
                
                
                % store that step (average later)
                FDA_KAM{fpa,step-1} = KAM;
                
                
            end
            

        end
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check number of heel cluster
        % compute FPA
        heelr = markerdata{subject, 1}.L_foot_2(events{subject,1}(step,1):events{subject,1}(step,2),:)/1000;
        toer = markerdata{subject, 1}.L_toe(events{subject,1}(step,1):events{subject,1}(step,2),:)/1000;
        % compute FPA vector
        footVec = toer-heelr;
        FPAr = zeros(length(heelr),1);
        for j = 1:length(footVec)
            FPAr(j) = atand(footVec(j,1)/footVec(j,3));
        end
        store_FPA{subject,step} = FPAr(ceil(end/2));%mean(FPAr);
        store_KAMb{subject,step} = KAMb;
        store_vGRF{subject,step}= grf_step3;
        
        KAMbrunsum = KAMbrunsum+KAMb;
        VGRFrunsum = VGRFrunsum+grf_step3;
        FPArunsum = FPArunsum + mean(FPAr);
        
    end
    
    % average sKAM for each FPA

    
    emptyCells = cellfun('isempty', LR_KAM); 

    LR_KAM(:,all(emptyCells,1)) = [];
    FDA_KAM(:,all(emptyCells,1)) = [];
    
    
    LR_KAM = reshape(mean(reshape(cell2mat(LR_KAM), [10,100,length(LR_KAM)]), 3),[10,100]); % [10,100]
    FDA_KAM = reshape(mean(reshape(cell2mat(FDA_KAM), [10,100,length(FDA_KAM)]), 3),[10,100]); % [10,100]
    for m = 1:1:10
        store_meanSKAM_LR{subject,m} = LR_KAM(m,:);
        store_meanSKAM_FDA{subject,m} = FDA_KAM(m,:);
    end
    

    store_meanKAMb{subject} = KAMbrunsum/stepct;
    store_meanVGRF{subject} = VGRFrunsum/stepct;
    store_meanFPA{subject} = FPArunsum/stepct;
    
    disp("completed subject: "+ subject)
end




%% 5. compute first peak KAM decrease
fpKAMd_LR = zeros(138,10);
fpKAMd_FDA = zeros(138,10);
for subject = 1:1:length(events)
    
    for fpa = 1:1:10
    
        [pksb,locsb] = findpeaks(store_meanKAMb{subject,1},'MinPeakDistance',40);
        [pkst,locst] = findpeaks(store_meanSKAM_LR{subject,fpa},'MinPeakDistance',40);
        [pkst2,locst2] = findpeaks(store_meanSKAM_FDA{subject,fpa},'MinPeakDistance',40);
        
        if subject==27
            [pksb,locsb] = findpeaks(store_meanKAMb{subject,1});
            [pkst,locst] = findpeaks(store_meanSKAM_LR{subject,fpa});
            [pkst2,locst2] = findpeaks(store_meanSKAM_FDA{subject,fpa});
            %fpKAMd_LR(subject,fpa) = pksb(2) - pkst(2);
            %fpKAMd_FDA(subject,fpa) = pksb(2) - pkst2(2);
            fpKAMd_LR(subject,fpa) = pksb(2) - store_meanSKAM_LR{subject,fpa}(locsb(2));
            fpKAMd_FDA(subject,fpa) = pksb(2) - store_meanSKAM_FDA{subject,fpa}(locsb(2));
        else
            %fpKAMd_LR(subject,fpa) = pksb(1) - pkst(1);
            %fpKAMd_FDA(subject,fpa) = pksb(1) - pkst2(1);
            fpKAMd_LR(subject,fpa) = pksb(1) - store_meanSKAM_LR{subject,fpa}(locsb(1));
            fpKAMd_FDA(subject,fpa) = pksb(1) - store_meanSKAM_FDA{subject,fpa}(locsb(1));
        end
    end
end



%% 6. potentially: regression?

% features:
ft_gender = zeros(138,1);
ft_height = grf_data_meta.sub_height_cm;
ft_weight = grf_data_meta.sub_weight;
ft_speed = grf_data_meta.session_gait_speed_1;
ft_vgrf = zeros(138,1);

for sub =1:1:138
    if string(grf_data_meta.gender(sub)) == 'Female'
        ft_gender(sub) = 1;
    end
    
    ft_vgrf(sub) = max(store_meanVGRF{sub,1});
end
ft_bfpa = cell2mat(store_meanFPA);
ft_align = store_alignment;


% calculate means and stds
g_mean = mean(ft_gender);
g_std = std(ft_gender);
h_mean = mean(ft_height);
h_std = std(ft_height);
w_mean = mean(ft_weight);
w_std = std(ft_weight);
s_mean = mean(ft_speed);
s_std = std(ft_speed);
b_mean = mean(ft_bfpa);
b_std = std(ft_bfpa);
a_mean = mean(ft_align);
a_std = std(ft_align);
v_mean = mean(ft_vgrf);
v_std = std(ft_vgrf);

%% normalize:

% NR commented all so I can do LASSO below:

% ft_gender = normalize(ft_gender);
% ft_height = normalize(ft_height);
% ft_weight = normalize(ft_weight);
% ft_speed = normalize(ft_speed);
% ft_vgrf = normalize(ft_vgrf);
% ft_bfpa = normalize(ft_bfpa);
% ft_align = normalize(ft_align);


% assemble with FPAd
% label_LR = reshape(fpKAMd_LR,[1380,1]); %all 138 subs at fpa1, then all 138 at fpa2, etc...
% label_FDA = reshape(fpKAMd_FDA,[1380,1]);
% ft_fpad = reshape(repmat([1:10],138,1), [1380,1]);
% ft = [ft_height,ft_weight,ft_speed,ft_bfpa,ft_align,ft_vgrf];
% 
% all_LR = [ repmat(ft,10,1), (ft_fpad-mean(1:10))/std(1:10), label_LR];
% all_FDA = [ repmat(ft,10,1), (ft_fpad-mean(1:10))/std(1:10), label_FDA];


% put all features together with label:
%%%%%% all = [ft_gender,ft_height,ft_weight,ft_speed,ft_vgrf,ft_bfpa,ft_align,fpKAMd];
% best performing model: simple linear regression using all features,
% RMSE=0.07 %BWHT and std of decrease is 0.0792


%% exported best model:
% LR: GPR, RMSE = 0.0187

% load CMU data
%load('D:\sKAM Classifier (2021)\CMUmeta16_final.mat') % CMUmeta delete s18 
CMU_height = (CMUmeta(:,2)*100 - h_mean)/h_std;
CMU_weight = (CMUmeta(:,3) - w_mean)/w_std;
CMU_speed = (CMUmeta(:,4) - s_mean)/s_std;
CMU_bfpa = (CMUmeta(:,5) - b_mean)/b_std;
CMU_align = (CMUmeta(:,7) - a_mean)/a_std;
CMU_fpad = (CMUmeta(:,6) - mean(1:10))/std(1:10);
CMU_vgrf = (CMUmeta(:,9) - v_mean)/v_std;

% arrange CMU data
all_CMU_LR = [CMU_height, CMU_weight,CMU_speed,CMU_bfpa, CMU_align,CMU_vgrf,CMU_fpad];
    
%gt = [0.559, 0.755, 0.491, 0.077,0.4668, 0.8447, 0.416, 0.306, 0.5422, 0.392, 0.497, 1.386,0.906,0.028, 0.365, 0.705, 0.65, 0.171];
gt = [0.913, 0.755, 0.491, 0.077, 0.791, 0.7522, 0.388, 0.298, 0.8154, 0.644, 0.629, 1.386, 0.906,0.028, 0.371, 0.705, 0.678, 0.065];

%% LASSO 8.13.2021

% CALGARY:

    % 1. split 80/10/10
calg_train_idx = 1:1:108;
calg_val_idx = 109:1:123;
calg_test_idx = 124:1:138;

% features:
ft_height = grf_data_meta.sub_height_cm;
ft_weight = grf_data_meta.sub_weight;
ft_speed = grf_data_meta.session_gait_speed_1;
ft_vgrf = zeros(138,1);
for sub =1:1:138    
    ft_vgrf(sub) = max(store_meanVGRF{sub,1});
end
ft_bfpa = cell2mat(store_meanFPA);
ft_align = store_alignment;

% calculate means and stds for training set
h_mean = mean(ft_height(calg_train_idx));
h_std = std(ft_height(calg_train_idx));
w_mean = mean(ft_weight(calg_train_idx));
w_std = std(ft_weight(calg_train_idx));
s_mean = mean(ft_speed(calg_train_idx));
s_std = std(ft_speed(calg_train_idx));
b_mean = mean(ft_bfpa(calg_train_idx));
b_std = std(ft_bfpa(calg_train_idx));
a_mean = mean(ft_align(calg_train_idx));
a_std = std(ft_align(calg_train_idx));
v_mean = mean(ft_vgrf);
v_std = std(ft_vgrf(calg_train_idx));

% normalize
calg_height = (ft_height - h_mean)/h_std;
calg_weight = (ft_weight - w_mean)/w_std;
calg_speed = (ft_speed - s_mean)/s_std;
calg_vgrf = (ft_vgrf - v_mean)/v_std;
calg_bfpa = (ft_bfpa - b_mean)/b_std;
calg_align = (ft_align - a_mean)/a_std;

% assemble into train/val/test
calg_train = [calg_height(calg_train_idx), calg_weight(calg_train_idx), calg_speed(calg_train_idx), ...
    calg_bfpa(calg_train_idx), calg_align(calg_train_idx), calg_vgrf(calg_train_idx)];
calg_val = [calg_height(calg_val_idx), calg_weight(calg_val_idx), calg_speed(calg_val_idx), ...
    calg_bfpa(calg_val_idx), calg_align(calg_val_idx), calg_vgrf(calg_val_idx)];
calg_test = [calg_height(calg_test_idx), calg_weight(calg_test_idx), calg_speed(calg_test_idx), ...
    calg_bfpa(calg_test_idx), calg_align(calg_test_idx), calg_vgrf(calg_test_idx)];

% add FPAd and label
calg_train_y = reshape(fpKAMd_FDA(calg_train_idx,:), [1080,1]);
calg_val_y = reshape(fpKAMd_FDA(calg_val_idx,:), [150,1]);
calg_test_y = reshape(fpKAMd_FDA(calg_test_idx,:), [150,1]);

fpad_train = reshape(repmat([1:10],108,1),[1080,1]);
fpad_valtest = reshape(repmat([1:10],15,1),[150,1]);

calg_train_all = [ repmat(calg_train,10,1), (fpad_train-mean(1:10))/std(1:10), calg_train_y ];
calg_val_all = [ repmat(calg_val,10,1), (fpad_valtest-mean(1:10))/std(1:10), calg_val_y ];
calg_test_all = [ repmat(calg_test,10,1), (fpad_valtest-mean(1:10))/std(1:10), calg_test_y ];


% 2. tune alpha parameter for LASSO

for alpha = 0.75
    [B, FitInfo] = lasso(calg_train_all(:,[1:5,7]), calg_train_all(:,8), 'Alpha', alpha, 'CV', 5);
    idxLambda1SE = FitInfo.Index1SE; % sparse = 52, full = 1
    idxLambda1SE = 1;%FitInfo.IndexMinMSE; % full set
    coef = B(:,idxLambda1SE);
    coef0 = FitInfo.Intercept(idxLambda1SE);

    yhat = calg_val_all(:,[1:5,7])*coef + coef0;

    rmse = sqrt(mean((yhat - calg_val_all(:,8)).^2));
    
    disp(["ALPHA: " + alpha + " RMSE: " + rmse]); % because best model uses all features
    FitInfo.IndexMinMSE
end

%% 3. calculate MAE 
% LASSO
yhat_lasso = calg_test_all(:,[1:5,7])*coef + coef0;
mae_calg = mean(abs(yhat_lasso - calg_test_all(:,8)))
mas_calg = std(abs(yhat_lasso - calg_test_all(:,8)))

% FPAd linreg only:
calg_LR = fitlm(calg_train_all(:,7),calg_train_all(:,8));
yhat_lr = predict(calg_LR, calg_test_all(:,7));
mae_calg = mean(abs(yhat_lr - calg_test_all(:,8)))
mas_calg = std(abs(yhat_lr - calg_test_all(:,8)))

%% 4. make scatter plot
subplot(1,2,1)
hold on
b = scatter( calg_train_all(:,[1:5,7])*coef+coef0, calg_train_all(:,8),6);
b.MarkerFaceAlpha = 0.3;
b.MarkerFaceColor = [0.47,0.72,0.77];
b.MarkerEdgeAlpha = 0.3;

a = scatter( calg_test_all(:, [1:5,7])*coef+coef0, calg_test_all(:,8));
a.MarkerFaceAlpha = 0.6;
a.MarkerFaceColor = [0.64,0.08,0.18];
a.MarkerEdgeAlpha = 0;

xlim([0,1])
ylim([0,1])
plot([0,1],[0,1], 'k--')
xlabel('Predicted KAM Decrease [%BW*HT]')
ylabel('Actual KAM Decrease [%BW*HT]')

% CI
subplot(1,2,2)
xlim([0,3])
hold on

%training:
error = calg_train_all(:,8) - (calg_train_all(:,[1:5,7])*coef+coef0);
error = mean(reshape(error,[10,108]),1);

a1 = scatter(ones(108,1),error);
a1.MarkerFaceAlpha = 0.3;
a1.MarkerFaceColor = [0.47,0.72,0.77];
a1.MarkerEdgeAlpha = 0.3;

plot([0.75,1.25],ones(1,2)*mean(error),'k:', 'LineWidth', 2)

sem1 = std(error)/sqrt(length(error));
ts1 = tinv([0.025  0.975],length(error)-1);
ci1 = mean(error) + ts1*sem1;
plot([1,1],[ci1(1),ci1(2)],'k', 'LineWidth', 2)
plot([0.875,1.125],[ci1(1),ci1(1)],'k','LineWidth', 2)
plot([0.875,1.125],[ci1(2),ci1(2)],'k','LineWidth', 2)
ylim([-0.35, 0.35])

ci1

% test: 

error1 =  calg_test_all(:,8) - (calg_test_all(:,[1:5,7])*coef+coef0);
error1 = mean(reshape(error1,[10,15]),1);
plot([1.75,2.25],ones(1,2)*mean(error1),'k:', 'LineWidth', 2)

sem1 = std(error1)/sqrt(length(error1));
ts1 = tinv([0.025  0.975],length(error1)-1);
ci1 = mean(error1) + ts1*sem1;
plot([2,2],[ci1(1),ci1(2)],'k', 'LineWidth', 2)
plot([1.875,2.125],[ci1(1),ci1(1)],'k','LineWidth', 2)
plot([1.875,2.125],[ci1(2),ci1(2)],'k','LineWidth', 2)


ci1

a1 = scatter(2*ones(15,1),error1);
a1.MarkerFaceAlpha = 0.6;
a1.MarkerFaceColor = [0.64,0.08,0.18];
a1.MarkerEdgeAlpha = 0;

ylabel('Error [Actual - Predicted]')

%% CMU TEST ACCURACY
gt = [0.913, 0.755, 0.5,  0.791, 0.74749, 0.398, 0.298, 0.8154, 0.656, 0.629, 0.906, 0.371, 0.705, 0.659, 0.164];

%load('D:\sKAM Classifier (2021)\CMUmeta16_final.mat') % CMUmeta delete s18 
CMU_height = (CMUmeta(:,2)*100 - h_mean)/h_std;
CMU_weight = (CMUmeta(:,3) - w_mean)/w_std;
CMU_speed = (CMUmeta(:,4) - s_mean)/s_std;
CMU_bfpa = (CMUmeta(:,5) - b_mean)/b_std;
CMU_align = (CMUmeta(:,7) - a_mean)/a_std;
CMU_fpad = (CMUmeta(:,6) - mean(1:10))/std(1:10);
CMU_vgrf = (CMUmeta(:,9) - v_mean)/v_std;

% arrange CMU data
all_CMU_LR = [CMU_height, CMU_weight,CMU_speed,CMU_bfpa, CMU_align,CMU_fpad];


% calculate MAE 
% LASSO
yhat_lasso = all_CMU_LR*coef + coef0;
mae_cmu = mean(abs(yhat_lasso - gt'))
mas_cmu = std(abs(yhat_lasso - gt'))

% FPAd linreg only:
yhat_lr = predict(calg_LR, all_CMU_LR(:,end));
mae_cmu = mean(abs(yhat_lr - gt'))
mas_cmu = std(abs(yhat_lr - gt'))


% 4. make scatter plot
figure;
subplot(1,2,1)
hold on
b = scatter( calg_train_all(:,[1:5,7])*coef+coef0, calg_train_all(:,8),6);
b.MarkerFaceAlpha = 0.3;
b.MarkerFaceColor = [0.47,0.72,0.77];
b.MarkerEdgeAlpha = 0.3;

a = scatter( all_CMU_LR*coef+coef0, gt');
a.MarkerFaceAlpha = 0.6;
a.MarkerFaceColor = [0.64,0.08,0.18];
a.MarkerEdgeAlpha = 0;
xlim([0,1])
ylim([0,1])
plot([0,1],[0,1], 'k--')

xlabel('Predicted KAM Decrease [%BW*HT]')
ylabel('Actual KAM Decrease [%BW*HT]')

% CI
subplot(1,2,2)

hold on
xlim([0,3])


% training
error = calg_train_all(:,8) - (calg_train_all(:,[1:5,7])*coef+coef0);
error = mean(reshape(error,[10,108]),1);
a1 = scatter(ones(108,1),error);
a1.MarkerFaceAlpha = 0.3;
a1.MarkerFaceColor = [0.47,0.72,0.77];
a1.MarkerEdgeAlpha = 0.3;

plot([0.75,1.25],ones(1,2)*mean(error),'k:', 'LineWidth', 2)

sem1 = std(error)/sqrt(length(error));
ts1 = tinv([0.025  0.975],length(error)-1);
ci1 = mean(error) + ts1*sem1;
plot([1,1],[ci1(1),ci1(2)],'k', 'LineWidth', 2)
plot([0.875,1.125],[ci1(1),ci1(1)],'k','LineWidth', 2)
plot([0.875,1.125],[ci1(2),ci1(2)],'k','LineWidth', 2)


% test

error2 =  gt' - (all_CMU_LR*coef+coef0);
plot([1.75,2.25],ones(1,2)*mean(error2),'k:', 'LineWidth', 2)



sem1 = std(error2)/sqrt(length(error2));
ts1 = tinv([0.025  0.975],length(error2)-1);
ci1 = mean(error2) + ts1*sem1;
plot([2,2],[ci1(1),ci1(2)],'k', 'LineWidth', 2)
plot([1.875,2.125],[ci1(1),ci1(1)],'k','LineWidth', 2)
plot([1.875,2.125],[ci1(2),ci1(2)],'k','LineWidth', 2)
ylim([-0.35, 0.35])

ci1
mean(error2)

a1 = scatter(2*ones(15,1),error2);
a1.MarkerFaceAlpha = 0.6;
a1.MarkerFaceColor = [0.64,0.08,0.18];
a1.MarkerEdgeAlpha = 0;


ylabel('Error [Actual - Predicted]')




%% compute ANOVA between MSE of synth_train, synth_test, gt_test
all_mse = [error,error1,error2']; % row vector
group_labels = [ones(1,108,1), 2*ones(1,15,1), 3*ones(1,15,1)];

%ANOVA
[p, table, stats] = anova1(all_mse, group_labels)
%Tukey HSD
[c,m,h,nms] = multcompare(stats, 'ctype','hsd')




