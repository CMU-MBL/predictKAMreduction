% Synthesizing toe-in gait
% nrokh 2021

clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   input: Calgary Running Clinic dataset from SimTK, toe-in FDA patterns  
%   output: synthetic toe-in KAM
%   utils: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. load data

disp('Loading Calgary data...')
% Calgary Running Clinic Data from SimTK
load('velocities.mat')
load('store_alignment.mat')
load('neutral.mat')
load('markerdata.mat')
load('joints.mat')
load('jointdata.mat')
load('grfdata.mat')
load('grf_data_meta.mat')
load('events.mat')
load('calibration.mat')
% toe-in FDA patterns from Stanford_synthToein or Github/Models
load('y_FDA_KJC.mat')
load('x_FDA_KJC.mat')
load('y_FDA_COP.mat')
load('x_FDA_COP.mat')
disp('Loaded')

%% 1. loop through each subject
% NOTE: three different marker cluster numberings were used; these were identified manually
markerset1 = [3,12,17,18,35,36,40,49,67,101,103,104,109,116,119,120,122,127,131,134,135];
markerset2 =[1,2,4,5,6,7,10,14,21,25,42,43,44,47,48,51,53,64,72,76,83,85,98,107,108,113,117];

% initialize parameters
store_KJC = cell(length(events),1);
store_KAMb = cell(length(events),56);
store_KAMt = cell(length(events),56);
store_FPA = cell(length(events),56);
store_vGRF = cell(length(events),56);
store_meanKAMb = cell(length(events),1);
store_meanSKAM_LR = cell(length(events),10); %%%%%%%%%
store_meanSKAM_FDA = cell(length(events),10);
store_meanVGRF = cell(length(events),1);
store_meanFPA = cell(length(events),1);

for subject = 1:1:length(events) 
    % a. get subject metadata
    BW_kg = grf_data_meta{subject,3};
    h_cm = grf_data_meta{subject,2};
    
    % b. compute KJC
        % i. reorder markers to match lab conventions
    L_med_knee = [joints{subject,1}.L_med_knee(3), joints{subject,1}.L_med_knee(1), joints{subject,1}.L_med_knee(2)];
    L_lat_knee = [joints{subject,1}.L_lat_knee(3), joints{subject,1}.L_lat_knee(1), joints{subject,1}.L_lat_knee(2)];
    L_hip = [joints{subject,1}.L_hip(3), joints{subject,1}.L_hip(1), joints{subject,1}.L_hip(2)];
    L_thigh_4 = [neutral{subject,1}.L_thigh_4(3), neutral{subject,1}.L_thigh_4(1), neutral{subject,1}.L_thigh_4(2)];
    L_thigh_3 = [neutral{subject,1}.L_thigh_3(3), neutral{subject,1}.L_thigh_3(1), neutral{subject,1}.L_thigh_3(2)];
    L_thigh_2 = [neutral{subject,1}.L_thigh_2(3), neutral{subject,1}.L_thigh_2(1), neutral{subject,1}.L_thigh_2(2)];
    L_thigh_1 = [neutral{subject,1}.L_thigh_1(3), neutral{subject,1}.L_thigh_1(1), neutral{subject,1}.L_thigh_1(2)];
    
        % ii. find T_lab_femur
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
    
        % iii. find static TCS depending on marker numbering set
    if ismember(subject,markerset1)
            TCS_o = L_thigh_3;
            TCS_x = L_thigh_1 - L_thigh_2;
            TCS_x = TCS_x/norm(TCS_x);
            temp = L_thigh_4 - L_thigh_1;
            temp = temp/norm(temp);
            TCS_z = cross(TCS_x, temp);
            TCS_z = temp/norm(TCS_z);
            TCS_y = cross(TCS_z, TCS_x);
        
        elseif ismember(subject,markerset2)
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
     
        % iv. assemble transformation matrix for the TCS:
    T_lab_thigh = [[TCS_x'; 0], [TCS_y'; 0], [TCS_z'; 0], [TCS_o'; 1]];
    
        % v. relate TCS to ACS
    T_TCS_ACS = inv(T_lab_thigh)*T_lab_femur; 
    
        % vi.  loop through each timestep, find dynamic TCS
        for j = 1:1:12000
            
            % reorder markers:
            % dynamic
            L_thigh_4 = [markerdata{subject,1}.L_thigh_4(j,3), markerdata{subject,1}.L_thigh_4(j,1), markerdata{subject,1}.L_thigh_4(j,2)];
            L_thigh_3 = [markerdata{subject,1}.L_thigh_3(j,3), markerdata{subject,1}.L_thigh_3(j,1), markerdata{subject,1}.L_thigh_3(j,2)];
            L_thigh_2 = [markerdata{subject,1}.L_thigh_2(j,3), markerdata{subject,1}.L_thigh_2(j,1), markerdata{subject,1}.L_thigh_2(j,2)];
            L_thigh_1 = [markerdata{subject,1}.L_thigh_1(j,3), markerdata{subject,1}.L_thigh_1(j,1), markerdata{subject,1}.L_thigh_1(j,2)];
            
            if ismember(subject,markerset1)
                dTCS_o = L_thigh_3;
                dTCS_x = L_thigh_1 - L_thigh_2;
                dTCS_x = dTCS_x/norm(dTCS_x);
                temp = L_thigh_4 - L_thigh_1;
                temp = temp/norm(temp);
                dTCS_z = cross(dTCS_x, temp);
                dTCS_z = temp/norm(dTCS_z);
                dTCS_y = cross(dTCS_z, dTCS_x);
                dT_lab_thigh = [[dTCS_x'; 0], [dTCS_y'; 0], [dTCS_z'; 0], [dTCS_o'; 1]];
            elseif ismember(subject,markerset2)
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
            % vii. compute dynamic ACS
            T_lab_femur = dT_lab_thigh * T_TCS_ACS;
            
            % viii. save origin of dACS as KJC
            KJCr(j,:) = T_lab_femur(1:3,4)';
        end
    store_KJC{subject} = KJCr;
    
    % c. make GRF convention consistent with markers, transform COP
   L_cop = (R_6to8{subject}*grfdata{subject}.L_C_mm' + repmat(D_6to8{subject},1,size(grfdata{subject}.L_C_mm,1)))'; 
   L_cop = [-L_cop(:,2) L_cop(:,1) L_cop(:,3)]; %AP,ML,V
   L_grf = [-grfdata{subject}.L_F_N(:,2) -grfdata{subject}.L_F_N(:,1)  -grfdata{subject}.L_F_N(:,3)]; 
    
   % d. loop through each step:
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
        
        % compute synthetic KAM
        if max(KAMb)>50
            ... # skip this step because it was no good
        else
            stepct = stepct+1;
            
            % compute FDA synthetic KAM at all FPAs
            
            for fpa = 1:1:10
                
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
    % e. compute FPA
        heelr = markerdata{subject, 1}.L_foot_2(events{subject,1}(step,1):events{subject,1}(step,2),:)/1000;
        toer = markerdata{subject, 1}.L_toe(events{subject,1}(step,1):events{subject,1}(step,2),:)/1000;

        footVec = toer-heelr;
        FPAr = zeros(length(heelr),1);
        for j = 1:length(footVec)
            FPAr(j) = atand(footVec(j,1)/footVec(j,3));
        end
        store_FPA{subject,step} = FPAr(ceil(end/2));
        store_KAMb{subject,step} = KAMb;
        store_vGRF{subject,step}= grf_step3;
        
        KAMbrunsum = KAMbrunsum+KAMb;
        VGRFrunsum = VGRFrunsum+grf_step3;
        FPArunsum = FPArunsum + mean(FPAr);
    end
    
    
    
    % f. average synthetic KAM for each FPA
    emptyCells = cellfun('isempty', FDA_KAM); 

    FDA_KAM(:,all(emptyCells,1)) = [];
    FDA_KAM = reshape(mean(reshape(cell2mat(FDA_KAM), [10,100,length(FDA_KAM)]), 3),[10,100]); % [10,100]
    
    for m = 1:1:10
        store_meanSKAM_FDA{subject,m} = FDA_KAM(m,:);
    end

    store_meanKAMb{subject} = KAMbrunsum/stepct;
    store_meanVGRF{subject} = VGRFrunsum/stepct;
    store_meanFPA{subject} = FPArunsum/stepct;
    
    disp("completed subject: "+ subject)
    
end

%% 2. Compute first peak KAM decrease (synthetic) and form dataset
% a. compute first peak KAM
fpKAMd_FDA = zeros(138,10);
for subject = 1:1:length(events)    
    for fpa = 1:1:10    
        [pksb,locsb] = findpeaks(store_meanKAMb{subject,1},'MinPeakDistance',40);
        [pkst2,locst2] = findpeaks(store_meanSKAM_FDA{subject,fpa},'MinPeakDistance',40);
        
        if subject==27 % this subject had a large data artefact
            [pksb,locsb] = findpeaks(store_meanKAMb{subject,1});
            [pkst2,locst2] = findpeaks(store_meanSKAM_FDA{subject,fpa});
            fpKAMd_FDA(subject,fpa) = pksb(2) - store_meanSKAM_FDA{subject,fpa}(locsb(2));

        else
            fpKAMd_FDA(subject,fpa) = pksb(1) - store_meanSKAM_FDA{subject,fpa}(locsb(1));
        end
    end
end

% b. form features:
ft_height = grf_data_meta.sub_height_cm;
ft_weight = grf_data_meta.sub_weight;
ft_speed = grf_data_meta.session_gait_speed_1;
ft_bfpa = cell2mat(store_meanFPA);
ft_align = store_alignment;

% c. split by 80/10/10
calg_train_idx = 1:1:108;
calg_val_idx = 109:1:123;
calg_test_idx = 124:1:138;

% d. calculate means and stds for training set
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


% e. standardize all features with training set values
calg_height = (ft_height - h_mean)/h_std;
calg_weight = (ft_weight - w_mean)/w_std;
calg_speed = (ft_speed - s_mean)/s_std;
calg_bfpa = (ft_bfpa - b_mean)/b_std;
calg_align = (ft_align - a_mean)/a_std;

% f. assemble into train/val/test
calg_train = [calg_height(calg_train_idx), calg_weight(calg_train_idx), calg_speed(calg_train_idx), ...
    calg_bfpa(calg_train_idx), calg_align(calg_train_idx)];
calg_val = [calg_height(calg_val_idx), calg_weight(calg_val_idx), calg_speed(calg_val_idx), ...
    calg_bfpa(calg_val_idx), calg_align(calg_val_idx)];
calg_test = [calg_height(calg_test_idx), calg_weight(calg_test_idx), calg_speed(calg_test_idx), ...
    calg_bfpa(calg_test_idx), calg_align(calg_test_idx)];

% g. add toe-in FPA as feature and label
calg_train_y = reshape(fpKAMd_FDA(calg_train_idx,:), [1080,1]);
calg_val_y = reshape(fpKAMd_FDA(calg_val_idx,:), [150,1]);
calg_test_y = reshape(fpKAMd_FDA(calg_test_idx,:), [150,1]);

fpad_train = reshape(repmat([1:10],108,1),[1080,1]);
fpad_valtest = reshape(repmat([1:10],15,1),[150,1]);

calg_train_all = [ repmat(calg_train,10,1), (fpad_train-mean(1:10))/std(1:10), calg_train_y ];
calg_val_all = [ repmat(calg_val,10,1), (fpad_valtest-mean(1:10))/std(1:10), calg_val_y ];
calg_test_all = [ repmat(calg_test,10,1), (fpad_valtest-mean(1:10))/std(1:10), calg_test_y ];

%% 3. Train LASSO and visualize accuracy:

% a. best fit model:
alpha = 0.7;
[B, FitInfo] = lasso(calg_train_all(:,[1:6]), calg_train_all(:,7), 'Alpha', alpha, 'CV', 5);
coef = B(:,1);
coef0 = FitInfo.Intercept(1);
yhat = calg_val_all(:,[1:6])*coef + coef0;
rmse = sqrt(mean((yhat - calg_val_all(:,7)).^2));
disp(["ALPHA: " + alpha + " RMSE: " + rmse]); 

% b. visualize scatter plot:
subplot(1,2,1)
hold on
b = scatter( calg_train_all(:,[1:6])*coef+coef0, calg_train_all(:,7),6);
b.MarkerFaceAlpha = 0.3;
b.MarkerFaceColor = [0.47,0.72,0.77];
b.MarkerEdgeAlpha = 0.3;

a = scatter( calg_test_all(:, [1:6])*coef+coef0, calg_test_all(:,7));
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
error = calg_train_all(:,7) - (calg_train_all(:,[1:6])*coef+coef0);
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


% test: 
error1 =  calg_test_all(:,7) - (calg_test_all(:,[1:6])*coef+coef0);
error1 = mean(reshape(error1,[10,15]),1);
plot([1.75,2.25],ones(1,2)*mean(error1),'k:', 'LineWidth', 2)

sem1 = std(error1)/sqrt(length(error1));
ts1 = tinv([0.025  0.975],length(error1)-1);
ci1 = mean(error1) + ts1*sem1;
plot([2,2],[ci1(1),ci1(2)],'k', 'LineWidth', 2)
plot([1.875,2.125],[ci1(1),ci1(1)],'k','LineWidth', 2)
plot([1.875,2.125],[ci1(2),ci1(2)],'k','LineWidth', 2)


a1 = scatter(2*ones(15,1),error1);
a1.MarkerFaceAlpha = 0.6;
a1.MarkerFaceColor = [0.64,0.08,0.18];
a1.MarkerEdgeAlpha = 0;

ylabel('Error [Actual - Predicted]')



