%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing synthetic KAM from postprocessed data
%    nrokh 2021.06.15
%
% input: "PP_ + {sub_num} + toein/walk + LEFT/RIGHT" processed .mat file
%
% output: single averaged synthetic KAM
%
% dependencies: FDA and LR offsets 
%
% purpose: generating and comparing synthetic KAM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% 0.0  META AND FILE LOADING


dataPath = 'C:\Users\nrokh\Box\CMU_MBL\Data\SageMotion_CTSleeve\02_Post-Processed Data\';
sub = '22';                 % de-id subject number

% subject parameters
wt_kg = 63.5;
ht_m = 1.59;

foot = 'R';                 % = 'L' or 'R'

% load files
left_walk = load(dataPath + "s" + sub + "\PP_" + sub + "walk_LEFT.mat");
right_walk = load(dataPath + "s" + sub + "\PP_" + sub + "walk_RIGHT.mat");

left_toein = load(dataPath + "s" + sub + "\PP_" + sub + "toein_LEFT.mat");
right_toein = load(dataPath + "s" + sub + "\PP_" + sub + "toein_RIGHT.mat");


right_walk.mean_RFPA
%%  0.1 LOADING OFFSETS

% load LR offsets
load('D:\Shull_FPA\FINAL\RESULTS_11272020\Shull_offsets.mat')

% load FDA offsets
% load('D:\sKAM Classifier (2021)\FDA\offsets\y_FDA_KJC.mat')
% load('D:\sKAM Classifier (2021)\FDA\offsets\x_FDA_KJC.mat')
% load('D:\sKAM Classifier (2021)\FDA\offsets\y_FDA_COP.mat')
% load('D:\sKAM Classifier (2021)\FDA\offsets\x_FDA_COP.mat')

% using raw splines instead of regressing on splines
load('D:\sKAM Classifier (2021)\FDA\KJCx_FDA.mat')
load('D:\sKAM Classifier (2021)\FDA\KJCy_FDA.mat')
load('D:\sKAM Classifier (2021)\FDA\COPx_FDA.mat')
load('D:\sKAM Classifier (2021)\FDA\COPy_FDA.mat')
x_FDA_KJC = KJCx_offset;
y_FDA_KJC = KJCy_offset;
x_FDA_COP = COPx_offset;
y_FDA_COP = COPy_offset;

%% 1. FIND IDS OF TOE-IN WRT BASELINE STEPS
if foot == 'R'
    % 5,9,13,21
    id = find(right_toein.store_RFPA(:,50)<right_walk.mean_RFPA)%-2 & right_toein.store_RFPA(:,50)>right_walk.mean_RFPA-8) % at least 2 deg toe-in
    
else
    id = find(left_toein.store_LFPA(:,50)<left_walk.mean_LFPA-2 & left_toein.store_LFPA(:,50)>left_walk.mean_LFPA-5); % at least 2 deg toe-in
    
end

disp("Number of steps: "+ length(id))

if strcmp(sub,'20')
    id(6) = [];
end
%% 2. COMPUTE FPA CHANGE

if foot == 'R'
    FPAdiff = round(right_walk.mean_RFPA - mean(right_toein.store_RFPA(id,50)));
else
    FPAdiff = round(left_walk.mean_LFPA - mean(left_toein.store_LFPA(id,50)));
end

if FPAdiff > 10
    disp('Toe-in angle outside of synthetic KAM range :( ')
    FPAdiff = 10;
end


%% 3. ADD OFFSETS

if foot == 'R'
    % linreg
    new_COP_ap = mean(right_walk.store_RCOP_ap) + ap_pred_COP(FPAdiff,:);
    new_COP_ml = mean(right_walk.store_RCOP_ml) - ml_pred_COP(FPAdiff,:);
    new_COP_v = mean(right_walk.store_RCOP_v);
    
    new_KJC_ap = mean(right_walk.store_RKJC_ap) + ap_pred_KJC(FPAdiff,:);
    new_KJC_ml = mean(right_walk.store_RKJC_ml) - ml_pred_KJC(FPAdiff,:);
    new_KJC_v = mean(right_walk.store_RKJC_v);
    
    
    % FDA
    new_COP_ap2 = mean(right_walk.store_RCOP_ap) + x_FDA_COP(FPAdiff,:)/1000;
    new_COP_ml2 = mean(right_walk.store_RCOP_ml) - y_FDA_COP(FPAdiff,:)/1000;
    new_COP_v2 = mean(right_walk.store_RCOP_v);
    
    new_KJC_ap2 = mean(right_walk.store_RKJC_ap) + x_FDA_KJC(FPAdiff,:)/1000;
    new_KJC_ml2 = mean(right_walk.store_RKJC_ml) - y_FDA_KJC(FPAdiff,:)/1000;
    new_KJC_v2 = mean(right_walk.store_RKJC_v);
else
    % linreg
    new_COP_ap = mean(left_walk.store_LCOP_ap) + ap_pred_COP(FPAdiff,:);
    new_COP_ml = mean(left_walk.store_LCOP_ml) + ml_pred_COP(FPAdiff,:);
    new_COP_v = mean(left_walk.store_LCOP_v);
    
    new_KJC_ap = mean(left_walk.store_LKJC_ap) + ap_pred_KJC(FPAdiff,:);
    new_KJC_ml = mean(left_walk.store_LKJC_ml) + ml_pred_KJC(FPAdiff,:);
    new_KJC_v = mean(left_walk.store_LKJC_v);
    
    % FDA
    new_COP_ap2 = mean(left_walk.store_LCOP_ap) + x_FDA_COP(FPAdiff,:)/1000;
    new_COP_ml2 = mean(left_walk.store_LCOP_ml) + y_FDA_COP(FPAdiff,:)/1000;
    new_COP_v2 = mean(left_walk.store_LCOP_v);
    
    new_KJC_ap2 = mean(left_walk.store_LKJC_ap) + x_FDA_KJC(FPAdiff,:)/1000;
    new_KJC_ml2 = mean(left_walk.store_LKJC_ml) + y_FDA_KJC(FPAdiff,:)/1000;
    new_KJC_v2 = mean(left_walk.store_LKJC_v);

end

%% 4. COMPUTE BASELINE KAM, REAL TOE-IN KAM, LR KAM, FDA KAM

if foot == 'R'
    %linreg
    GRF = [mean(right_walk.store_RGRF_ap); mean(right_walk.store_RGRF_ml); mean(right_walk.store_RGRF_v)];
    COP = [new_COP_ap; new_COP_ml; new_COP_v];
    KJC = [new_KJC_ap; new_KJC_ml; new_KJC_v];
    %
    r = COP - KJC;
    KM_LR = cross(r, GRF);
    KAM_LR = (KM_LR(1,:)/(wt_kg*ht_m*9.8))*100;
    
    %fda
    GRF2 = [mean(right_walk.store_RGRF_ap); mean(right_walk.store_RGRF_ml); mean(right_walk.store_RGRF_v)];
    COP2 = [new_COP_ap2; new_COP_ml2; new_COP_v2];
    KJC2 = [new_KJC_ap2; new_KJC_ml2; new_KJC_v2];
    %
    r = COP2 - KJC2;
    KM_FDA = cross(r, GRF2);
    KAM_FDA = (KM_FDA(1,:)/(wt_kg*ht_m*9.8))*100;
else
    %linreg
    GRF = [mean(left_walk.store_LGRF_ap); mean(left_walk.store_LGRF_ml); mean(left_walk.store_LGRF_v)];
    COP = [new_COP_ap; new_COP_ml; new_COP_v];
    KJC = [new_KJC_ap; new_KJC_ml; new_KJC_v];
    %
    r = COP - KJC;
    KM_LR = cross(r, GRF);
    KAM_LR = (KM_LR(1,:)/(wt_kg*ht_m*9.8))*100;
    
    %fda
    GRF2 = [mean(left_walk.store_LGRF_ap); mean(left_walk.store_LGRF_ml); mean(left_walk.store_LGRF_v)];
    COP2 = [new_COP_ap2; new_COP_ml2; new_COP_v2];
    KJC2 = [new_KJC_ap2; new_KJC_ml2; new_KJC_v2];
    %
    r = COP2 - KJC2;
    KM_FDA = cross(r, GRF2);
    KAM_FDA = (KM_FDA(1,:)/(wt_kg*ht_m*9.8))*100;
end


%% 5. PLOT EVERYTHING

if foot == 'R'
    tKAM = mean(right_toein.store_RKAM(id,:));
    bKAM = mean(right_walk.store_RKAM);
    
    figure;
    title("Right KAM at " + FPAdiff  +" degrees toe-in")
    hold on;
    plot(tKAM);
    plot(bKAM);
    plot(KAM_LR);
    plot(KAM_FDA);
    ylim([-2,4])
    legend('real toein', 'baseline', 'sKAM linreg', 'sKAM FDA')
else
    tKAM = mean(left_toein.store_LKAM(id,:));
    bKAM = mean(left_walk.store_LKAM);
    
    figure;
    title("Left KAM at " + FPAdiff  +" degrees toe-in")
    hold on;
    plot(tKAM);
    plot(bKAM);
    plot(-KAM_LR);
    plot(-KAM_FDA);
    legend('real toein', 'baseline', 'sKAM linreg', 'sKAM FDA')
end



%% 6. DISP ACCURACY
FP_R_b = max(mean(right_walk.store_RKAM(:,10:40)));
FP_R_t = max(mean(right_toein.store_RKAM(id,10:40)));
FP_R_s = max(KAM_LR(10:40));
FP_R_s2 = max(KAM_FDA(10:40));

real_RFPdecrease = FP_R_b - FP_R_t
synthLR_RFPdecrease = FP_R_b - FP_R_s
synthFDA_RFPdecrease = FP_R_b - FP_R_s2

acc_LR = synthLR_RFPdecrease - real_RFPdecrease
acc_FDA = synthFDA_RFPdecrease - real_RFPdecrease


%% 7. SAVE KAMS

tCOP = [mean(right_toein.store_RCOP_ap(id,:)); mean(right_toein.store_RCOP_ml(id,:))]...
    - [mean(right_toein.store_pelvis_ap(id,:)); mean(right_toein.store_pelvis_ml(id,:))];
FDACOP = [new_COP_ap2; new_COP_ml2]...
    - [mean(right_walk.store_pelvis_ap); mean(right_walk.store_pelvis_ml)];
tKJC = [mean(right_toein.store_RKJC_ap(id,:)); mean(right_toein.store_RKJC_ml(id,:))]...
    - [mean(right_toein.store_pelvis_ap(id,:)); mean(right_toein.store_pelvis_ml(id,:))];
FDAKJC = [new_KJC_ap2; new_KJC_ml2]...
    - [mean(right_walk.store_pelvis_ap); mean(right_walk.store_pelvis_ml)];
save(dataPath + "s" + sub + "\meanKAMs_NOLR" + sub, 'tKAM','bKAM','KAM_LR','KAM_FDA','FPAdiff','tCOP','FDACOP','tKJC','FDAKJC')

