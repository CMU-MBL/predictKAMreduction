% Checking synthetic toe-in Stanford
% nrokh 2021


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   input: preprocessed Stanford dataset, FDA fit gait patterns
%   output: synthetic KAM visualization, toe-in accuracy 
%   utils: none; run Stanford_preproc.m first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. Compute accuracy of toe-in gait patterns and synthetic KAM using LOOCV


for i =1:1:12
    fpa = target_FPA(i);
    if fpa>size(KJCx_offset,1) % set to max if outside FDA range
        fpa = size(KJCx_offset,1); 
    end
    
    % a. load toe-in patterns learned from all but left-out subject
    load("KJCy_FDA" + i) 
    load("KJCx_FDA" + i)
    load("COPy_FDA" + i)
    load("COPx_FDA" + i)
       

    % b. FDA pred:
        % i. baseline and toe-in GRF and Pelvis are same:
    pred_grf = [baseline_GRF_x(i,:); baseline_GRF_y(i,:); baseline_GRF_z(i,:)]; 
    pred_pelvis = [baseline_pelvis_x(i,:); baseline_pelvis_y(i,:); baseline_pelvis_z(i,:)];
    
        % ii. add KJC and COP offsets to baseline trajectory:
    pred_kjc = [KJCx_offset(fpa,:) + baseline_KJC_x(i,:); KJCy_offset(fpa,:) + baseline_KJC_y(i,:);  baseline_KJC_z(i,:)];
    pred_cop = [COPx_offset(fpa,:) + baseline_COP_x(i,:); COPy_offset(fpa,:) + baseline_COP_y(i,:); baseline_COP_z(i,:)];
        % iii. relate all back to global by adding back pelvis:
    pred_kjc_g = pred_kjc + pred_pelvis;
    pred_cop_g = pred_cop + pred_pelvis;
    
       
    % c. calculate RMSE:
    real_kjc_g = [trial_KJC_x(i,:); trial_KJC_y(i,:); trial_KJC_z(i,:)] ;
    real_cop_g = [trial_COP_x(i,:); trial_COP_y(i,:); trial_COP_z(i,:)] ;
    
    rmseKJCap(i,1) = sqrt(mean((real_kjc_g(1,[10:95])-pred_kjc(1,[10:95])).^2));
    rmseKJCml(i,1) = sqrt(mean((real_kjc_g(2,[10:95])-pred_kjc(2,[10:95])).^2));
    rmseCOPap(i,1) = sqrt(mean((real_cop_g(1,[10:95])-pred_cop(1,[10:95])).^2));
    rmseCOPml(i,1) = sqrt(mean((real_cop_g(2,[10:95])-pred_cop(2,[10:95])).^2));
    
     % d. calculate synthetic KAM:
    r = pred_cop_g/1000 - pred_kjc_g/1000;
    KM = cross(r, pred_grf);
    KAM_fda(i,:) = (-KM(1,:)/(subjectHeightsWeights(i,2)*subjectHeightsWeights(i,3)/100*9.8))*100;

    % e. compare to ground-truth baseline and toe-in KAM:
        % i. toe-in
    real_kjc = [trial_KJC_x(i,:); trial_KJC_y(i,:); trial_KJC_z(i,:)];
    real_cop = [trial_COP_x(i,:); trial_COP_y(i,:); trial_COP_z(i,:)];
    real_grf = [trial_GRF_x(i,:); trial_GRF_y(i,:); trial_GRF_z(i,:)];
    real_pelvis = [trial_pelvis_x(i,:); trial_pelvis_y(i,:); trial_pelvis_z(i,:)];
    real_kjc_g = real_kjc + real_pelvis;
    real_cop_g = real_cop + real_pelvis;
    
    r = real_cop_g/1000 - real_kjc_g/1000;
    KM2 = cross(r, real_grf);
    KAM_toein(i,:) = (-KM2(1,:)/(subjectHeightsWeights(i,2)*subjectHeightsWeights(i,3)/100*9.8))*100;
    
        % ii. baseline
    base_kjc = [baseline_KJC_x(i,:); baseline_KJC_y(i,:); baseline_KJC_z(i,:)];
    base_cop = [baseline_COP_x(i,:); baseline_COP_y(i,:); baseline_COP_z(i,:)];
    base_grf = [baseline_GRF_x(i,:); baseline_GRF_y(i,:); baseline_GRF_z(i,:)];
    base_pelvis = [baseline_pelvis_x(i,:); baseline_pelvis_y(i,:); baseline_pelvis_z(i,:)];
    base_kjc_g = base_kjc + base_pelvis;
    base_cop_g = base_cop + base_pelvis;
     
    r = base_cop_g/1000 - base_kjc_g/1000;
    KM3 = cross(r, base_grf);
    KAM_base(i,:) = (-KM3(1,:)/(subjectHeightsWeights(i,2)*subjectHeightsWeights(i,3)/100*9.8))*100;

end

%% 2. Visualization:

figure;
subplot(1,4,1)
sgtitle('LOOCV accuracy')
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

figure()
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

legend('synthetic toein', 'real toe-in', 'baseline')