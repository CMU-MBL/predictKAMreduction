% Validating predictive model
% nrokh 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   input: Carnegie Mellon input features
%   output: predicted vs real KAM reduction with toe-in 
%   utils: run after Calgary_main.m to compare CMU results to Calg results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. load data
load('CMUmeta.mat')
load('meanstd_params.mat')
load('coefficients.mat')
load('offset.mat')

%% 1. create model input and test:
% a. standardize:
CMU_height = (CMUmeta(:,2)*100 - h_mean)/h_std;
CMU_weight = (CMUmeta(:,3) - w_mean)/w_std;
CMU_speed = (CMUmeta(:,4) - s_mean)/s_std;
CMU_bfpa = (CMUmeta(:,5) - b_mean)/b_std;
CMU_fpad = (CMUmeta(:,6) - mean(1:10))/std(1:10);
CMU_align = (CMUmeta(:,7) - a_mean)/a_std;

% b. arrange to match coefficient format:
all_CMU = [CMU_height, CMU_weight,CMU_speed,CMU_bfpa, CMU_align,CMU_fpad];
yhat_lasso = all_CMU*coef + coef0;
y_CMU = CMUmeta(:,8);

%% 2. If comparing to Calgary results, visualize scatter plot:
figure;
subplot(1,2,1)
hold on
b = scatter( calg_train_all(:,[1:6])*coef+coef0, calg_train_all(:,7),6);
b.MarkerFaceAlpha = 0.3;
b.MarkerFaceColor = [0.47,0.72,0.77];
b.MarkerEdgeAlpha = 0.3;

a = scatter( yhat_lasso, y_CMU);
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

% test
error2 =  y_CMU - (yhat_lasso);
plot([1.75,2.25],ones(1,2)*mean(error2),'k:', 'LineWidth', 2)

sem1 = std(error2)/sqrt(length(error2));
ts1 = tinv([0.025  0.975],length(error2)-1);
ci1 = mean(error2) + ts1*sem1;
plot([2,2],[ci1(1),ci1(2)],'k', 'LineWidth', 2)
plot([1.875,2.125],[ci1(1),ci1(1)],'k','LineWidth', 2)
plot([1.875,2.125],[ci1(2),ci1(2)],'k','LineWidth', 2)
ylim([-0.35, 0.35])

a1 = scatter(2*ones(15,1),error2);
a1.MarkerFaceAlpha = 0.6;
a1.MarkerFaceColor = [0.64,0.08,0.18];
a1.MarkerEdgeAlpha = 0;

ylabel('Error [Actual - Predicted]')

%% 3. If comparing to Calgary results, compute ANOVA between MSE:

all_mse = [error,error1,error2]; 
group_labels = [ones(1,108,1), 2*ones(1,15,1), 3*ones(1,15,1)];

% ANOVA
[p, table, stats] = anova1(all_mse, group_labels)
% Tukey HSD
[c,m,h,nms] = multcompare(stats, 'ctype','hsd')