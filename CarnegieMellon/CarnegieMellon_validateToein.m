% Validating gait patterns
% nrokh 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   input: Carnegie Mellon dataset from SimTK (Subject and Meta)
%   output: synthetic KAM visualization, toe-in accuracy 
%   utils: none; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPath = 'D:\sKAM Classifier (2021)\CMUdata_toSimTK\';

%% 0. load data
% make empty store for KAMs:
bKAMs = zeros(1,100);
tKAMs = zeros(1,100);
fdaKAMs = zeros(1,100);

% make empty store for COP and KJCs
tCOPs = zeros(2,100); %ap, ml
tKJCs = zeros(2,100);
fdaCOPs = zeros(2,100); %ap, ml
fdaKJCs = zeros(2,100);

sublist = 1:1:15;
for i =1:1:length(sublist)
    sub = string(sublist(i));
    
    if sublist(i) < 10
        load(dataPath + "s0" + sub + "\pp_s0" + sub)
    else
        load(dataPath + "s" + sub + "\pp_s" + sub)
    end
    bKAMs = bKAMs+bKAM;
    tKAMs = tKAMs+tKAM;
    fdaKAMs = KAM_FDA+fdaKAMs;
    
    tCOPs = tCOPs + tCOP;
    tKJCs = tKJCs + tKJC;
    
    fdaCOPs = fdaCOPs + FDACOP;
    fdaKJCs = fdaKJCs + FDAKJC;
    
    
    rmse1(i) =  sqrt(mean((tCOP(1,[10:95])-FDACOP(1,[10:95])).^2)); %ap
    rmse2(i) =  sqrt(mean((tCOP(2,[10:95])-FDACOP(2,[10:95])).^2));
    rmse3(i) =  sqrt(mean((tKJC(1,:)-FDAKJC(1,:)).^2)); %ap
    rmse4(i) =  sqrt(mean((tKJC(2,:)-FDAKJC(2,:)).^2));
    
    store_COPap_real(i,:) = tCOP(1,[10:95]);
    store_COPap_pred(i,:) = FDACOP(1,[10:95]);
    store_COPml_real(i,:) = tCOP(2,[10:95]);
    store_COPml_pred(i,:) = FDACOP(2,[10:95]);
    
    store_KAM_real(i,:) = tKAM;
    store_KAM_fda(i,:) = KAM_FDA;
    store_KAM_base(i,:) = bKAM;

end

%% 1. visualize accuracy of toe-in patterns

figure;
sgtitle('CMU toe-in accuracy')
subplot(1,4,1)
boxplot(rmse3*1000);
hold on
plot(ones(15,1), rmse3*1000,  '*')
title(["RMSE KJC AP", "mean = " + mean(rmse3*1000), "std = " + std(rmse3*1000)])
subplot(1,4,2)
boxplot(rmse4*1000);
hold on
plot(ones(15,1), rmse4*1000,  '*')
title(["RMSE KJC ML", "mean = " + mean(rmse4*1000), "std = " + std(rmse4*1000)])
subplot(1,4,3)
boxplot(rmse1*1000);
hold on
plot(ones(15,1), rmse1*1000,  '*')
title(["RMSE COP AP", "mean = " + mean(rmse1*1000), "std = " + std(rmse1*1000)])
subplot(1,4,4)
boxplot(rmse2*1000);
hold on
plot(ones(15,1), rmse2*1000,  '*')
title(["RMSE COP ML", "mean = " + mean(rmse2*1000), "std = " + std(rmse2*1000)])

%% 2. visualize synthetic vs real toe-in KAM
figure()

hold on

sem1 = std(store_KAM_fda)/sqrt(length(store_KAM_fda));
ts1 = tinv([0.025  0.975],length(store_KAM_fda)-1);
ci1 = mean(store_KAM_fda) + ts1'.*sem1;
x = 1:length(ci1);
x2 = [x, fliplr(x)];
inbetween = [ci1(1,:), fliplr(ci1(2,:))];
color = [0.64,0.08,0.18];
a2 = fill(x2, inbetween, color);
a2.FaceAlpha = 0.3;
a2.EdgeAlpha = 0;

sem2 = std(store_KAM_real)/sqrt(length(store_KAM_real));
ci2 = mean(store_KAM_real) + ts1'.*sem2;
inbetween = [ci2(1,:), fliplr(ci2(2,:))];
color = [0.47,0.72,0.77];
a2 = fill(x2, inbetween, color);
a2.FaceAlpha = 0.3;
a2.EdgeAlpha = 0;

sem3 = std(store_KAM_base)/sqrt(length(store_KAM_base));
ci3 = mean(store_KAM_base) + ts1'.*sem3;
inbetween = [ci3(1,:), fliplr(ci3(2,:))];
color = [0.80,0.80,0.80];
a2 = fill(x2, inbetween, color);
a2.FaceAlpha = 0.8;
a2.EdgeAlpha = 0;

p1 = plot(mean(store_KAM_base), 'LineWidth', 3);
p2 = plot(mean(store_KAM_real), 'LineWidth', 3);
p3 = plot(mean(store_KAM_fda), 'LineWidth', 3);

p2.Color = [0.47,0.72,0.77];
p1.Color = [0 0 0];
p3.Color = [0.64,0.08,0.18];

ylim([-2,4])