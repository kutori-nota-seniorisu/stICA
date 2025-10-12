% 比较十个数据集的分解结果
clear; clc; close all;
addpath('./Func');
%% 绘制两种方式分解结果
weight = 0.9;
load(['F:/EEEMG/stICA/Results/compo12_NC_NMFm' num2str(weight) '_10sets.mat']);
for i = 1:10
    space_mean_NMFm(i) = mean(matchresult_final_all{i}.space);
    time_mean_NMFm(i) = mean(matchresult_final_all{i}.time);
    No_NMFm(i) = size(matchresult_final_all{i}, 1);
    % time_No_NMFm(i) = size(matchresult_time_all{i}, 1);
    % space_No_NMFm(i) = size(matchresult_space_all{i}, 1);
end
match_NMFm = matchresult_final_all;

load(['F:/EEEMG/stICA/Results/compo12_NC_NMFmXt' num2str(weight) '_10sets.mat']);
for i = 1:10
    space_mean_NMFmXt(i) = mean(matchresult_final_all{i}.space);
    time_mean_NMFmXt(i) = mean(matchresult_final_all{i}.time);
    No_NMFmXt(i) = size(matchresult_final_all{i}, 1);
    % time_No_NMFmXt(i) = size(matchresult_time_all{i}, 1);
    % space_No_NMFmXt(i) = size(matchresult_space_all{i}, 1);
end
match_NMFmXt = matchresult_final_all;

figure;
subplot(3,1,1);
for i = 1:10
    boxplot(match_NMFm{i}.space, 'Positions', i-0.1, 'Colors', 'r');
    hold on;
    boxplot(match_NMFmXt{i}.space, 'Positions', i+0.1, 'Colors', 'g');
    hold on;
end
plot((1:10)-0.1,space_mean_NMFm,'r-*');
hold on;
plot((1:10)+0.1,space_mean_NMFmXt,'g-*');
legend('NMFm','NMFmXt');
% xlabel('dataset');
ylabel('Corr2 (Space)');
% title('Space Fit');
tick = num2cell(1:10);
xticks(1:10)
xticklabels(tick)

subplot(3,1,2);
for i = 1:10
    boxplot(match_NMFm{i}.time, 'Positions', i-0.1, 'Colors', 'r');
    hold on;
    boxplot(match_NMFmXt{i}.time, 'Positions', i+0.1, 'Colors', 'g');
    hold on;
end
plot((1:10)-0.1,time_mean_NMFm,'r-*');
hold on;
plot((1:10)+0.1,time_mean_NMFmXt,'g-*');
legend('NMFm','NMFmXt');
% xlabel('dataset');
ylabel('RoA (Time)');
% title('Time Fit');
xticks(1:10)
xticklabels(tick)

subplot(3,1,3);
plot(1:10, No_NMFm, 'm-*');
hold on;
plot(1:10, No_NMFmXt, 'k-*');
ylabel('No. MU');
legend('NMFm', 'NMFmXt');
sgtitle(['\alpha = ' num2str(weight)])
set(gcf,'unit','normalized','position',[0.05,0.1,0.9,0.6]);


%% NMFm在不同权重下的情况
figure;
weight = 0.1;
load(['F:/EEEMG/stICA/Results/compo12_NC_NMFm' num2str(weight) '_10sets.mat']);
for i = 1:10
    space_mean_NMFm(i) = mean(matchresult_final_all{i}.space);
    time_mean_NMFm(i) = mean(matchresult_final_all{i}.time);
    No_NMFm(i) = size(matchresult_final_all{i}, 1);
    % time_No_NMFm(i) = size(matchresult_time_all{i}, 1);
    % space_No_NMFm(i) = size(matchresult_space_all{i}, 1);
end
subplot(3,1,1);
plot(1:10, No_NMFm, 'm-*');
ylabel('No. MU');
hold on;
subplot(3,1,2);
plot(1:10, space_mean_NMFm, 'm-*');
ylabel('Corr2 mean');
hold on;
subplot(3,1,3);
plot(1:10, time_mean_NMFm, 'm-*');
ylabel('RoA mean')
hold on;

weight = 0.5;
load(['F:/EEEMG/stICA/Results/compo12_NC_NMFm' num2str(weight) '_10sets.mat']);
for i = 1:10
    space_mean_NMFm(i) = mean(matchresult_final_all{i}.space);
    time_mean_NMFm(i) = mean(matchresult_final_all{i}.time);
    No_NMFm(i) = size(matchresult_final_all{i}, 1);
    % time_No_NMFm(i) = size(matchresult_time_all{i}, 1);
    % space_No_NMFm(i) = size(matchresult_space_all{i}, 1);
end
subplot(3,1,1);
plot(1:10, No_NMFm, 'b-*');
hold on;
subplot(3,1,2);
plot(1:10, space_mean_NMFm, 'b-*');
hold on;
subplot(3,1,3);
plot(1:10, time_mean_NMFm, 'b-*');
hold on;

weight = 0.9;
load(['F:/EEEMG/stICA/Results/compo12_NC_NMFm' num2str(weight) '_10sets.mat']);
for i = 1:10
    space_mean_NMFm(i) = mean(matchresult_final_all{i}.space);
    time_mean_NMFm(i) = mean(matchresult_final_all{i}.time);
    No_NMFm(i) = size(matchresult_final_all{i}, 1);
    % time_No_NMFm(i) = size(matchresult_time_all{i}, 1);
    % space_No_NMFm(i) = size(matchresult_space_all{i}, 1);
end
subplot(3,1,1);
plot(1:10, No_NMFm, 'g-*');
hold on;
subplot(3,1,2);
plot(1:10, space_mean_NMFm, 'g-*');
hold on;
subplot(3,1,3);
plot(1:10, time_mean_NMFm, 'g-*');
hold on;
legend('0.1', '0.5', '0.9');
sgtitle('不同权重情况下的结果比较 拟合sech^2')
set(gcf,'unit','normalized','position',[0.05,0.1,0.9,0.6]);

%% 比较CBSS与BPM
load('Results/compo12_NC_NMFm0.9_BPM_10sets.mat');
for i = 1:10
    space_mean_BPM(i) = mean(matchresult_final_all{i}.space);
    time_mean_BPM(i) = mean(matchresult_final_all{i}.time);
    No_BPM(i) = size(matchresult_final_all{i}, 1);
    % time_No_NMFm(i) = size(matchresult_time_all{i}, 1);
    % space_No_NMFm(i) = size(matchresult_space_all{i}, 1);
end
match_BPM = matchresult_final_all;

load('Results/compo12_NC_NMFm0.9_CBSS_10sets.mat');
for i = 1:10
    space_mean_CBSS(i) = mean(matchresult_final_all{i}.space);
    time_mean_CBSS(i) = mean(matchresult_final_all{i}.time);
    No_CBSS(i) = size(matchresult_final_all{i}, 1);
    % time_No_NMFmXt(i) = size(matchresult_time_all{i}, 1);
    % space_No_NMFmXt(i) = size(matchresult_space_all{i}, 1);
end
match_CBSS = matchresult_final_all;

figure;
subplot(3,1,1);
for i = 1:10
    boxplot(match_BPM{i}.space, 'Positions', i-0.1, 'Colors', 'r');
    hold on;
    boxplot(match_CBSS{i}.space, 'Positions', i+0.1, 'Colors', 'g');
    hold on;
end
plot((1:10)-0.1,space_mean_BPM,'r-*');
hold on;
plot((1:10)+0.1,space_mean_CBSS,'g-*');
legend('BPM','CBSS');
% xlabel('dataset');
ylabel('Corr2 (Space)');
% title('Space Fit');
tick = num2cell(1:10);
xticks(1:10)
xticklabels(tick)

subplot(3,1,2);
for i = 1:10
    boxplot(match_BPM{i}.time, 'Positions', i-0.1, 'Colors', 'r');
    hold on;
    boxplot(match_CBSS{i}.time, 'Positions', i+0.1, 'Colors', 'g');
    hold on;
end
plot((1:10)-0.1,time_mean_BPM,'r-*');
hold on;
plot((1:10)+0.1,time_mean_CBSS,'g-*');
legend('BPM','CBSS');
% xlabel('dataset');
ylabel('RoA (Time)');
% title('Time Fit');
xticks(1:10)
xticklabels(tick)

subplot(3,1,3);
plot(1:10, No_BPM, 'm-*');
hold on;
plot(1:10, No_CBSS, 'k-*');
ylabel('No. MU');
legend('BPM', 'CBSS');
sgtitle('BPM vs CBSS')
set(gcf,'unit','normalized','position',[0.05,0.1,0.9,0.6]);

%% 统计iEMG的特征
for set_i = 1:10
    load(['F:/EEEMG/stICA/Data/simulation/MU_time_response/TimeCompoDatasets' num2str(set_i) '/ipulses.mat']);
    for i = 1:length(ipulses)
        ipulse = ipulses{i}(ipulses{i}<=4000);
        counts(set_i,i) = length(ipulse);
        interval = diff(ipulse);
        % figure;
        % plot(interval,'g-s');
        % title(['MU #' num2str(i)])
        CoV(set_i,i) = std(interval)/mean(interval);
    end

    plotDecomps(ipulses,[],2000,0,0,[]);
    xlim([0,2]);
end