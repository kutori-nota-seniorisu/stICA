clear; clc; close all;
%% 使用时间源成分进行CBSS
datasets_num = '2';
% 拓展因子
exFactor = 20;
% nMAD
nMAD = 1.5;
% MPD(ms)
MPD = 20;
for i = 1:10
    load(['./Data/simulation/MU_time_response/TimeCompoDatasets' datasets_num '/Time_component' num2str(i) '.mat']);
    [source, ~, PT, ~, w_new] = blindDeconvPeakFinding(MU_noisy_conv, exFactor, nMAD, MPD*2);
    decompo_pulses{i} = PT;
end

%% 绘制iEMG的pulse
data = importdata('./Data/iEMG/24-06-21/iEMG_S1_M1_level1_trial1_24-06-21_UUS.eaf');
muNum = max(data.data(:, 2));
iPulses = {};
for mu = 1:muNum
    iPulses{mu} = round(data.data(find(data.data(:, 2) == mu), 1)' * 1000);
end
plotDecomps(iPulses, [], 1000, 0, 0, []);
% xlim([3, 13])
% plotDecomps(decompoPulseAll, [], 1000, 0, 0, []);
%% 绘制ICData的pulse，截取trigger后的30s数据，转化为对应于1000Hz的pulse
load('./Data/experiment/ICdata/R10/R10.mat');
pulses = struct2cell(Data{2, 2});
for iii = 1:length(pulses)
    % 需要对原始数据转置一下，变成行的形式，调用plotDecomps时才不会绘制出错
    tmp = pulses{iii};
    tmp = tmp(tmp>=12397 & tmp<=73661);
    tmp = tmp - 12397;
    tmp = tmp/2048*1000;
    pulses_new{iii} = tmp';
end
plotDecomps(pulses_new, [], 1000, 0, 0, []);
% xlim([12397,73661]/2048);
yticks(1:length(pulses_new))
%% 绘制Data{1,2}每条通道的波形
for i = 1:13
    figure;
    for j = 1:10
        subplot(5,2,j);
        plot(Data{1,2}((i-1)*10+j,:));
        title(num2str((i-1)*10+j));
    end
end
%%
d1 = 8;
d2 = 4;
ref = 9;
figure;
subplot(5,1,1);
plot(-T(:, d1));
title('BPM')
subplot(5,1,2);
plot(-T_norm(:, d1));
hold on;
scatter(decompo_pulses{d1},zeros(length(decompo_pulses{d1})), 1000, 'black', '|');
title('BPM normalization')
subplot(5,1,3);
plot(-decompoSourceFirstAll(:, d2));
% hold on;
% scatter(decompoPulseAll{d2}, zeros(length(decompoPulseAll{d2})), 1000, "black", '|');
title('USCBSS First Step')
subplot(5,1,4);
plot(decompoSourceAll(:, d2));
hold on;
scatter(decompoPulseAll{d2}, zeros(length(decompoPulseAll{d2})), 1000, "black", '|');
title('USCBSS Second Step')
subplot(5,1,5);
plot(Xt(:, ref));
hold on;
scatter(ipulses{ref}, zeros(length(ipulses{ref})), 1000, "black", '|');
title('Ref Source')
set(gcf,'unit','normalized','position',[0.05,0.1,0.9,0.6]);

%% 将USConvBSS的结果与BPM进行对比
figure;
subplot(2,1,1);
hold on;
for i=1:10
    % load(['Results/result' num2str(i) '_Update.mat']);
    load(['Results/datasets' num2str(i) '_result.mat']);
    boxplot(matchresult_time.time, 'Positions', i-0.2, 'Colors', 'r');
    hold on;
    timeMean(i) = mean(matchresult_time.time);
    numMU(i) = size(matchresult_time, 1);
end
plot((1:10)-0.2, timeMean, 'r-*');
hold on;

for i=1:10
    % load(['Results/result' num2str(i) '_Update.mat']);
    load(['Results/datasets' num2str(i) '_resultnew.mat']);
    boxplot(matchresult_time.time, 'Positions', i, 'Colors', 'b');
    hold on;
    timeMean2(i) = mean(matchresult_time.time);
    numMU2(i) = size(matchresult_time, 1);
end
plot((1:10), timeMean2, 'b-*');
hold on;

load('Results/compo12_NC_NMFm0.9_BPM_10sets.mat');
for i = 1:10
    % space_mean_BPM(i) = mean(matchresult_final_all{i}.space);
    time_mean_BPM(i) = mean(matchresult_final_all{i}.time);
    No_BPM(i) = size(matchresult_final_all{i}, 1);
    % time_No_NMFm(i) = size(matchresult_time_all{i}, 1);
    % space_No_NMFm(i) = size(matchresult_space_all{i}, 1);
end
match_BPM = matchresult_final_all;
for i = 1:10
    boxplot(match_BPM{i}.time, 'Positions', i+0.2, 'Colors', 'g');
    hold on;
end
plot((1:10)+0.2,time_mean_BPM,'g-*');
ylabel('RoA (Time)');
xticks(1:10)
xticklabels(1:10);
legend('USCBSS', 'USCBSS2', 'BPM');

subplot(2,1,2);
plot(1:10, numMU, 'r-*');
hold on;
plot(1:10, numMU2, 'b-*');
hold on;
plot(1:10, No_BPM, 'g-*');
ylabel('No. MU');
legend('USCBSS', 'USCBSS2', 'BPM');

sgtitle('USCBSS vs BPM')
set(gcf,'unit','normalized','position',[0.05,0.1,0.9,0.6]);

%% 尝试对ICData的结果进行BPM提取脉冲串
fsampu = 1000;
[B, A] = butter(4, [5, 25]/fsampu*2);
T_filter = filtfilt(B, A, decompoSourceFirstAll);
% z-score标准化
T_mean = mean(T_filter);
T_std = std(T_filter);
T_norm = (T_filter-T_mean)./T_std;
for iii = 1:size(T_norm, 2)
    figure;
    subplot(2,1,1);
    plot(decompoSourceFirstAll(:, iii));
    subplot(2,1,2);
    plot(T_norm(:, iii));
end
decompo_pulses = {};
for i = 1:size(T_norm, 2)
    [~, locs] = findpeaks(T_norm(:,i),'MinPeakDistance',50,'MinPeakHeight',0.5);
    decompo_pulses{end+1} = locs';
end
plotDecomps(decompo_pulses, [], 1000, 0, 0, []);

%%
close all;
figure;
for iii = 1:length(pulses_new)
    subplot(length(pulses_new),1,iii);
    % figure;
    plot(decompoSourceAll(:, 4));
    % plot(T_norm(:, 1))
    hold on;
    scatter(pulses_new{iii}, zeros(length(pulses_new{iii})), 3000, 'red', '|');
    ylim([-1, 1])
    xlim([10000, 15000])
end
%% 计算RoA
matchresult_time_raw = [];
for i = 1:length(decompo_pulses)
    for j = 1:length(pulses_new)
        [Array1, Array2] = meshgrid(decompo_pulses{i}, pulses_new{j});
        diff_values = Array1 - Array2;
        valid_elements = diff_values <= 50 & diff_values >= -50;
        count = sum(valid_elements(:));
        r = count/(length(decompo_pulses{i})+length(pulses_new{j})-count);
        if r > 1
            r = 1;
        end
        spike_ROA_matrix(i,j) = r;
        matchresult_time_raw(end+1,:) = [i, j, r];
    end
end

%%
% 检查异常值
outlier_ratio = zeros(1, size(TVI_data, 2));
for i = 1:size(TVI_data, 2)
    Q1 = quantile(TVI_data(:, i), 0.25);
    Q3 = quantile(TVI_data(:, i), 0.75);
    IQR = Q3 - Q1;
    outliers = TVI_data(:, i) < (Q1 - 1.5*IQR) | TVI_data(:, i) > (Q3 + 1.5*IQR);
    outlier_ratio(i) = sum(outliers) / length(TVI_data(:, i));
end

if any(outlier_ratio > 0.05)
    fprintf('检测到异常值，建议使用按列归一化或Robust归一化\n');
end