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
xlim([3, 13])
% plotDecomps(decompoPulseAll, [], 1000, 0, 0, []);
%% 绘制ICData的pulse，截取trigger后的30s数据，转化为对应于1000Hz的pulse
% 时间成分滤波
% fsampu = 1000;
% [B,A] = butter(4,[5,25]/fsampu*2);
% T_filter = filtfilt(B,A,T);
T_filter = T;
% z-score标准化
T_mean = mean(T_filter);
T_std = std(T_filter);
T_norm = (T_filter-T_mean)./T_std;

% figure;
% subplot(2,1,1);
% plot(T(:,1));
% subplot(2,1,2);
% plot(T_norm(:,1));

% 把处理后的时间成分画出来
figure;
for i = 1:size(T,2)
    subplot(2,6,i);
    plot(-T_norm(:,i));
    title(['Time #' num2str(i)]);
end
sgtitle('标准化 时间成分');
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

decompo_pulses = {};
% 提取脉冲串
for i = 1:size(T_norm,2)
    [~, locs] = findpeaks(-T_norm(:,i),'MinPeakDistance',100,'MinPeakHeight',1);
    decompo_pulses{i} = locs';
end
plotDecomps(decompo_pulses, [], fsampu, 0, 0, []);

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
plotDecomps(pulses_new, [], fsampu, 0, 0, []);
% xlim([12397,73661]/2048);
yticks(1:length(pulses_new))
%%
matchresult_time_raw = [];
% 计算RoA
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
% plotDecomps(pulses_new, [], 1000, 0, 0, []);
plotDecomps({decompo_pulses{1},pulses_new{4}}, [], 1000, 0, 0, []);
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
        valid_elements = diff_values <= 15*2 & diff_values >= 0;
        count = sum(valid_elements(:));
        r = count/(length(decompo_pulses{i})+length(pulses_new{j})-count);
        if r > 1
            r = 1;
        end
        spike_ROA_matrix(i,j) = r;
        matchresult_time_raw(end+1,:) = [i, j, r];
    end
end


%% 统计ipulses的CoV
for i = 1:length(ipulses)
    CoV(i) = std(diff(ipulses{i})) / mean(diff(ipulses{i}));
end

%%
ddd = 4;
rrr = 9;
plotDecomps({ipulses{rrr}, decompoPulseAll{ddd}}, [], 2000, 0, 0, []);