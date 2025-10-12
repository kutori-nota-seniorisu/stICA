% 使用BPM的效果比move median的效果好一点
clear; clc; close all;
addpath('./Func')
%% 导入数据
load('F:/EEEMG/stICA_simple/Results_UU/S1M1L1T1_compo12_NMFm_skew_stICA_1.mat');
%% spike train extraction
% move median
T_smooth = smoothdata(T, "movmedian", 8);

% z-score标准化
T_mean1 = mean(T_smooth);
T_std1 = std(T_smooth);
T_norm1 = (T_smooth-T_mean1)./T_std1;

% 5-25Hz带通滤波
fsampu = 1000;
[B,A] = butter(4,[5,25]/fsampu*2);
T_filter = filtfilt(B,A,T);

% z-score标准化
T_mean2 = mean(T_filter);
T_std2 = std(T_filter);
T_norm2 = (T_filter-T_mean2)./T_std2;

% HWM哈尔小波变换
scale = 24;
for i = 1:size(T,2)
    tmp = cwt(T(:,i),scale,'haar');
    T_HWM(:,i) = tmp';
end

% z-score标准化
T_mean3 = mean(T_HWM);
T_std3 = std(T_HWM);
T_norm3 = (T_HWM-T_mean3)./T_std3;

figure;
subplot(4,1,1);
plot(T(:,2)); title('原始分量')
subplot(4,1,2);
plot(T_norm1(:,2)); title('move median & z-score')
subplot(4,1,3);
plot(T_norm2(:,2)); title('band pass filter & z-score')
subplot(4,1,4);
plot(T_norm3(:,2)); title('HWM & z-score')

% 采样率1000Hz，故1ms=1样本点。设置MinPeakDistance时需要注意ms与样本点的转换
decompo1_pulses = {};
for i = 1:size(T_norm1,2)
    [~, locs] = findpeaks(T_norm1(:,i),'MinPeakDistance',25,'MinPeakHeight',0.1);
    decompo1_pulses{end+1} = locs';
end

decompo2_pulses = {};
for i = 1:size(T_norm2,2)
    [~, locs] = findpeaks(T_norm2(:,i),'MinPeakDistance',25,'MinPeakHeight',0.1);
    decompo2_pulses{end+1} = locs';
end

decompo3_pulses = {};
for i = 1:size(T_norm3,2)
    [~, locs] = findpeaks(T_norm3(:,i),'MinPeakDistance',25,'MinPeakHeight',0.94);
    decompo3_pulses{end+1} = locs';
end

%% 导入iEMG
data = importdata('F:/EEEMG/stICA_simple/Data/iEMG/24-06-21/iEMG_S1_M1_level1_trial1_24-06-21_UUS.eaf');
muNum = max(data.data(:,2)); % MU的个数
iPulses = {};
for mu = 1:muNum
    ipulse = round(data.data(find(data.data(:,2)==mu),1)'*fsampu);
    ipulse(find(ipulse<=3000)) = [];
    ipulse(find(ipulse>13000)) = [];
    ipulse = ipulse - 3000;
    iPulses{mu} = ipulse;
end

%% 绘制spike train
plotDecomps(iPulses,[],fsampu,0,0,[]);
title('iPulses');
plotDecomps(decompo1_pulses,[],fsampu,0,0,[]);
title('decompose pulses 1');
plotDecomps(decompo2_pulses,[],fsampu,0,0,[]);
title('decompose pulses 2');
plotDecomps(decompo3_pulses,[],fsampu,0,0,[]);
title('decompose pulses 3');

%% 脉冲串匹配
[match_result1, match_result_raw1] = PulseMatch(decompo1_pulses, iPulses, 0.3, fsampu);
[match_result2, match_result_raw2] = PulseMatch(decompo2_pulses, iPulses, 0.3, fsampu);
[match_result3, match_result_raw3] = PulseMatch(decompo3_pulses, iPulses, 0.3, fsampu);
match_result1
match_result2
match_result3
