% 使用BPM的效果比move median的效果好一点
clear; clc; close all;
addpath('./Func')
%% 导入数据
load('./Data/experiment/24-06-21/UUS-iEMG/S1M1L1T1_C12_NMFm_1.mat');
% load('F:/EEEMG/stICA_simple/Results_UU/S1M1L1T1_compo12_NMFm_skew_stICA_1.mat');
%% spike train extraction
T=(T(2:end,:)-mean(T(2:end,:),1))./std(T(2:end,:),1); %~ 幅值标准化，第一行的值过大？
% move median
T_smooth = smoothdata(T, "movmedian", 8);
decompo_pulses = {};
for i = 1:size(T_smooth,2)
    [~, locs] = findpeaks(T_smooth(:,i),'MinPeakDistance',95); % 最后一个参数视情况调节，对匹配准确度影响较大
    decompo_pulses{end+1} = locs';
end

fsampu = 1000;
data = importdata('./Data/EMG/24-06-21/iEMG_S1_M1_level1_trial1_24-06-21_UUS.eaf');
muNum = max(data.data(:, 2));
pulses_2s = {};
iPulses = {};
for mu = 1:muNum
    tmp = round(data.data(find(data.data(:,2)==mu),1)'*fsampu);
    iPulses{mu} = tmp;
    tmp(tmp<=5000) = [];
    tmp(tmp>7000) = [];
    tmp = tmp-5000;
    pulses_2s{mu} = tmp;
end
plotDecomps(iPulses, [], fsampu, 0, 0, []);
plotDecomps(pulses_2s, [], fsampu, 0, 0, []);

% spike匹配ROA
for i = 1:length(decompo_pulses)
    for j = 1:length(pulses_2s)
        [Array1, Array2] = meshgrid(decompo_pulses{i}, pulses_2s{j});
        diff_values = Array1 - Array2;
        % 相较于iEMG，在-10ms~30ms之内都算匹配上
        valid_elements = diff_values <= 30 & diff_values >= -10;
        count = sum(valid_elements(:));
        % spike_ROA_matrix(i,j) = count/size(pulses_2s{j},2);
        spike_ROA_matrix(i,j) = count/(length(decompo_pulses{i})+length(pulses_2s{j})-count);
    end
end

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

figure;
subplot(3,1,1);
plot(T(:,2)); title('原始分量')
subplot(3,1,2);
plot(T_norm1(:,2)); title('move median & z-score')
subplot(3,1,3);
plot(T_norm2(:,2)); title('band pass filter & z-score')

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

%% 导入iEMG
data = importdata('F:/EEEMG/stICA_simple/Data/iEMG/24-06-21/iEMG_S1_M1_level1_trial1_24-06-21_UUS.eaf');
muNum = max(data.data(:,2)); % MU的个数
pulses_2s = {};
for mu = 1:muNum
    ipulse = round(data.data(find(data.data(:,2)==mu),1)'*fsampu);
    ipulse(find(ipulse<=3000)) = [];
    ipulse(find(ipulse>13000)) = [];
    ipulse = ipulse - 3000;
    pulses_2s{mu} = ipulse;
end

%% 绘制spike train
plotDecomps(pulses_2s,[],fsampu,0,0,[]);
title('iPulses');
plotDecomps(decompo1_pulses,[],fsampu,0,0,[]);
title('decompose pulses 1');
plotDecomps(decompo2_pulses,[],fsampu,0,0,[]);
title('decompose pulses 2');
plotDecomps(decompo3_pulses,[],fsampu,0,0,[]);
title('decompose pulses 3');

%% 脉冲串匹配
[match_result1, match_result_raw1] = PulseMatch(decompo1_pulses, pulses_2s, 0.3, fsampu);
[match_result2, match_result_raw2] = PulseMatch(decompo2_pulses, pulses_2s, 0.3, fsampu);
[match_result3, match_result_raw3] = PulseMatch(decompo3_pulses, pulses_2s, 0.3, fsampu);
match_result1
match_result2
match_result3
