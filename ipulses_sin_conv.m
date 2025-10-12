% ipulse和正弦函数卷积，作为仿真数据
% 使用的ipulses是来自stICA_simple/Data/simulation/MU_time_response_sin/TimeCompoDatasets/ipulses.mat
clear; clc; close all;
addpath('F:/EEEMG/stICA/Func');

fsampu = 2000; % 采样率

datasets_num = '10';
time_datasets = ['TimeCompoDatasets' datasets_num];% 时间分量数据集所在的文件夹
load(['F:/EEEMG/stICA/Data/simulation/MU_time_response/' time_datasets '/ipulses.mat']);
plotDecomps(ipulses,[],fsampu,0,0,[]);

%% 生成一个单周期的正弦函数
x = linspace(0, 200, 200);
y = sin(2*pi*x/200);
figure
plot(x, y);

%% 绘制机械单位MU的响应MUAPT
figure;
for i = 1:length(ipulses)
    pulse = zeros(2*fsampu,1);
    pulse(ipulses{i}) = 1;
    MU_conv = conv(pulse, y);

    % 添加高斯噪声
    snr = 20; % 信噪比为~dB
    MU_noisy_conv = awgn(MU_conv(1:end-length(x)+1), snr, 'measured');

    subplot(2,5,i);
    plot(MU_noisy_conv);
    title(['Time #' num2str(i)])
    save(['F:/EEEMG/stICA/Data/simulation/MU_time_response/' time_datasets '/Time_component' num2str(i) '.mat'],'MU_noisy_conv')
end
clear MU_noisy_conv
