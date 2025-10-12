% ipulse和6PTM模型卷积，作为仿真数据
% 模型ref: Variability of successive contractions subtracted from unfused tetanus of fast and slow motor units
clc; clear; close all;
addpath('F:/EEEMG/stICA/Func');

fsampu = 2000; % 采样率

datasets_num = '1';
time_datasets = ['TimeCompoDatasets' datasets_num];% 时间分量数据集所在的文件夹
load(['F:/EEEMG/stICA/Data/simulation/MU_time_response/' time_datasets '/ipulses.mat']);
% plotDecomps(ipulses,[],fsampu,0,0,[]);

muNum = length(ipulses);

for mu = 5%1:muNum
    % 脉冲串长度
    L = length(ipulses{mu});

    % 6PTM模型，单位ms
    % 考虑收缩曲线的变异性，对于每一个放电时刻都要生成一条收缩曲线
    Tlead = 0;
    Thc = unifrnd(20,25,L,1);
    Tc = unifrnd(50,75,L,1);
    Thr = unifrnd(100,130,L,1);
    Ttw = unifrnd(300,350,L,1);
    Dmax = unifrnd(40,70,L,1);

    MU_conv_all = [];

    for i = 1:L
        t = 0:Ttw(i);

        e = exp(1);
        c1 = log(2)*Tc(i) / (Thc(i) - Tc(i) + Tc(i)*log(Tc(i)/Thc(i)));
        c2 = log(2)*Tc(i) / (Thr(i) - Tc(i) + Tc(i)*log(Tc(i)/Thr(i)));

        P1 = (t - Tlead) / Tc(i);
        P2 = 1 + exp(2*e*(P1-1));
        P3 = (t - 0.5*(Ttw(i)+Thr(i))) / (Ttw(i)-Thr(i));

        f1 = Dmax(i) * (P1.^c1 .* exp(c1-c1*P1) + (P2-1) .* P1.^c2 .* exp(c2-c2*P1));
        f2 = P2 + P2.*exp(4*e*P3);

        F = f1 ./ f2;
        
        pulse = zeros(1,2*fsampu);
        pulse(ipulses{mu}(i)) = 1;

        tmp = conv(pulse, F);
        MU_conv_all(i,:) = tmp(1:4000);
    end

    MU_conv = sum(MU_conv_all);
    % 差分，得到速度曲线
    MU_conv_diff = diff(MU_conv);
    MU_conv_diff(end+1) = 0;

    % 添加高斯噪声
    snr = 20; % 信噪比为~dB
    MU_noisy_conv = awgn(MU_conv_diff, snr, 'measured');

    % figure;
    % subplot(3,1,1);
    % plot(MU_conv);
    % title(['MU' num2str(mu) ' twitch']);
    % subplot(3,1,2);
    % plot(MU_conv_diff);
    % title(['MU' num2str(mu) ' twitch(diff)']);
    % subplot(3,1,3);
    % plot(MU_noisy_conv);
    % title(['MU' num2str(mu) ' twitch(diff,noisy)']);

    save(['F:/EEEMG/stICA/Data/simulation/MU_time_response/' time_datasets '/Time_component' num2str(mu) '.mat'],'MU_noisy_conv')
end
