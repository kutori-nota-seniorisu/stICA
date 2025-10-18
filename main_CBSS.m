%% 对stICA后的单个成分进行卷积盲源分离
clear; clc; close all;
addpath('./Func');
datasets_num = '1';
load(['Data/simulation/datasets' datasets_num '/UV_compo12_NC_NMFm0.9.mat']);
% 拓展因子
exFactor = 20;
% nMAD
nMAD = 4;
% MPD(ms)
MPD = 20;

decompo_pulses = {};
for i = 1:12
    % 是否要在T的前面加上负号？我认为不需要，迭代结果不受符号的影响
    [source, PT, ~, w_new] = blindDeconvPeakFinding(Xt(:, i)', exFactor, nMAD, MPD*2, 1);
    decompo_pulses{i} = PT;
    T_norm(:, i) = source';
    W(:, i) = w_new;
end
