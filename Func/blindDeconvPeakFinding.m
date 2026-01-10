function [source, PT, CoV, w_new] = blindDeconvPeakFinding(S, fs, R, nMAD, MPD, mode)
%
% 这是一种采用卷积盲源分离获取spike train的方法，
% 输入为估计源，可以是ICA后得到的单个时间成分，也可以是对USS进行CBSS后的粗略估计源
%
% Inputs:
%   S - estimated twitch
%   fs - sample rate
%   R - extend factor
%   nMAD - the number of mean absolute distances
%   MPD - the minimum peak distance (ms)
%   mode - the mode of contrast function
% Outputs:
%   source - the estimated source
%   PT - the times of pulses
%   CoV - the Coefficient of Variation
%   w_new - the separation vector
%
% ref: Estimating the neural spike train from an unfused tetanic signal of 
% low‑threshold motor units using convolutive blind source separation

% 迭代容差
Tolx = 1e-4;

% 拓展信号
eS = extend(S, R);
eS = eS(:, 1:size(S, 2));

% 减去均值
eS = eS - mean(eS, 2);
% 协方差矩阵特征值分解
[V, D] = eig(cov(eS'));
[d, idx] = sort(diag(D), 'descend');
V = V(:, idx);
D = diag(d);%+1e-3*max(d));
% 白化矩阵WM，采用PCA白化格式
WM = sqrt(D)\V';
% 白化后的数据
Z = WM*eS;
% 迭代更新
w_new = randn(R, 1);
% 迭代计数器
iterCount = 0;
while true
    w_old = w_new;
    % 固定点迭代
    if mode == 1
        w_new = mean(Z.*log(cosh(w_old'*Z)), 2) - mean(tanh(w_old'*Z)).*w_old;
    elseif mode == 2
        % w_new = mean(Z.*(1/6*(w_old'*Z).^3), 2) - mean(1/2*(w_old'*Z).^2).*w_old;
        w_new = mean(Z.*(1/2.*(w_old'*Z).^2), 2) - mean(w_old'*Z).*w_old;
        % w_new = mean(Z.*(w_old'*Z).^2, 2) - mean(2*w_old'*Z).*w_old;
    end
    % 归一化处理
    w_new = w_new / norm(w_new);
    % 记录迭代次数
    iterCount = iterCount + 1;
    if abs(w_new'*w_old - 1) < Tolx || iterCount >= 100
        disp(['迭代完成，本次迭代' num2str(iterCount) '次']);
        break;
    end
end

% 提取放电时刻
source = w_new' * Z;
[~, PT] = findpeaks(source, 'MinPeakHeight', nMAD*mad(source), 'MinPeakDistance', MPD*(fs/1000));

% 计算变异系数
CoV = std(diff(PT)) / mean(diff(PT));