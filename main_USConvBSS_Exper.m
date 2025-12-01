% 对超声速度图像进行卷积盲源分离，使用了仿真数据集，采用并行计算
clear; clc; close all;
addpath('./Func');
% 拓展因子
exFactor = 10;
% 迭代容差
Tolx = 1e-4;
% 潜在成分个数
numCompo = 25;
% 超声采样率
fsampu = 2000;

%% 开启并行池
if isempty(gcp('nocreate'))
    parpool;
end

for motion = 1:2
    for trial = 1:2

%% TVI数据预处理
disp('导入数据');
% 导入TVI数据
tviFile = ['./Data/experiment/25-07-04/TVIData_15000_S_wrl_M' num2str(motion) '_level1_trial' num2str(trial) '_Single_25-07-04.mat'];
% tviFile = ['./Data/experiment/24-06-21/UUS-iEMG/TVIData_S1_M1_level' num2str(level) '_trial' num2str(trial) '_Dual_24-06-21_' num2str(pp) '.mat'];
load(tviFile);

% 数据预处理
disp('开始数据预处理');
tic;
TVIData = cat(3, zeros(395, 128, 20), TVIData);

% filter the TVI data
TVIDataFilter = TVIData;

% 轴向0.5MHz低通滤波
[Be1, Ae1] = butter(4, 0.5/(7.7*4)*2, 'low');
parfor i = 1:size(TVIDataFilter, 3)
    tmp = TVIDataFilter(:, :, i);
    tmp = filtfilt(Be1, Ae1, tmp);
    TVIDataFilter(:, :, i) = tmp;
end

% 时间5-100Hz带通滤波
% [Be2, Ae2] = butter(4, [5, 100]/fsampu*2);
% for r = 1:size(TVIDataFilter, 1)
%     parfor c = 1:size(TVIDataFilter, 2)
%         tmp = squeeze(TVIDataFilter(r, c, :));
%         tmp = filtfilt(Be2, Ae2, tmp);
%         TVIDataFilter(r, c, :) = tmp;
%     end
% end

% 对每一列降采样
parfor i = 1:size(TVIDataFilter, 3)
    tmp = TVIDataFilter(:, :, i);
    tmp = resample(tmp, 128, 395);
    TVITmp(:, :, i) = tmp;
end
TVIDataFilter = TVITmp;
clear TVITmp;
toc;
disp(['数据预处理用时' num2str(toc)]);

%% USConvBSS算法主体
disp('开始数据迭代');
tic;
[M, N, L] = size(TVIDataFilter);
% 窗口大小
Row = 8; Col = 8;
% 窗口移动距离
dRow = 4; dCol = 4;
% 行数与列数
numRows = (M-Row)/dRow+1;
numCols = (N-Col)/dCol+1;

% 存储空间预分配
tmpB = cell(numRows*numCols, 1);
tmpSources = cell(numRows*numCols, 1);
tmpDecompoPulses = cell(numRows*numCols, 1);
tmpCoV = cell(numRows*numCols, 1);
tmpTwitchesFinal = cell(numRows*numCols, 1);
tmpTwitches = cell(numRows*numCols, 1);

parfor kkk = 1:(numRows*numCols)
    r = ceil(kkk / numCols);
    c = mod(kkk-1, numCols) + 1;

    % 划取处理窗口
    disp(['row=' num2str(r) ',col=' num2str(c)]);
    winRow = (1:Row)+(r-1)*dRow;
    winCol = (1:Col)+(c-1)*dCol;

    winRow(winRow>M) = [];
    winCol(winCol>N) = [];

    TVIDataWin = TVIDataFilter(winRow, winCol, :);
    TVIDataWin = reshape(TVIDataWin, length(winRow)*length(winCol), []);

    % 2.沿着时间进行Z-score
    TVIDataWin = (TVIDataWin - mean(TVIDataWin, 2)) ./ std(TVIDataWin, 0, 2);

    % 3.数据拓展
    eY = extend(TVIDataWin, exFactor);
    eY = eY(:, 1:L);

    % 4.在每个维度上减去均值
    eY = stripmean(eY, 'st');

    % 5.白化
    % 协方差矩阵特征值分解
    [V, D] = eig(cov(eY'));
    [d, idx] = sort(diag(D), 'descend');
    V = V(:, idx);
    D = diag(d);
    % 选取贡献占比70%的特征值
    d = d ./ sum(d);
    cumuSum = cumsum(d);
    ii = find(cumuSum > 0.7, 1);
    % 生成新的特征向量与特征值
    D_new = D(1:ii, 1:ii) - mean(diag(D(ii+1:end, ii+1:end))) * eye(ii);
    V_new = V(:, 1:ii);
    % 白化矩阵WM，采用PCA白化格式
    WM = sqrt(D_new)\V_new';
    % 白化后的数据
    Z = WM * eY;

    % 6.初始化矩阵B
    B = zeros(ii, numCompo);
    twitches = zeros(size(eY, 2), numCompo);
    twitchesFinal = zeros(size(eY, 2), numCompo);
    sources = zeros(size(eY, 2), numCompo);
    decompo_pulses = cell(1, numCompo);
    CoV = zeros(1, numCompo);

    % 7.迭代更新
    for i = 1:numCompo
        iterCount = 0;
        w_new = randn(size(Z, 1), 1);

        while true
            w_old = w_new;
            % 固定点迭代
            w_new = mean(Z.*tanh(w_old'*Z), 2) - mean(sech(w_old'*Z).^2).*w_old;
            % 正交化处理
            w_new = w_new - B*B'*w_new;
            % 归一化处理
            w_new = w_new / norm(w_new);
            % 记录迭代次数
            iterCount = iterCount + 1;
            if abs(w_new'*w_old - 1) < Tolx || iterCount >= 1000
                disp(['r' num2str(r) 'c' num2str(c) ' #' num2str(i) ' 一阶段迭代' num2str(iterCount) '次']);
                break;
            end
        end
        % 一阶段结果存储
        twitches(:, i) = Z' * w_new;

        CoV_new = Inf;
        countcount = 0;
        while true
            CoV_old = CoV_new;
            s = w_new' * Z;
            [source_new, PT, CoV_new, ~] = blindDeconvPeakFinding(s, fsampu, 20, 2, 50, 2);
            w_new = mean(Z(:, PT), 2);
            countcount = countcount + 1;
            if CoV_new > CoV_old
                disp(['r' num2str(r) 'c' num2str(c) ' #' num2str(i) ' 二阶段迭代' num2str(countcount) '次']);
                break;
            end
        end

        % 存储结果
        B(:, i) = w_new;
        sources(:, i) = source_new;
        twitchesFinal(:, i) = s;
        decompo_pulses{i} = PT;
        CoV(i) = CoV_new;
    end

    tmpB{kkk} = B;
    tmpSources{kkk} = sources;
    tmpDecompoPulses{kkk} = decompo_pulses;
    tmpCoV{kkk} = CoV;
    tmpTwitches{kkk} = twitches;
    tmpTwitchesFinal{kkk} = twitchesFinal;
end

toc;
disp(['数据迭代用时' num2str(toc)]);

%% 结果保存

DecompoResults.B = reshape(tmpB, numRows, numCols)';
DecompoResults.sources = reshape(tmpSources, numRows, numCols)';
DecompoResults.twitches = reshape(tmpTwitches, numRows, numCols)';
DecompoResults.twitchesFinal = reshape(tmpTwitchesFinal, numRows, numCols)';
DecompoResults.decompo_pulses = reshape(tmpDecompoPulses, numRows, numCols)';
DecompoResults.CoV = reshape(tmpCoV, numRows, numCols)';

save(['./Data/experiment/25-07-04/M' num2str(motion) 'L1T' num2str(trial) '_USCBSS_compo25.mat'], 'DecompoResults', '-v7.3');

% save(['./Data/experiment/24-06-21/UUS-iEMG/S1M1L' num2str(level) 'T' num2str(trial) '_USCBSS_compo' num2str(numCompo) '_' num2str(pp) '_2s1.mat'], 'DecompoResults', '-v7.3');
    end
end