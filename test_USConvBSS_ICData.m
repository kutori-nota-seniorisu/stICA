% 对超声速度图像进行卷积盲源分离，使用了仿真数据集，采用并行计算
clear; clc; close all;
addpath('./Func');
% 拓展因子
exFactor = 20;
% 迭代容差
Tolx = 1e-4;
% 潜在成分个数
numCompo = 10;
% 超声采样率
fsampu = 1000;

%% 读取数据并进行卷积盲源分离
sub = 16;

%% TVI数据预处理
disp('导入数据');
% 导入TVI数据
disp(['Sub=' num2str(sub)]);
tviFile = ['./Data/experiment/ICdata/R' num2str(sub) '/v_2d_all.mat'];

% disp(['M' num2str(motion) 'L1T' num2str(trial)]);
% tviFile = ['./Data/experiment/25-07-04/TVIData_15000_S_wrl_M' num2str(motion) '_level1_trial' num2str(trial) '_Single_25-07-04.mat'];

% disp(['M1L' num2str(level) 'T' num2str(trial) 'P' num2str(probe)]);
% tviFile = ['./Data/experiment/24-06-21/UUS-iEMG/TVIData_S1_M1_level' num2str(level) '_trial' num2str(trial) '_Dual_24-06-21_' num2str(probe) '.mat'];

load(tviFile);

% 数据预处理
disp('开始数据预处理');
tic;
TVIData = cat(3, zeros(119, 128, 2), v_2d_all);
% TVIData = cat(3, zeros(395, 128, 20), TVIData);

% filter the TVI data
TVIDataFilter = TVIData;
% TVIDataFilter = TVIData(:, :, 2001:end);

% 轴向0.5MHz低通滤波
% [Be1, Ae1] = butter(4, 0.5/(7.7*4)*2, 'low');
% parfor i = 1:size(TVIDataFilter, 3)
%     tmp = TVIDataFilter(:, :, i);
%     tmp = filtfilt(Be1, Ae1, tmp);
%     TVIDataFilter(:, :, i) = tmp;
% end

% 时间5-100Hz带通滤波
[Be2, Ae2] = butter(4, [5, 100]/fsampu*2);
for r = 1:size(TVIDataFilter, 1)
    for c = 1:size(TVIDataFilter, 2)
        tmp = squeeze(TVIDataFilter(r, c, :));
        tmp = filtfilt(Be2, Ae2, tmp);
        TVIDataFilter(r, c, :) = tmp;
    end
end

% 对每一列降采样
% parfor i = 1:size(TVIDataFilter, 3)
%     tmp = TVIDataFilter(:, :, i);
%     tmp = resample(tmp, 128, 395);
%     TVITmp(:, :, i) = tmp;
% end
% TVIDataFilter = TVITmp;
% clear TVITmp;

toc;
disp(['数据预处理用时' num2str(toc)]);

%% USConvBSS算法主体
disp('开始数据迭代');
tic;
[M, N, L] = size(TVIDataFilter);
% 窗口大小
Row = 10; Col = 10;
% 窗口移动距离
dRow = 5; dCol = 5;
% 行数与列数
numRows = ceil((M-Row)/dRow+1);
numCols = ceil((N-Col)/dCol+1);

% 存储空间预分配
tmpB = cell(numRows*numCols, 1);
tmpSources = cell(numRows*numCols, 1);
tmpDecompoPulses = cell(numRows*numCols, 1);
tmpCoV = cell(numRows*numCols, 1);
tmpTwitches = cell(numRows*numCols, 1);

tmpB1 = cell(numRows*numCols, 1);
tmpB2 = cell(numRows*numCols, 1);

for kkk = 1:(numRows*numCols)
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

    figure;
    t=tiledlayout('flow', 'TileSpacing', 'tight', 'Padding', 'compact');
    for iii = 1:size(Z,1)
        nexttile;
        plot(Z(iii,:));
        % xticks(0:2000:20000);
        % xticklabels(0:1:10);
        % kstest(Z(iii,:));
    end
    xlabel(t,'t (s)');
    ylabel(t, 'amplitude')
    sgtitle('白化后的数据，方差占比70%')

    % 6.初始化
    % B对应于numCompo个分离向量
    B = zeros(ii, numCompo);
    % twitches对应于一阶段迭代完成后的收缩曲线
    twitches = zeros(L, numCompo);
    % sources对应于二阶段迭代完成后的IPT曲线
    sources = zeros(L, numCompo);
    % 提取得到的放电脉冲串
    decompo_pulses = cell(1, numCompo);
    % 放电变异率
    CoV = zeros(1, numCompo);
    % B1用来存储一阶段所有的分离向量
    B1 = cell(1, numCompo);
    % B2用来存储二阶段所有的分离向量
    B2 = cell(1, numCompo);

    % 7.迭代更新
    for i = 1:numCompo
        B1{i} = [];
        B2{i} = [];

        iterCount = 0;
        w_new = randn(size(Z, 1), 1);
        B1{i}(:, end+1) = w_new;

        while true
            w_old = w_new;
            % 固定点迭代
            w_new = mean(Z.*tanh(w_old'*Z), 2) - mean(sech(w_old'*Z).^2).*w_old;
            % 正交化处理
            w_new = w_new - B*B'*w_new;
            % 归一化处理
            w_new = w_new / norm(w_new);
            B1{i}(:, end+1) = w_new;
            % 记录迭代次数
            iterCount = iterCount + 1;
            if abs(w_new'*w_old - 1) < Tolx || iterCount >= 1000
                disp(['r' num2str(r) 'c' num2str(c) ' #' num2str(i) ' 一阶段迭代' num2str(iterCount) '次']);
                break;
            end
        end

        figure;
        tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'compact');
        for iii = 1:min(100,size(B1{i}, 2))
            www = B1{i}(:, iii);
            sss = www' * Z;
            nexttile;
            plot(sss);
        end


        % 一阶段结果存储
        twitches(:, i) = Z' * w_new;
        B2{i}(:, end+1) = w_new;

        CoV_new = Inf;
        countcount = 0;
        while true
            CoV_old = CoV_new;
            s = w_new' * Z;
            % MPD=50ms,R=20,nMAD=2
            [source_new, PT, CoV_new, ~] = blindDeconvPeakFinding(s, fsampu, 20, 2, 50, 2);
            w_new = mean(Z(:, PT), 2);
            B2{i}(:, end+1) = w_new;
            countcount = countcount + 1;
            if CoV_new > CoV_old
                disp(['r' num2str(r) 'c' num2str(c) ' #' num2str(i) ' 二阶段迭代' num2str(countcount) '次']);
                break;
            end
        end

        figure;
        tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'compact');
        for iii = 1:min(100,size(B2{i}, 2))
            www = B2{i}(:, iii);
            sss = www' * Z;
            nexttile;
            plot(sss);
        end

        % 存储结果
        B(:, i) = w_new;
        sources(:, i) = source_new;
        decompo_pulses{i} = PT;
        CoV(i) = CoV_new;
    end

    tmpB{kkk} = B;
    tmpSources{kkk} = sources;
    tmpDecompoPulses{kkk} = decompo_pulses;
    tmpCoV{kkk} = CoV;
    tmpTwitches{kkk} = twitches;
    tmpB1{kkk} = B1;
    tmpB2{kkk} = B2;
end

toc;
disp(['数据迭代用时' num2str(toc)]);

%% 结果保存
DecompoResults.B = reshape(tmpB, numRows, numCols)';
DecompoResults.sources = reshape(tmpSources, numRows, numCols)';
DecompoResults.twitches = reshape(tmpTwitches, numRows, numCols)';
DecompoResults.decompo_pulses = reshape(tmpDecompoPulses, numRows, numCols)';
DecompoResults.CoV = reshape(tmpCoV, numRows, numCols)';

% DecompoResults.B1 = reshape(tmpB1, numRows, numCols)';
% DecompoResults.B2 = reshape(tmpB2, numRows, numCols)';

save(['./Data/experiment/ICdata/R' num2str(sub) '/USCBSS_C' num2str(numCompo) 'F.mat'], 'DecompoResults', '-v7.3');
% save(['./Data/experiment/24-06-21/UUS-iEMG/S1M1L' num2str(level) 'T' num2str(trial) '_USCBSS_C' num2str(numCompo) '_' num2str(probe) '.mat'], 'DecompoResults', '-v7.3');
% save(['./Data/experiment/25-07-04/M' num2str(motion) 'L1T' num2str(trial) '_USCBSS_R' num2str(exFactor) 'C' num2str(numCompo) '.mat'], 'DecompoResults', '-v7.3');
disp('数据保存完成！');
