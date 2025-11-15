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
fsampu = 1000;

%% 开启并行池
if isempty(gcp('nocreate'))
    parpool;
end

% 导入TVI数据
for Sub = [16]

    disp(['Sub=' num2str(Sub)]);
    tviFile = ['./Data/experiment/ICdata/R' num2str(Sub) '/v_2d_all.mat'];
    load(tviFile);

    % TVI数据预处理
    disp('开始数据预处理');
    tic;
    TVIData = cat(3, zeros(119, 128, 2), v_2d_all);

    % filter the TVI data
    TVIDataFilter = TVIData;
    % TVIDataFilter = v_2d_all;

    % 时间5-100Hz带通滤波
    [Be2, Ae2] = butter(4, [5, 100]/fsampu*2);
    for r = 1:size(TVIDataFilter, 1)
        parfor c = 1:size(TVIDataFilter, 2)
            tmp = squeeze(TVIDataFilter(r, c, :));
            tmp = filtfilt(Be2, Ae2, tmp);
            TVIDataFilter(r, c, :) = tmp;
        end
    end

    toc;
    disp(['数据预处理用时' num2str(toc)]);

    %%
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
    tmpWFirst = cell(numRows*numCols, 1);
    tmpSourceFirst = cell(numRows*numCols, 1);
    tmpMsg = cell(numRows*numCols, 1);

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
        % figure;
        % plot(TVIData');

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
        % WM = sqrt(inv(D)) * V';
        WM = sqrt(D_new)\V_new';
        % 白化后的数据
        Z = WM * eY;

        % 6.初始化矩阵B
        B = zeros(ii, numCompo);
        source = zeros(size(eY, 2), numCompo);
        decompo_pulses = cell(1, numCompo);
        CoV = zeros(1, numCompo);
        sourcesFirst = zeros(size(eY, 2), numCompo);
        wFirst = zeros(ii, numCompo);

        Msg = {};

        % 7.迭代更新
        for i = 1:numCompo
            iterCount = 0;
            w_new = randn(size(D_new, 1), 1);

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
                    disp(['r=' num2str(r) ',c=' num2str(c) ',i=' num2str(i) '，一阶段迭代完成，本次迭代' num2str(iterCount) '次']);
                    Msg{end+1} = sprintf('r=%d,c=%d,i=%d，一阶段迭代完成，本次迭代%d次', r, c, i, iterCount);
                    break;
                end
            end
            % 一阶段结果存储
            sourcesFirst(:, i) = Z' * w_new;
            wFirst(:, i) = w_new;

            CoV_new = Inf;
            countcount = 0;
            while true
                CoV_old = CoV_new;
                s = w_new' * Z;
                [source_new, PT, CoV_new, ~] = blindDeconvPeakFinding(s, fsampu, 20, 2, 20, 2);
                w_new = mean(Z(:, PT), 2);
                countcount = countcount + 1;
                disp(['r=' num2str(r) ',c=' num2str(c) ',i=' num2str(i) '，二阶段迭代' num2str(countcount) '次']);
                if CoV_new > CoV_old
                    Msg{end+1} = sprintf('r=%d,c=%d,i=%d，二阶段迭代%d次', r, c, i, countcount);
                    break;
                end
            end

            % 存储结果
            B(:, i) = w_new;
            source(:, i) = source_new;
            decompo_pulses{i} = PT;
            CoV(i) = CoV_new;
        end

        tmpB{kkk} = B;
        tmpSources{kkk} = source;
        tmpDecompoPulses{kkk} = decompo_pulses;
        tmpCoV{kkk} = CoV;
        tmpWFirst{kkk} = wFirst;
        tmpSourceFirst{kkk} = sourcesFirst;
        tmpMsg{kkk} = Msg;
    end

    toc;
    disp(['数据迭代用时' num2str(toc)]);

    DecompoResults.B = reshape(tmpB, numCols, numRows)';
    DecompoResults.sources = reshape(tmpSources, numCols, numRows)';
    DecompoResults.decompo_pulses = reshape(tmpDecompoPulses, numCols, numRows)';
    DecompoResults.CoV = reshape(tmpCoV, numCols, numRows)';
    DecompoResults.wFirst = reshape(tmpWFirst, numCols, numRows)';
    DecompoResults.sourceFirst = reshape(tmpSourceFirst, numCols, numRows)';
    DecompoResults.Msg = reshape(tmpMsg, numCols, numRows)';

    savepath = ['./Data/experiment/ICdata/R' num2str(Sub)];
    if ~exist(savepath, 'dir')
        mkdir(savepath);
        disp('创建路径！');
    end
    save([savepath '/USCBSS_compo' num2str(numCompo) '_filter.mat'], 'DecompoResults', '-v7.3');
    disp('数据保存完成！');

end