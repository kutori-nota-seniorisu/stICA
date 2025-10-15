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

% 开启并行池
if isempty(gcp('nocreate'))
    parpool;
end

%% TVI数据预处理
for level = 1:2
    for trial = 1:2
        for pp = 1:2
            % 导入TVI数据
            tviFile = ['./Data/experiment/24-06-21/UUS-iEMG/TVIData_S1_M1_level' num2str(level) '_trial' num2str(trial) '_Dual_24-06-21_' num2str(pp) '.mat'];
            load(tviFile);
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
            [Be2, Ae2] = butter(4, [5, 100]/fsampu*2);
            for r = 1:size(TVIDataFilter, 1)
                parfor c = 1:size(TVIDataFilter, 2)
                    tmp = squeeze(TVIDataFilter(r, c, :));
                    tmp = filtfilt(Be2, Ae2, tmp);
                    TVIDataFilter(r, c, :) = tmp;
                end
            end
            % 对每一列降采样
            TVITmp = zeros(128, 128, 15000);
            parfor i = 1:size(TVIDataFilter, 3)
                tmp = TVIDataFilter(:, :, i);
                tmp = resample(tmp, 128, 395);
                TVITmp(:, :, i) = tmp;
            end
            TVIDataFilter = TVITmp;
            clear TVITmp;

            %%
            [M, N, ~] = size(TVIDataFilter);
            % 窗口大小
            Row = 8; Col = 8;
            % 窗口移动距离
            dRow = 4; dCol = 4;
            % 行数与列数
            numRows = (M-Row)/dRow+1;
            numCols = (N-Col)/dCol+1;

            % 存储空间预分配
            % DecompoResults = struct();
            % DecompoResults.B = cell(numRows, numCols);
            % DecompoResults.sources = cell(numRows, numCols);
            % DecompoResults.decompo_pulses = cell(numRows, numCols);
            % DecompoResults.CoV = cell(numRows, numCols);
            % DecompoResults.wFirst = cell(numRows, numCols);
            % DecompoResults.sourceFirst = cell(numRows, numCols);
            tmpB = cell(numRows, numCols);
            tmpSources = cell(numRows, numCols);
            tmpDecompoPulses = cell(numRows, numCols);
            tmpCoV = cell(numRows, numCols);
            tmpWFirst = cell(numRows, numCols);
            tmpSourceFirst = cell(numRows, numCols);

            for r = 1:numRows
                parfor c = 1:numCols
                    disp(['row=' num2str(r) ',col=' num2str(c)]);
                    winRow = (1:Row)+(r-1)*dRow;
                    winCol = (1:Col)+(c-1)*dCol;

                    TVIDataWin = TVIDataFilter(winRow, winCol, :);
                    TVIDataWin = reshape(TVIDataWin, Row*Col, []);

                    % 2.沿着时间进行Z-score
                    TVIDataWin = (TVIDataWin - mean(TVIDataWin, 2)) ./ std(TVIDataWin, 0, 2);
                    % figure;
                    % plot(TVIData');

                    % 3.数据拓展
                    eY = extend(TVIDataWin, exFactor);
                    eY = eY(:, 1:size(TVIDataWin, 2));

                    % 4.在每个维度上减去均值
                    eY = stripmean(eY, 'st');
                    % L = size(eY, 2);

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
                    D_new = D(1:ii, 1:ii) - mean(diag(D(ii+1:end, ii+1:end)));
                    V_new = V(:, 1:ii);
                    % 白化矩阵WM，采用PCA白化格式
                    % WM = sqrt(inv(D)) * V';
                    WM = sqrt(inv(D_new)) * V_new';
                    % 白化后的数据
                    Z = WM * eY;

                    % 6.初始化矩阵B
                    B = zeros(ii, numCompo);
                    source = zeros(size(eY, 2), numCompo);
                    decompo_pulses = cell(1, numCompo);
                    CoV = zeros(1, numCompo);
                    sourcesFirst = zeros(size(eY, 2), numCompo);
                    wFirst = zeros(ii, numCompo);

                    % 7.迭代更新
                    for i = 1:numCompo
                        iterCount = 0;
                        w_new = randn(size(D_new, 1), 1);

                        while true
                            w_old = w_new;
                            % 固定点迭代
                            w_new = Z * tanh(w_old' * Z)' / size(eY, 2) - mean(sech(w_old' * Z).^2) * w_old;
                            % 正交化处理
                            w_new = w_new - B * B' * w_new;
                            % 归一化处理
                            w_new = w_new / norm(w_new);
                            % 记录迭代次数
                            iterCount = iterCount + 1;
                            if abs(w_new'*w_old - 1) < Tolx
                                disp(['compo=' num2str(i) '，一阶段迭代完成，本次迭代' num2str(iterCount) '次']);
                                break;
                            end
                            if iterCount == 10000
                                disp(['compo=' num2str(i) '，一阶段迭代达到上限，迭代终止']);
                                break;
                            end
                        end
                        % 一阶段结果存储
                        sourcesFirst(:, i) = Z' * w_new;
                        wFirst(:, i) = w_new;
                        % tmpSourceFirst{r, c}(:, i) = Z' * w_new;
                        % tmpWFirst{r, c}(:, i) = w_new;

                        CoV_new = Inf;
                        while true
                            CoV_old = CoV_new;
                            s = w_new' * Z;
                            [source_new, PT, CoV_new, ~] = blindDeconvPeakFinding(s, 20, 4, 20*2, 2);
                            w_new = mean(Z(:, PT), 2);
                            if CoV_new > CoV_old
                                break;
                            end
                        end

                        % 存储结果
                        B(:, i) = w_new;
                        source(:, i) = source_new;
                        decompo_pulses{i} = PT;
                        CoV(i) = CoV_new;
                        % tmpB{r, c}(:, i) = w_new;
                        % tmpSources{r, c}(:, i) = source_new;
                        % tmpDecompoPulses{r, c}{i} = PT;
                        % tmpCoV{r, c}(i) = CoV_new;
                    end

                    tmpB{r, c} = B;
                    tmpSources{r, c} = source;
                    tmpDecompoPulses{r, c} = decompo_pulses;
                    tmpCoV{r, c} = CoV;
                    tmpWFirst{r, c} = wFirst;
                    tmpSourceFirst{r, c} = sourcesFirst;

                    % close all;
                end
            end

            DecompoResults.B = tmpB;
            DecompoResults.sources = tmpSources;
            DecompoResults.decompo_pulses = tmpDecompoPulses;
            DecompoResults.CoV = tmpCoV;
            DecompoResults.wFirst = tmpWFirst;
            DecompoResults.sourceFirst = tmpSourceFirst;

            save(['./Data/experiment/24-06-21/UUS_iEMG/S1M1L' num2str(level) 'T' num2str(trial) 'P' num2str(pp) '_USCBSS_compo' num2str(numCompo) '.mat'], 'DecompoResults')
        end
    end
end