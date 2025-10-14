% 对超声速度图像进行卷积盲源分离，使用了仿真数据集，采用并行计算
clear; clc; close all;
addpath('./Func');
% 拓展因子
exFactor = 10;
% 迭代容差
Tolx = 1e-4;
% 潜在成分个数
numCompo = 15;

% 开启并行池
if isempty(gcp('nocreate'))
    parpool;
end

% 1.三维数组转二维数组
for setset = 1%:10
    tic;
    datasets_num = num2str(setset);
    % 导入空间源成分
    folderPath = ['./Data/simulation/figure/image_mat' datasets_num '/'];
    fileFormat = '*.mat';
    % 使用dir函数获取文件夹中符合文件格式的文件信息
    fileList = dir(fullfile(folderPath, fileFormat));
    image_all = [];
    % 导入十个空间源成分
    for i = 1:numel(fileList)
        load(fullfile(folderPath, fileList(i).name));
        image_all(:,:,i) = image_data;
    end
    % 高斯噪声图像
    imageGS1 = randn(400,128); image_all(:,:,end+1) = imageGS1*3;
    imageGS2 = randn(400,128); image_all(:,:,end+1) = imageGS2*3;
    Xs = reshape(image_all,size(image_all,1)*size(image_all,2),size(image_all,3));

    % 导入时间源成分
    waveGS1 = randn(4000,1);
    waveGS2 = randn(4000,1);
    for i = 1:10
        load(['./Data/simulation/MU_time_response/TimeCompoDatasets' datasets_num '/Time_component' num2str(i) '.mat']);
        Xt_temp(:,i) = MU_noisy_conv(1:4000);
        clear MU_noisy_conv;
    end
    Xt = [Xt_temp waveGS1 waveGS2];

    % 生成仿真数据
    TVIdata = Xs * Xt';
    TVIdata = reshape(TVIdata, 400, 128, []);

    [M, N, ~] = size(TVIdata);
    % 窗口大小
    Row = 16; Col = 16;
    % 窗口移动距离
    dRow = 8; dCol = 8;
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
            TVIData = TVIdata(winRow, winCol, :);
            TVIData = reshape(TVIData, Row*Col, []);

            % 2.沿着时间进行Z-score
            TVIData = (TVIData - mean(TVIData, 2)) ./ std(TVIData, 0, 2);
            % figure;
            % plot(TVIData');

            % 3.数据拓展
            eY = extend(TVIData, exFactor);
            eY = eY(:, 1:size(TVIData, 2));

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

    save(['./Data/simulation/datasets' datasets_num '/USCBSS_compo' num2str(numCompo) '.mat'], 'DecompoResults');
    % save(['result' datasets_num '.mat'], 'DecompoResults');
    toc;
    disp(['程序用时：' num2str(toc)])
end