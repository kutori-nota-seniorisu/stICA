%% 对超声速度图像进行卷积盲源分离
clear; clc; close all;
addpath('./Func');
% 拓展因子
exFactor = 10;
% 迭代容差
Tolx = 1e-4;
% 潜在成分个数
numCompo = 15;

%% 1.三维数组转二维数组
% load('Data/simulation/datasets1/TVIdata_compo12_NC.mat');
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

    %%
    for r = 1:400/8-1
        for c = 1:128/8-1
            disp(['row=' num2str(r) ',col=' num2str(c)]);
            winRow = (1:16)+(r-1)*8;
            winCol = (1:16)+(c-1)*8;
            TVIData = TVIdata(winRow, winCol, :);
            TVIData = reshape(TVIData, length(winRow)*length(winCol), []);

            %% 2.沿着时间进行Z-score
            TVIData = (TVIData - mean(TVIData, 2)) ./ std(TVIData, 0, 2);
            % figure;
            % plot(TVIData');

            %% 3.数据拓展
            eY = extend(TVIData, exFactor);
            eY = eY(:, 1:size(TVIData, 2));

            %% 4.在每个维度上减去均值
            eY = stripmean(eY, 'st');
            % L = size(eY, 2);

            %% 5.白化
            % 协方差矩阵特征值分解
            [V, D] = eig(cov(eY'));
            [d, idx] = sort(diag(D), 'descend');
            V = V(:, idx);
            D = diag(d);
            % 选取贡献占比70%的特征值
            d = d ./ sum(d);
            d_accu = 0;
            for ii = 1:length(d)
                d_accu = d_accu + d(ii);
                if d_accu > 0.7
                    break;
                end
            end
            D_new = D(1:ii, 1:ii) - mean(diag(D(ii+1:end, ii+1:end)));
            V_new = V(:, 1:ii);
            % 白化矩阵WM，采用PCA白化格式
            % WM = sqrt(inv(D)) * V';
            WM = sqrt(inv(D_new)) * V_new';
            % 白化后的数据
            Z = WM * eY;

            % figure;
            % plot(Z');

            %% 6.初始化矩阵B
            B = zeros(size(D_new));

            %% 7.迭代更新
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
                sourcesFirst(:, i) = Z' * w_new;
                wFirst(:, i) = w_new;
                % w_new

                CoV_new = Inf;
                while true
                    CoV_old = CoV_new;
                    s = w_new' * Z;
                    [source_new, PT, CoV_new, ~] = blindDeconvPeakFinding(s, 20, 4, 20*2, 2);

                    % figure;
                    % subplot(2,1,1); plot(s); title('s');
                    % subplot(2,1,2); plot(source_new); title('source');
                    % sgtitle(['row=' num2str(r) ' col=' num2str(c) ' cp=' num2str(i)]);
                    % drawnow;

                    w_new = mean(Z(:, PT), 2);
                    if CoV_new >= CoV_old
                        break;
                    end
                end

                % 存储结果
                B(:, i) = w_new;
                source(:, i) = source_new;
                decompo_pulses{i} = PT;
                CoV(i) = CoV_new;
            end
            DecompoResults.B{r, c} = B;
            DecompoResults.source{r, c} = source;
            DecompoResults.decompo_pulses{r, c} = decompo_pulses;
            DecompoResults.CoV{r, c} = CoV;
            DecompoResults.wFirst{r, c} = wFirst;
            DecompoResults.sourceFirst{r, c} = sourcesFirst;

            % close all;
        end
    end
    save(['./Data/simulation/datasets' datasets_num '/USCBSS_compo' num2str(numCompo) '.mat'], 'DecompoResults');
    % save(['result' datasets_num '.mat'], 'DecompoResults');
        toc;
    disp(['程序用时：' num2str(toc)])
end