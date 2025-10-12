% 对CBSS的参数进行优化
clear; clc; close all;
addpath('./Func');
% 选出最优 R - 40 nMAD - 1.5 MPD - 10
% 使用下来这个组合的效果并不好，牺牲了匹配个数达到的高RoA不具有说服力。
% R - 20 nMAD - 1.5 MPD - 20的组合还可以
% 拓展因子
exFactor = 10:5:60;%10:5:40;
% nMAD
nMAD = 1:0.5:4;
% MPD(ms)
MPD = 10:5:40;

for ii = 1:length(exFactor)
    for jj = 1:length(nMAD)
        for kk = 1:length(MPD)
            for setset = 1:10

                for i = 1:10
                    load(['./Data/simulation/MU_time_response/TimeCompoDatasets' num2str(setset) '/Time_component' num2str(i) '.mat']);
                    [source, PT, ~, w_new] = blindDeconvPeakFinding(MU_noisy_conv, exFactor(ii), nMAD(jj), MPD(kk)*2);
                    decompo_pulses{i} = PT;
                end
                %% 导入数据
                % load(['Data/simulation/datasets' num2str(setset) '/UV_compo12_NC_NMFm.mat']);
                
                %% spike train extract
                % decompo_pulses = {};
                % for i = 1:12
                %     [source, ~, PT, ~, w_new] = blindDeconvPeakFinding(T(:, i)', exFactor(ii), nMAD(jj), MPD(kk)*2);
                %     decompo_pulses{i} = PT;
                %     % T_new(:, i) = source';
                %     % W(:, i) = w_new;
                % end

                %% 时间成分匹配
                load(['Data/simulation/MU_time_response/TimeCompoDatasets' num2str(setset) '/ipulses.mat']);

                matchresult_time_raw = [];
                % 计算RoA
                for i = 1:length(decompo_pulses)
                    for j = 1:length(ipulses)
                        [Array1, Array2] = meshgrid(decompo_pulses{i}, ipulses{j});
                        diff_values = Array1 - Array2;
                        valid_elements = diff_values <= (exFactor(ii)+10)*2 & diff_values >= 0;
                        count = sum(valid_elements(:));
                        r = count/(length(decompo_pulses{i})+length(ipulses{j})-count);
                        if r > 1
                            r = 1;
                        end
                        spike_ROA_matrix(i, j) = r;
                        matchresult_time_raw(end+1,:) = [i, j, r];
                    end
                end

                matchresult_time = matchresult_time_raw;
                for mu = 1:length(decompo_pulses)
                    tmpInd = find(matchresult_time(:,1) == mu);
                    if length(tmpInd) > 1
                        [~, tmpInd2] = max(matchresult_time(tmpInd,3));
                        deleteInd = setdiff(1:length(tmpInd), tmpInd2);
                        matchresult_time(tmpInd(deleteInd), :) = [];
                    end
                end
                for mu = 1:length(ipulses)
                    tmpInd = find(matchresult_time(:,2) == mu);
                    if length(tmpInd) > 1
                        [~, tmpInd2] = max(matchresult_time(tmpInd,3));
                        deleteInd = setdiff(1:length(tmpInd), tmpInd2);
                        matchresult_time(tmpInd(deleteInd), :) = [];
                    end
                end
                % matchresult_time_raw = array2table(matchresult_time_raw, 'VariableNames', {'decomp', 'ref', 'time'});
                matchresult_time = array2table(matchresult_time, 'VariableNames', {'decomp', 'ref', 'time'});

                %% 空间成分匹配
                % % 计算相关系数并选出最大匹配
                % corr_matrix_space = corr(S,Xs);
                % corr_matrix_space(:,end-1:end) = [];
                % [max_values_space, indices_space] = max(abs(corr_matrix_space),[],2);
                % 
                % % 在重复匹配中找最大匹配
                % matchresult_space_raw = [(1:12)',indices_space, max_values_space];
                % matchresult_space = matchresult_space_raw;
                % for mu = 1:size(S,2)
                %     tmpInd = find(matchresult_space(:,2) == mu);
                %     if length(tmpInd) > 1
                %         [~,tmpInd2] = max(matchresult_space(tmpInd,3));
                %         deleteInd = setdiff(1:length(tmpInd),tmpInd2);
                %         matchresult_space(tmpInd(deleteInd),:) = [];
                %     end
                % end
                % matchresult_space = array2table(matchresult_space, 'VariableNames', {'decomp', 'ref', 'space'});

                %% 时空结果匹配
                % matchresult_final = innerjoin(matchresult_space, matchresult_time, 'Keys', {'decomp', 'ref'});
                % matchresult_final_all{ii, jj, kk}{setset} = matchresult_final;
                matchresult_final_all{ii, jj, kk}{setset} = matchresult_time;
            end
        end
    end
end

%% 随着R的增大，RoA mean增大，但是这种增大是由可以匹配上的MU变少，甚至无法匹配上导致的。
for ii = 1:size(matchresult_final_all, 1)
    exFactorMatch = [];
    for jj = 1:size(matchresult_final_all, 2)
        for kk = 1:size(matchresult_final_all, 3)
            for setset = 1:10
                % matchresult_final_all{ii, jj, kk}{setset}
                exFactorMatch = vertcat(exFactorMatch, matchresult_final_all{ii, jj, kk}{setset});
            end
        end
    end
    exFactorMatchNum(ii) = size(exFactorMatch, 1);
    exFactorMatchMean(ii) = mean(exFactorMatch.time*100);
    exFactorMatchSE(ii) = std(exFactorMatch.time*100)/sqrt(size(exFactorMatch, 1));
end
figure;
yyaxis left;
errorbar(exFactorMatchMean, exFactorMatchSE);
xticklabels(10:5:60);
hold on;
yyaxis right;
plot(exFactorMatchNum);

%% 在某个确定的R情况下，绘制nMAD与MPD对识别准确率的影响热力图
ii = size(matchresult_final_all, 1);
for jj = 1:size(matchresult_final_all, 2)
    for kk = 1:size(matchresult_final_all, 3)
        mmMatch = [];
        for setset = 1:10
            mmMatch = vertcat(mmMatch, matchresult_final_all{ii, jj, kk}{setset});
        end
        mmMatchNum(jj, kk) = size(mmMatch, 1);
        mmMatchMean(jj, kk) = mean(mmMatch.time*100);
        mmMatchSE(jj, kk) = std(mmMatch.time*100)/sqrt(size(mmMatch, 1));
    end
end
figure('Position', [100, 100, 1200, 500]);
subplot(1, 3, 1);
imagesc(mmMatchMean);
set(gca, 'YDir', 'normal')
colormap hot; colorbar;
yticklabels(nMAD);
xticklabels(MPD);
title('Mean Rate of Agreement (%)');
xlabel('Minimal peak distance (MPD)');
ylabel('Number of MAD (nMAD)');

subplot(1, 3, 2);
imagesc(mmMatchSE);
set(gca, 'YDir', 'normal')
colormap hot; colorbar;
yticklabels(nMAD);
xticklabels(MPD);
title('Standard Error of RoA (%)');
xlabel('Minimal peak distance (MPD)');
ylabel('Number of MAD (nMAD)');

subplot(1, 3, 3);
imagesc(mmMatchNum);
set(gca, 'YDir', 'normal')
colormap hot; colorbar;
yticklabels(nMAD);
xticklabels(MPD);
title('Num of MUs');
xlabel('Minimal peak distance (MPD)');
ylabel('Number of MAD (nMAD)');

%% 精细化热力图，类似于等高线
% 2. 创建精细的网格用于插值（关键步骤！）
[MPD_coarse, nMAD_coarse] = meshgrid(MPD, nMAD); % 原始粗糙网格
[MPD_fine, nMAD_fine] = meshgrid(linspace(min(MPD), max(MPD), 50), ...
                                 linspace(min(nMAD), max(nMAD), 50)); % 精细网格

% 3. 插值：将粗糙网格的数据插值到精细网格上
mean_roa_fine = interp2(MPD_coarse, nMAD_coarse, mmMatchMean, MPD_fine, nMAD_fine, 'spline');
std_error_fine = interp2(MPD_coarse, nMAD_coarse, mmMatchSE, MPD_fine, nMAD_fine, 'spline');
num_fine = interp2(MPD_coarse, nMAD_coarse, mmMatchNum, MPD_fine, nMAD_fine, 'spline');

% 4. 绘制 Fig. 4B: Mean RoA 地形图
figure('Position', [100, 100, 1200, 500]);

% 绘制填充等高线图
subplot(1, 3, 1);
contourf(MPD_fine, nMAD_fine, mean_roa_fine, 20, 'LineColor', 'none'); % 20条等高线，隐藏线
colormap('parula');
colorbar;
title('Mean Rate of Agreement (%)');
xlabel('Minimal peak distance (MPD)');
ylabel('Number of MAD (nMAD)');

subplot(1, 3, 2);
contourf(MPD_fine, nMAD_fine, std_error_fine, 20, 'LineColor', 'none');
colormap('hot');
colorbar;
title('Standard Error of RoA (%)');
xlabel('Minimal peak distance (MPD)');
ylabel('Number of MAD (nMAD)');

subplot(1, 3, 3);
contourf(MPD_fine, nMAD_fine, num_fine, 20, 'LineColor', 'none');
colormap('hot');
colorbar;
title('Num of MUs');
xlabel('Minimal peak distance (MPD)');
ylabel('Number of MAD (nMAD)');

sgtitle('Parameter Evaluation (Contourf Method)');

%% 非常哇塞的图。创建三维切片图和等值面图。阅读起来比较麻烦
% 对每一种参数组合进行统计计算
for ii = 1:size(matchresult_final_all, 1)
    for jj = 1:size(matchresult_final_all, 2)
        for kk = 1:size(matchresult_final_all, 3)
            match = [];
            for setset = 1:10
                % matchresult_final_all{ii, jj, kk}{setset}
                match = vertcat(match, matchresult_final_all{ii, jj, kk}{setset});
            end
            matchNum(ii,jj,kk) = size(match, 1);
            matchMean(ii,jj,kk) = mean(match.time*100);
            matchSE(ii,jj,kk) = std(match.time*100)/sqrt(size(match, 1));
        end
    end
end
matchMean(isnan(matchMean)) = 100;

param1_values = exFactor;
param2_values = nMAD;
param3_values = MPD;
% 1. 创建网格
[P1, P2, P3] = meshgrid(param2_values, param1_values, param3_values);

% 2. 绘制等值面图 (Isosurface) - 显示特定性能阈值的“外壳”
figure;
hold on;
patch(isosurface(P1, P2, P3, matchMean, 95), 'FaceColor', [0.8,0.8,1.0], 'FaceAlpha', 0.4); % 绘制性能=85的等值面
% ylabel('Param1 (R)');
% xlabel('Param2 (nMAD)');
% zlabel('Param3 (MPD)');
% title('Isosurface of Performance = 85%');
% axis equal;
% grid on;
% colorbar;

% 3. 绘制切片图 (Slice Planes) - 查看内部结构
% figure;
slice(P1, P2, P3, matchMean, [], [], param3_values);
patch(isosurface(P1, P2, P3, matchNum, 90), 'FaceColor', [0.4,0.4,0.9], 'FaceAlpha', 0.9); % 绘制性能=85的等值面

% 在 param3 的每个值处创建一个切片
ylabel('Param1 (R)');
xlabel('Param2 (nMAD)');
zlabel('Param3 (MPD)');
title('Slice Planes through the Parameter Space');
% shading interp; % 平滑着色
colormap('jet');
colorbar;
% 3D视角
view(3)

% 3. 绘制切片图 (Slice Planes) - 查看内部结构
figure;
slice(P1, P2, P3, matchNum, [], [], param3_values);
% 在 param3 的每个值处创建一个切片
ylabel('Param1 (R)');
xlabel('Param2 (nMAD)');
zlabel('Param3 (MPD)');
title('Slice Planes through the Parameter Space');
shading interp; % 平滑着色
colormap('jet');
colorbar;

figure;
plot(matchMean(:));
hold on;
plot(matchNum(:));
xline(1:77:539);
legend('Mean', 'Num');