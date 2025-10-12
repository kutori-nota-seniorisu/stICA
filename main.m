% 提取出激活区域的US进行卷积盲源分离
clear; clc; close all;
addpath('./Func');

% 拓展因子
exFactor = 10;
% 迭代容差
Tolx = 1e-4;

datasets_num = '1';
load(['Data/simulation/datasets' datasets_num '/UV_compo12_NC_NMFm.mat']);
%% 1.三维数组转二维数组
% load('Data/simulation/datasets1/TVIdata_compo12_NC.mat');

% 导入空间源成分
folderPath = ['F:/EEEMG/stICA/Data/simulation/figure/image_mat' datasets_num '/'];
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

%% 
for i = 1:size(S, 2)
    % thresh = 0.5 * max(abs(S(:, i)));
    tmpS = abs(S(:, i));
    thresh = 0.5 * (max(tmpS) - min(tmpS)) + min(tmpS);
    tmpInd = find(tmpS > thresh);

    % 没有转换成二维坐标
    actArea{i} = tmpInd;

    % 空间图像二值化
    tmpS(find(tmpS <= thresh)) = 0;
    tmpS(find(tmpS > thresh)) = 1;
    S_bin(:, i) = tmpS;
end

% load('Data/simulation/datasets1/TVIdata_compo12_NC.mat');

for na = 1:length(actArea)
%% 1.提取出激活区域的原始速度曲线
TVIData = TVIdata(actArea{na}, :);

%% 2.沿着时间进行Z-score
TVIData = (TVIData - mean(TVIData, 2)) ./ std(TVIData, 0, 2);

%% 3.数据拓展
eY = extend(TVIData, exFactor);

%% 4.在每个维度上减去均值
eY = stripmean(eY, 'st');
L = size(eY, 2);

%% 5.白化
% 协方差矩阵特征值分解
[V, D] = eig(cov(eY'));
[d, idx] = sort(diag(D), 'ascend');
V = V(:, idx);
D = diag(d);
% 选取贡献占比70%的特征值
d = d ./ sum(d);
d_accu = 0;
for ii = length(d):-1:1
    d_accu = d_accu + d(ii);
    if d_accu > 0.7
        break;
    end
end
D_new = D(ii:end, ii:end) - mean(diag(D(1:ii-1, 1:ii-1)));
V_new = V(:, ii:end);
% 白化矩阵WM，采用PCA白化格式
% WM = sqrt(inv(D)) * V';
WM = sqrt(inv(D_new)) * V_new';
% 白化后的数据
Z = WM * eY;

%% 6.固定点迭代
iterCount = 0;
w_new = randn(size(D_new, 1), 1);

while true
    w_old = w_new;
    % 固定点迭代
    w_new = Z * tanh(w_old' * Z)' / L - mean(sech(w_old' * Z).^2) * w_old;
    % 正交化处理
    % w_new = w_new - B * B' * w_new;
    % 归一化处理
    w_new = w_new / norm(w_new);
    % 记录迭代次数
    iterCount = iterCount + 1;
    if abs(w_new'*w_old - 1) < Tolx
        disp(['第' num2str(na) '个成分一阶段迭代完成，本次迭代' num2str(iterCount) '次']);
        break;
    end
    if iterCount == 1000
        disp(['第' num2str(na) '个成分迭代次数达到上限，迭代终止']);
        break;
    end
end

CoV_new = Inf;
while true
    CoV_old = CoV_new;
    s = w_new' * Z;
    [~, spike_new, PT, CoV_new, ~] = blindDeconvPeakFinding(s, 20, 1.5, 20*2);
    w_new = mean(Z(:, PT), 2);
    if CoV_new > CoV_old
        break;
    end
end

% 存储结果
% B(:, i) = w_new;
% spike{i} = spike_new;
decompo_pulses{na} = PT;
% CoV(i) = CoV_old;

end

