%% 对stICA后的单个成分进行卷积盲源分离
% clear; clc; close all;
% addpath('./Func');
% datasets_num = '1';
% load(['Data/simulation/datasets' datasets_num '/UV_compo12_NC_NMFm0.9.mat']);
% 拓展因子
R = 20;
% nMAD
nMAD = 2;
% MPD(ms)
MPD = 50;
% 超声采样率
fsampu = 2000;
% 代价函数
mode = 2;
%%
www = B2(:, 1);
S = www'*Z;
figure;
plot(S);
% 迭代容差
Tolx = 1e-4;
% 拓展信号
eS = extend(S, R);
eS = eS(:, 1:size(S, 2));
% 减去均值
eS = eS - mean(eS, 2);
% figure;
% tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'compact');
% for iii=1:size(eS,1)
%     nexttile;
%     plot(eS(iii,:));
% end
% 协方差矩阵特征值分解
[V, D] = eig(cov(eS'));
[d,idx]=sort(diag(D),'descend');
V=V(:,idx);
D=diag(d);%+1e-6*max(d)

% d = d ./ sum(d);
% cumuSum = cumsum(d);
% ii = find(cumuSum > 0.7, 1);
% % 生成新的特征向量与特征值
% D_new = D(1:ii, 1:ii) - mean(diag(D(ii+1:end, ii+1:end))) * eye(ii);
% V_new = V(:, 1:ii);
% WM = sqrt(D_new)\V_new';

% 白化矩阵WM，采用PCA白化格式
WM = sqrt(D)\V';
% 白化后的数据
Z_CBSS = WM*eS;
figure;
t=tiledlayout('flow', 'TileSpacing', 'tight', 'Padding', 'compact');
for iii=1:size(Z_CBSS,1)
    nexttile;
    plot(Z_CBSS(iii,:));
end
% xlabel(t,'t (s)');
% ylabel(t, 'amplitude')
sgtitle('白化后数据')
%%
% 迭代更新
w_new = randn(R, 1);
% 迭代计数器
iterCount = 0;
figure;
tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'compact');
while true
    w_old = w_new;
    % 固定点迭代
    if mode == 1
        w_new = mean(Z_CBSS.*log(cosh(w_old'*Z_CBSS)), 2) - mean(tanh(w_old'*Z_CBSS)).*w_old;
    elseif mode == 2
        % w_new = Z * 1/6 * (w_old' * Z).^3' / size(eS, 2) - mean(1/2*(w_old' * Z).^2) * w_old;
        w_new = mean(Z_CBSS.*1/2.*(w_old'*Z_CBSS).^2, 2) - mean(w_old'*Z_CBSS).*w_old;
    end
    % 归一化处理
    w_new = w_new / norm(w_new);
    sss=w_new'*Z_CBSS;
    nexttile;
    plot(sss);
    % 记录迭代次数
    iterCount = iterCount + 1;
    if abs(w_new'*w_old - 1) < Tolx || iterCount >= 100
        disp(['迭代完成，本次迭代' num2str(iterCount) '次']);
        break;
    end
end

% 提取放电时刻
source = w_new' * Z_CBSS;
[~, PT] = findpeaks(source, 'MinPeakHeight', nMAD*mad(source), 'MinPeakDistance', MPD*(fsampu/1000));
figure;
plot(source);
hold on;
plot(PT, source(PT), 'ro');

% 计算变异系数
CoV = std(diff(PT)) / mean(diff(PT));