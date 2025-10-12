%% 测试多概率密度分布函数的拟合效果
% 导入数据并画出原始图
clear; clc; close all;
addpath('./Func')

%% 矩阵分解结果
flag_decomp = 'NMF';
% 导入数据
filepath = 'F:/EEEMG/stICA/Data/simulation/datasets9/UV_compo12_NMFm.mat';
load(filepath);

%% U与V分布
figure;
for i = 1:size(U_hat,2)
    [f1,x1] = ksdensity(U_hat(:,i));
    % 去除0元素的目的是为了让概率密度分布函数更光滑
    zero_count = find(U_hat(:,i)==0);
    U_temp = U_hat(:,i);
    U_temp(zero_count) = [];
    if ~isempty(U_temp)
        [f,xi] = ksdensity(U_temp);
    else
        f = zeros(21,1);
        xi = -10:10;
    end
    subplot(2,6,i);
    plot(x1,f1);
    hold on;
    plot(xi,f);
    title(['Space #' num2str(i)]);
end
legend('raw','process');
sgtitle(['method=' flag_decomp ',data=U']);
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

% figure;
% for i = 1:size(U_hat,2)
%     subplot(2,6,i);
%     imagesc(reshape(U_hat(:,i),400,128));
% end
% sgtitle(['method=' flag_decomp ',data=U']);
% set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

figure;
for i = 1:size(V_hat,2)
    [f1,x1] = ksdensity(V_hat(:,i));
    % 去除0元素的目的是为了让概率密度分布函数更光滑
    zero_count = find(V_hat(:,i)==0);
    V_temp = V_hat(:,i);
    V_temp(zero_count) = [];
    if ~isempty(V_temp)
        [f,xi] = ksdensity(V_temp);
    else
        f = zeros(21,1);
        xi = -10:10;
    end
    subplot(2,6,i);
    plot(x1,f1);
    hold on;
    plot(xi,f);
    title(['Time #' num2str(i)]);
end
legend('raw','process');
sgtitle(['method=' flag_decomp ',data=V']);
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

%% sech^2函数拟合
figure;
i=1;
zero_count = find(V_hat(:,i)==0);
V_temp = V_hat(:,i);
V_temp(zero_count) = [];
V_temp = V_temp-mean(V_temp);
if ~isempty(V_temp)
    [f,xi] = ksdensity(V_temp);
    % if ~isempty(find(f>1.5))
    %     f = normalize(f,'range');
    % end
else
    f = zeros(21,1);
    xi = -10:10;
end
y_sech = sech(xi).^2;
% subplot(2,6,i);
plot(xi,f);
hold on;
plot(xi,y_sech);
title(['Time #' num2str(i)]);

%% S与T分布
figure;
for i = 1:size(S,2)
    [f,xi] = ksdensity(S(:,i));
    subplot(2,12,i);
    plot(xi,f);
    title(['Space #' num2str(i)]);
    subplot(2,12,i+12);
    imagesc(reshape(S(:,i),400,128));
end
sgtitle(['method=' flag_decomp ',data=S']);
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

figure;
[B,A] = butter(4,[5,25]/2000*2);
T_filter = filtfilt(B,A,T);
T_mean = mean(T_filter);
T_std = std(T_filter);
T_norm = (T_filter-T_mean)./T_std;
for i = 1:size(T,2)
    [f,xi] = ksdensity(T(:,i));
    subplot(2,12,i);
    plot(xi,f);
    title(['Time #' num2str(i)]);
    subplot(2,12,i+12);
    plot(T_norm(:,i));
    % if ismember(i,match_result{8}.Pulses1)
    %     legend('match');
    % end
    %
end
sgtitle(['method=' flag_decomp ',data=T']);
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

% [r1,p1] = corr(S,U_hat,'type','Pearson');
% [r2,p2] = corr(T,V_hat,'type','Pearson');

%% Xs与Xt分布
figure;
for i = 1:10
    [f,xi] = ksdensity(Xs(:,i));
    subplot(2,10,i);
    plot(xi,f);
    title(['Space #' num2str(i)]);
    subplot(2,10,i+10);
    imagesc(reshape(Xs(:,i),400,128));
end
sgtitle(['method=' flag_decomp ',data=Xs']);
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

figure;
for i = 1:10
    [f,xi] = ksdensity(Xt(:,i));
    subplot(2,10,i);
    plot(xi,f);
    title(['Time #' num2str(i)]);
    subplot(2,10,i+10);
    plot(Xt(:,i));
end
sgtitle(['method=' flag_decomp ',data=Xt']);
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

%% 空间pdf拟合
pts = -10:0.05:50;
for i=1:size(U_hat,2)
    % 去除0元素的目的是为了让概率密度分布函数更光滑
    zero_count = find(U_hat(:,i)==0);
    U_temp = U_hat(:,i);
    U_temp(zero_count) = [];
    if ~isempty(U_temp)
        [fs_temp(:,i), xs(:,i)] = ksdensity(U_temp,pts);
        if ~isempty(find(fs_temp(:,i) > 1.5))
            fs(:,i) = normalize(fs_temp(:,i),'range');
        else
            fs(:,i) = fs_temp(:,i);
        end
        [fs_max(i), idx_s(i)] = max(fs(:,i));
    else
        fs(:,i) = 0;
        xs(:,i) = 0;
        fs_max(i) = 0;
        idx_s(i) = 1;
    end
end

% 计算pdf的参数
x = -10:0.05:50;
for i = 1:size(U_hat,2)
    syms b p
    %~ 确定参数a
    a = 2;
    %~ 求解u
    u = floor(xs(idx_s(i),i)/4)*4;
    if u == xs(idx_s(i),i)
        u = u-1;%~ 防止u的值和最大值的横坐标相等导致a=0，即无偏概率分布函数；
    end
    if xs(idx_s(i),i) < params.noise
        x_eqn = a/(sqrt(b^2-a^2)) + u == xs(idx_s(i),i);
        %~ 求解b
        b_temp = double(solve(x_eqn));
        b = abs(b_temp(1)); %~ 取其中一个解的正数即可
        %~ 求解p
        x_eqn_value = a/(sqrt(b^2-a^2)) + u;%~ 概率密度分布函数的峰值对应的横坐标的值
        p_eqn = p*exp(a*(x_eqn_value-u) - b*sqrt((x_eqn_value-u).^2+1)) == fs_max(i);
        p = double(solve(p_eqn));
    else
        a = 1;
        b = 2;
        u = 0;
        p = 1;
    end

    ys(:,i) = p*exp(a*(x-u) - b*sqrt((x-u).^2+1));

    %~ 参数整合
    % params.p(i) = p;
    % params.a(i) = a;
    % params.b(i) = b;
    % params.u(i) = u;
end

%~ 绘制概率密度分布函数图
figure;
for i=1:size(U_hat,2)
    % hold on
    subplot(2,6,i);
    plot(xs(:,i),fs(:,i));
    hold on;
    plot(x,ys(:,i));
    legend('raw data','fit data');
end
sgtitle('Space Fit');
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

%% 时间pdf拟合
pts = -0.5:0.001:1;
for i=1:size(V_hat,2)
    zero_count = find(V_hat(:,i)==0);
    V_temp = V_hat(:,i);
    V_temp(zero_count) = [];
    if ~isempty(V_temp)
        % ksdensity核心平滑密度估计
        [ft_temp(:,i), xt(:,i)] = ksdensity(V_temp,pts);
        % 如果存在大于1.5的密度值，就对估计的概率密度函数进行归一化。
        % if ~isempty(find(ft_temp(:,i) > 1.5))
        %     ft(:,i) = normalize(ft_temp(:,i),'range');
        % else
        ft(:,i) = ft_temp(:,i);
        % end
        % 找到成分i的概率密度函数最大值以及其索引
        [ft_max(i), idx_t(i)] = max(ft(:,i));
    else
        ft(:,i) = 0;
        xt(:,i) = 0;
        ft_max(i) = 0;
        idx_t(i) = 1;
    end
end

x = -0.5:0.001:1;
for i = 1:size(V_hat,2)
    % 拉普拉斯分布
    m = xt(idx_t(i),i);
    c = ft_max(i);

    yt(:,i) = c * exp(-2*c*abs(x-m));

    %~ 参数整合
    % params.m(i) = m;
    % params.c(i) = c;
end

% y = sech(x).^2;

%~ 绘制概率密度分布函数图
figure;
for i=1:size(V_hat,2)
    % hold on
    subplot(2,6,i);
    plot(xt(:,i),ft(:,i));
    hold on;
    plot(x,yt(:,i));
end
legend('raw data','fit data');
sgtitle('Time Fit');
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

%% 量化评估"空间分量直方图分布和拟合函数"之间的相关性
[rho_s,pval_s] = corr(ys,fs,'type','Pearson'); %~ 采用皮尔逊相关系数，rho表示相关系数值，pval<0.05表明拒绝两列之间不存在相关性的假设（即相关）。
[rho_t,pval_t] = corr(yt,ft,'type','Pearson'); %~ 采用皮尔逊相关系数，rho表示相关系数值，pval<0.05表明拒绝两列之间不存在相关性的假设（即相关）。

%%
figure;
for i = 1:10
    zero_count = find(Xt(:,i)==0);
    Xt_temp = Xt(:,i);
    Xt_temp(zero_count) = [];
    [f, xi] = ksdensity(Xt_temp);

    % 峰值的位置
    [~, id] = max(f);
    m = xi(id);
    % 左右比例参数的估计
    bl = sqrt(mean((Xt_temp(Xt_temp<m)).^2));
    br = sqrt(mean((Xt_temp(Xt_temp>=m)).^2));
    % 比值
    gamma_hat = bl/br;
    % 一阶绝对矩
    m1_hat = mean(abs(Xt_temp));
    % 二阶原点矩
    u2_hat = mean(Xt_temp.^2);
    % 比值
    r_hat = (m1_hat-m)^2 / (u2_hat-2*m*m1_hat+m^2);
    % ρ(α)的估计值
    R_hat = r_hat*(gamma_hat^3 +1)*(gamma_hat+1)/((gamma_hat^2+1)^2);
    % 给定候选的形状参数
    alpha = 0:0.001:10;
    alpha(1) = [];
    % 计算所有可能的ρ(α)
    rho_alpha = ((gamma(2./alpha)).^2)./(gamma(1./alpha).*gamma(3./alpha));
    % 寻找距离最近的α
    [~, idx] = min((rho_alpha - R_hat).^2);
    a = alpha(idx);

    x1 = xi(xi<0);
    x2 = xi(xi>=0);
    % [a, b1, b2, m] = AGGD_param_estimate(Xt(:,i));
    f11 = a/((bl+br)*gamma(1/a)) * exp(-(abs(x1-m)/bl).^a);
    f22 = a/((bl+br)*gamma(1/a)) * exp(-(abs(x2-m)/br).^a);
    f2 = [f11,f22];

    subplot(2,5,i);
    plot(xi,f);
    hold on;
    plot(xi,f2);
    title(['Time #' num2str(i)]);
    r(i)=corr(f',f2','Type','Pearson');
    ppp(i,:) = [a,bl,br,m];
end
%%
figure;
pts = -0.05:0.0001:0.05;
for i = 1:12
    zero_count = find(V_hat(:,i)==0);
    V_temp = V_hat(:,i);
    V_temp(zero_count) = [];
    mean_V(i) = mean(V_temp);
    [f_temp,x_temp] = ksdensity(V_temp,pts);
    [~,id] = max(f_temp);
    max_V(i) = x_temp(id);
    rang_V(i) = max(V_temp) - min(V_temp);
    V_temp = V_temp - max_V(i);
    [f,xi] = ksdensity(V_temp,pts);

    x1 = xi(xi<0);
    x2 = xi(xi>=0);
    [a, bl, br] = AGGD_param_estimate(V_temp);
    f11 = a/((bl+br)*gamma(1/a)) * exp(-(abs(x1)/bl).^a);
    f22 = a/((bl+br)*gamma(1/a)) * exp(-(abs(x2)/br).^a);
    f2 = [f11,f22];

    subplot(2,6,i);

    % [counts,edge] = histcounts(V_temp);
    % e = edge(2)-edge(1);
    % new_counts = counts / sum(counts*e);
    % histogram('BinEdges',edge,'BinCounts',new_counts);
    % hold on;

    plot(xi,f);
    hold on;
    plot(xi,f2);

    r(i)=corr(f',f2','Type','Pearson');
    ppp(i,:) = [a,bl,br];

    tt = diff(xi);
    delta = tt(1);
    s(i)=sum(f2*delta);
end
legend('raw data','fit data');