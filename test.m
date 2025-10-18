clear; clc; close all;
%%
datasets_num = '2';
% 拓展因子
exFactor = 20;
% nMAD
nMAD = 1.5;
% MPD(ms)
MPD = 20;
for i = 1:10
    load(['./Data/simulation/MU_time_response/TimeCompoDatasets' datasets_num '/Time_component' num2str(i) '.mat']);
    [source, ~, PT, ~, w_new] = blindDeconvPeakFinding(MU_noisy_conv, exFactor, nMAD, MPD*2);
    decompo_pulses{i} = PT;
end

%% 绘制单次实验的一个拟合曲线
compo = 1;
figure;
subplot(2,1,1);
plot(params.x{compo},params.f{compo},'Color','b');
dif = diff(params.x{compo});
delta = dif(1);
s = sum(params.f{compo}*delta);
xx = -5:0.001:5;
yy = sech(xx).^2;
hold on;
plot(xx,yy,'Color','r');
ss = sum(yy*0.001);
% figure;
hold on;
[counts,edges] = histcounts(Xt(:,compo),'BinWidth',delta);
histogram('BinCounts',counts/(4000*delta),'BinEdges',edges);
ylabel('概率密度');
legend('ksdensity','sech^2','hist');
title('Xt vs sech^2');

% figure;
subplot(2,1,2);
tt = 0:0.0005:2;
tt(1) = [];
plot(tt,Xt(:,compo)');
xlabel('time');
ylabel('velocity');
title('twitch curve');

%% 初步筛选
% 保留的脉冲串
decompoPulseAll = {};
% 保留的估计源
decompoSourceAll = [];
% 保留的一阶段估计源
decompoSourceFirstAll = [];
% 保留的CoV
decompoCoVAll = [];
for r = 1:400/8-1
    for c = 1:128/8-1

        % disp(['row=' num2str(r) ',col=' num2str(c)]);
        tmpPulses = DecompoResults.decompo_pulses{r, c};
        tmpSourcesFirst = DecompoResults.sourceFirst{r, c};
        tmpSources = DecompoResults.sources{r, c};
        % tmpCoV = DecompoResults.CoV{r, c};
        % rho = corr(sources);
        for mu = 1:10
            if isempty(tmpPulses{mu})
                continue;
            end
            decompoPulseAll{end+1} = tmpPulses{mu};
            decompoSourceAll(:, end+1) = tmpSources(:, mu);
            decompoSourceFirstAll(:, end+1) = tmpSourcesFirst(:, mu);
            % decompoCoVAll(end+1) = tmpCoV(mu);
        end
    end
end

for i = 1%:size(decompoSourceAll, 2)
    for j = i:size(decompoSourceAll, 2)
        [corrVals, ~] = xcorr(decompoSourceAll(:, i), decompoSourceAll(:, j), 'coeff');
        [~, idx] = max(abs(corrVals));
        ccc(i, j) = corrVals(idx);
    end
end
%%
d1 = 8;
d2 = 4;
ref = 9;
figure;
subplot(5,1,1);
plot(-T(:, d1));
title('BPM')
subplot(5,1,2);
plot(-T_norm(:, d1));
hold on;
scatter(decompo_pulses{d1},zeros(length(decompo_pulses{d1})), 1000, 'black', '|');
title('BPM normalization')
subplot(5,1,3);
plot(-decompoSourceFirstAll(:, d2));
% hold on;
% scatter(decompoPulseAll{d2}, zeros(length(decompoPulseAll{d2})), 1000, "black", '|');
title('USCBSS First Step')
subplot(5,1,4);
plot(decompoSourceAll(:, d2));
hold on;
scatter(decompoPulseAll{d2}, zeros(length(decompoPulseAll{d2})), 1000, "black", '|');
title('USCBSS Second Step')
subplot(5,1,5);
plot(Xt(:, ref));
hold on;
scatter(ipulses{ref}, zeros(length(ipulses{ref})), 1000, "black", '|');
title('Ref Source')
set(gcf,'unit','normalized','position',[0.05,0.1,0.9,0.6]);

%%
plotDecomps({ipulses{6}, decompoPulseAll{12}}, [], 2000, 0, 0, []);
[c, l] = xcorr(decompoSourceAll(:, 7), decompoSourceAll(:, 10), 'coeff');
figure;
plot(c);
%% 将USConvBSS的结果与BPM进行对比
figure;
subplot(2,1,1);
hold on;
for i=1:10
    % load(['Results/result' num2str(i) '_Update.mat']);
    load(['Results/datasets' num2str(i) '_result.mat']);
    boxplot(matchresult_time.time, 'Positions', i-0.2, 'Colors', 'r');
    hold on;
    timeMean(i) = mean(matchresult_time.time);
    numMU(i) = size(matchresult_time, 1);
end
plot((1:10)-0.2, timeMean, 'r-*');
hold on;

for i=1:10
    % load(['Results/result' num2str(i) '_Update.mat']);
    load(['Results/datasets' num2str(i) '_resultnew.mat']);
    boxplot(matchresult_time.time, 'Positions', i, 'Colors', 'b');
    hold on;
    timeMean2(i) = mean(matchresult_time.time);
    numMU2(i) = size(matchresult_time, 1);
end
plot((1:10), timeMean2, 'b-*');
hold on;

load('Results/compo12_NC_NMFm0.9_BPM_10sets.mat');
for i = 1:10
    % space_mean_BPM(i) = mean(matchresult_final_all{i}.space);
    time_mean_BPM(i) = mean(matchresult_final_all{i}.time);
    No_BPM(i) = size(matchresult_final_all{i}, 1);
    % time_No_NMFm(i) = size(matchresult_time_all{i}, 1);
    % space_No_NMFm(i) = size(matchresult_space_all{i}, 1);
end
match_BPM = matchresult_final_all;
for i = 1:10
    boxplot(match_BPM{i}.time, 'Positions', i+0.2, 'Colors', 'g');
    hold on;
end
plot((1:10)+0.2,time_mean_BPM,'g-*');
ylabel('RoA (Time)');
xticks(1:10)
xticklabels(1:10);
legend('USCBSS', 'USCBSS2', 'BPM');

subplot(2,1,2);
plot(1:10, numMU, 'r-*');
hold on;
plot(1:10, numMU2, 'b-*');
hold on;
plot(1:10, No_BPM, 'g-*');
ylabel('No. MU');
legend('USCBSS', 'USCBSS2', 'BPM');

sgtitle('USCBSS vs BPM')
set(gcf,'unit','normalized','position',[0.05,0.1,0.9,0.6]);
%% 梯度1：双曲正割
y_t = V_hat * Wt;
% sigma_dd = -2 * sech(y_t).^3 .* sinh(y_t);
% sigma_d = sech(y_t).^2;
tmp = -2*sech(y_t)'.*sinh(y_t)';
sigma = reshape(tmp, params.k, 1, []);
x = reshape(V_hat', 1, params.k, []);
% tmp3 = sigma .* x;
phi_yx = mean(sigma .* x, 3);
% phi_y = mean(sigma_dd ./ sigma_d);
% x = mean(V_hat);
% nabla_h = inv(Wt') + phi_y' * x;
%% 梯度2：带有偏度
p = params.p;
a = params.a;
b = params.b;
u = params.u;
y_s = U_hat * Ws - u;
% sigma_d = p .* exp(a .* y_s - b .* sqrt(y_s.^2 + 1));
% sigma_dd = sigma_d .* (a - b .* y_s ./ sqrt(y_s.^2 + 1));
tmp = (a - b .* y_s ./ sqrt(y_s.^2 + 1))';
sigma = reshape(tmp, params.k, 1, []);
x = reshape(U_hat', 1, params.k, []);
phi_yx = mean(sigma .* x, 3);
% phi_y = mean(sigma_dd ./ sigma_d);
% x = mean(U_hat);
% nabla_h = inv(Ws') + phi_y' * x;
%% 梯度3：插值函数
pts = -3 : 0.005 : 3;
for i = 1:size(Xt,2)
    [f, x] = ksdensity(Xt(:,i), pts);
    pp = spline(x, f);
    pp_d = fnder(pp, 1);
    params.pp(i) = pp;
    params.pp_d(i) = pp_d;
end

y_t = V_hat * Wt;
idx1 = find(y_t(:,1)<params.pp(1).breaks(1));
p_t(:,1) = ppval(params.pp(1), y_t(:,1));
p_t(idx1,1) = 0;
idx2 = find(y_t(:,1)>params.pp(1).breaks(end))

% pts = -3 : 0.001 : 3;
% [y,x] = ksdensity(Xt(:,1),pts);
figure;
scatter(x,y(:,1),'*');
hold on;
scatter(y_t(:,1),squeeze(s(1,1,:)));
legend('Sample Points','spline');

figure;
scatter(x,y(:,2),'*');
hold on;
scatter(y_t(:,2),squeeze(s(2,2,:)));
legend('Sample Points','spline');

for i = 1:12
    pp(i) = spline(x,y(:,i)');
    s2(i,:) = ppval(pp(i),y_t(:,i));
    pp_deriv1(i) = fnder(pp(i), 1);
    s2_deriv1(i,:) = ppval(pp_deriv1(i), y_t(:,i));
    % s2(i,:) = spline(x,y(:,i)',y_t(:,i)');
end
figure;
scatter(x,y(:,1),'*');
hold on;
scatter(y_t(:,1),s2(1,:));
legend('Sample Points','spline');

figure;
scatter(x,y(:,2),'*');
hold on;
scatter(y_t(:,2),s2(2,:));
legend('Sample Points','spline');

phi_y = mean(s2_deriv1 ./ s2);
x = mean(V_hat);
nabla_h = int(Wt') + phi_y' * x;
%% 绘制匹配成分图
decomp = 3;
ref = 5;

figure;
subplot(1,2,1);
imagesc(reshape(S(:,decomp),400,128));
title(['S' num2str(decomp)]);
subplot(1,2,2);
imagesc(reshape(Xs(:,ref),400,128));
title(['Xs' num2str(ref)]);

decompo_pulses = {};
% 5-25Hz带通滤波
fsampu = 2000;
[B,A] = butter(4,[5,25]/fsampu*2);
T_filter = filtfilt(B,A,T);
% z-score标准化
T_mean = mean(T_filter);
T_std = std(T_filter);
T_norm = (T_filter-T_mean)./T_std;
for i = 1:size(T_norm,2)
    [~, locs] = findpeaks(-T_norm(:,i),'MinPeakDistance',50,'MinPeakHeight',0.5);
    decompo_pulses{end+1} = locs';
end

figure;
subplot(2,1,1);
plot(-T_norm(:,decomp))
title(['T' num2str(decomp)]);
subplot(2,1,2);
plot(Xt(:,ref))
title(['Xt' num2str(ref)]);

plotDecomps({decompo_pulses{decomp}, ipulses{ref}}, [], 2000, 0, 0, []);
yticks(1:2);
yticklabels({['T' num2str(decomp)], ['Xt' num2str(ref)]})

%% 统计ipulses的CoV
figure;
pulseDiff = {};
for i = 1:length(ipulses)
    pulseDiff{i} = diff(ipulses{i})/2;
    subplot(2,5,i);
    plot(pulseDiff{i});
    ylabel('ms');
    CoV(i) = std(pulseDiff{i}) / mean(pulseDiff{i});
end