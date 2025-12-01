clear; clc; close all;

%% 对ICData的结果进行统计
EMGMU = [1,2,2,4,4,4,4,5,5,7,6,6];
recordLabel = [3,4,5,7,10,11,12,14,15,16,17,18];
countAll = 0;
for i = 1:length(recordLabel)
    sub = recordLabel(i);
    load(['./Data/experiment/ICdata/R' num2str(sub) '/USCBSS_decomp_result.mat']);
    countMU(i) = length(decompoMUFiltered);
    countAll = countAll + countMU(i);
end
figure;
plot(countMU);
hold on;
plot(EMGMU);
xticks(1:length(recordLabel));
xticklabels(recordLabel);
ylabel('the num of MU');
xlabel('record')
title('各trial分解得到的MU个数');
legend('US', 'EMG');
set(gcf,'unit','normalized','position',[0.15,0.1,0.7,0.6]);


TPR = [];
figure;
for i = 1:length(recordLabel)
    sub = recordLabel(i);
    load(['./Data/experiment/ICdata/R' num2str(sub) '/USCBSS_decomp_result.mat']);
    if isempty(matchResultRaw)
        continue;
    end
    boxplot(matchResultRaw.("match ratio"), 'Positions', i, 'Colors', 'r');
    TPR = [TPR; matchResultRaw.("match ratio")];
    hold on;
end
xticks(1:length(recordLabel));
xticklabels(recordLabel);
ylabel('true positive ratio (%)');
xlabel('record')
title('真阳性率')
set(gcf,'unit','normalized','position',[0.15,0.1,0.7,0.6]);
mean(TPR)
std(TPR)
median(TPR)

RoA = [];
figure;
for i = 1:length(recordLabel)
    sub = recordLabel(i);
    load(['./Data/experiment/ICdata/R' num2str(sub) '/USCBSS_decomp_result.mat']);
    if isempty(matchResultRaw)
        continue;
    end
    boxplot(matchResultRaw.RoA*100, 'Positions', i, 'Colors', 'r');
    RoA = [RoA; matchResultRaw.RoA*100];
    hold on;
end
xticks(1:length(recordLabel));
xticklabels(recordLabel);
ylabel('RoA (%)');
xlabel('record')
title('一致率RoA')
set(gcf,'unit','normalized','position',[0.15,0.1,0.7,0.6]);
mean(RoA)
std(RoA)
median(RoA)

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
mu = 43;
twitch = decompoMUFiltered.Twitch(:, mu);
twitch = twitch * (1 - 2 * (abs(max(twitch)) <= abs(min(twitch))));
source = decompoMUFiltered.Source(:, mu);
fsampu = 2000;
L = 30000;

Y1 = fft(twitch);
P21 = abs(Y1/L);
P11 = P21(1:L/2+1);
P11(2:end-1) = 2*P11(2:end-1);

Y2 = fft(source);
P22 = abs(Y2/L);
P12 = P22(1:L/2+1);
P12(2:end-1) = 2*P12(2:end-1);

f = (0:1:L/2)*fsampu/L;

figure;
subplot(4,1,1);
plot(twitch);
title('twitch')
subplot(4,1,2);
plot(f, P12);
subplot(4,1,3);
plot(source);
title('source')
subplot(4,1,4);
plot(f, P11);
% xlim([0,100]);
%%
mu = 5;
fsampu = 2000;
L = 30000;

source = decompoMUFiltered.Source(:, mu);
Y1 = fft(source);
P21 = abs(Y1/L);
P11 = P21(1:L/2+1);
P11(2:end-1) = 2*P11(2:end-1);

[Be, Ae] = butter(4, [5,50]/fsampu*2);
sourceF = filtfilt(Be, Ae, source);
Y2 = fft(sourceF);
P22 = abs(Y2/L);
P12 = P22(1:L/2+1);
P12(2:end-1) = 2*P12(2:end-1);

f = (0:1:L/2)*fsampu/L;

figure;
subplot(4,1,1);
plot(source);
title('source')
subplot(4,1,2);
plot(f, P11);

subplot(4,1,3);
plot(sourceF);
title('source after filter')
subplot(4,1,4);
plot(f, P12);
% xlim([0,100]);
%%
TVIData = cat(3, zeros(395, 128, 20), TVIData);
% TVISum = squeeze(sum(sum(TVIData)));
TVISum = squeeze(TVIData(200,64,:));
fsampu = 2000;
L = 30000;
f = (0:1:L/2)*fsampu/L;
Y = fft(TVISum);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure;
subplot(2,1,1);
plot(TVISum);
subplot(2,1,2);
plot(f, P1);

%%
fff = zeros(128,128);
Row = 8; Col = 8;
dRow = 4; dCol = 4;
for mu=75%[4,5,11,20,21,22,54,55]
    r = decompoMUFiltered.Row(mu);
    c = decompoMUFiltered.Col(mu);
    winRow = (1:Row)+(r-1)*dRow;
    winCol = (1:Col)+(c-1)*dCol;
    fff(winRow,winCol) = fff(winRow,winCol)+1;
end
figure;
imagesc(fff);
colorbar

%% 取round会造成最大约0.5个样本点的误差，对应于0.25ms的误差。
t = 0:1e-4:2;
y1 = t*2000/2048;
y2 = round(t*2000/2048);
figure;
plot(t,y1);
hold on;
plot(t,y2);
legend('y1', 'y2');
figure;
plot(t, abs(y1-y2));

%% ipulse和6PTM模型卷积，作为仿真数据
% 模型ref: Variability of successive contractions subtracted from unfused tetanus of fast and slow motor units
addpath('./Func');

fsampu = 2000; % 采样率

datasets_num = '1';
time_datasets = ['TimeCompoDatasets' datasets_num];% 时间分量数据集所在的文件夹
load(['./Data/simulation/MU_time_response/' time_datasets '/ipulses.mat']);

muNum = length(ipulses);

for mu = 1:muNum
    % 脉冲串长度
    L = length(ipulses{mu});

    % 6PTM模型，单位ms
    % 考虑收缩曲线的变异性，对于每一个放电时刻都要生成一条收缩曲线
    Tlead = 0;
    Thc = unifrnd(20,25,L,1);
    Tc = unifrnd(50,75,L,1);
    Thr = unifrnd(100,130,L,1);
    Ttw = unifrnd(300,350,L,1);
    Dmax = unifrnd(40,70,L,1);

    MU_conv_all = [];

    for i = 1:L
        t = 0:Ttw(i);

        e = exp(1);
        c1 = log(2)*Tc(i) / (Thc(i) - Tc(i) + Tc(i)*log(Tc(i)/Thc(i)));
        c2 = log(2)*Tc(i) / (Thr(i) - Tc(i) + Tc(i)*log(Tc(i)/Thr(i)));

        P1 = (t - Tlead) / Tc(i);
        P2 = 1 + exp(2*e*(P1-1));
        P3 = (t - 0.5*(Ttw(i)+Thr(i))) / (Ttw(i)-Thr(i));

        f1 = Dmax(i) * (P1.^c1 .* exp(c1-c1*P1) + (P2-1) .* P1.^c2 .* exp(c2-c2*P1));
        f2 = P2 + P2.*exp(4*e*P3);

        F = f1 ./ f2;
        
        pulse = zeros(1,2*fsampu);
        pulse(ipulses{mu}(i)) = 1;

        tmp = conv(pulse, F);
        MU_conv_all(i,:) = tmp(1:4000);
    end

    MU_conv = sum(MU_conv_all);
    % 差分，得到速度曲线
    MU_conv_diff = diff(MU_conv);
    MU_conv_diff(end+1) = 0;

    % 添加高斯噪声
    snr = 20; % 信噪比为~dB
    MU_noisy_conv = awgn(MU_conv_diff, snr, 'measured');

    figure;
    subplot(3,1,1);
    plot(MU_conv);
    title(['MU' num2str(mu) ' twitch']);
    subplot(3,1,2);
    plot(MU_conv_diff);
    title(['MU' num2str(mu) ' twitch(diff)']);
    subplot(3,1,3);
    plot(MU_noisy_conv);
    title(['MU' num2str(mu) ' twitch(diff,noisy)']);

    % save(['F:/EEEMG/stICA/Data/simulation/MU_time_response/' time_datasets '/Time_component' num2str(mu) '.mat'],'MU_noisy_conv')
end


%%
for ni = 1:2
    IPTs = decomps{ni}.IPTs;
    PNRs = decomps{ni}.PNRs;
    Pulses = decomps{ni}.MUPulses;
    for i = 1:size(IPTs, 1)
        figure;
        % subplot(2,1,1);
        plot(IPTs(i, :));
        hold on;
        plot(Pulses{i}, IPTs(i, Pulses{i}), 'ro');
        % hold on;
        % plot(edges(2), IPTs(i, edges(2)), 'ro');
        % plot(edges(1), IPTs(i, edges(1)), 'ro');
        set(gcf,'unit','normalized','position',[0.05,0.1,0.9,0.6]);
        title(['ni=' num2str(ni) ',mu=' num2str(i) ',PNR=' num2str(PNRs(i))]);
        % subplot(2,1,2);
        % scatter(Pulses{i}, zeros(length(Pulses{i})), 1000, "black", '|');
        % xlim([edges(2), edges(1)])
        % ylim([0, 1e-3])
    end
end

%% 绘制TVI滤波前后的时域与频域图像
s=squeeze(TVIData(59,70,:));
fsampu = 1000;
L=length(s);
Y=fft(s);
P2=abs(Y/L);
P1=P2(1:L/2+1);
P1(2:end-1)=2*P1(2:end-1);
f=(0:L/2)*fsampu/L;
figure;
subplot(2,1,1);
plot(s);
xticklabels(0:5:30);
xlabel('Time (s)');
ylabel('Amplitude');
title('TVI (滤波前)');
subplot(2,1,2);
plot(f,P1);
xlabel('frequency (Hz)')
ylabel('amplitude');
title('幅度谱');

s2=squeeze(TVIDataFilter(59,70,:));
fsampu = 1000;
L2=length(s2);
Y2=fft(s2);
P22=abs(Y2/L2);
P12=P22(1:L2/2+1);
P12(2:end-1)=2*P12(2:end-1);
f2=(0:L2/2)*fsampu/L2;
figure;
subplot(2,1,1);
plot(s2);
xticklabels(0:5:30);
xlabel('Time (s)');
ylabel('Amplitude');
title('TVI (滤波后)');
subplot(2,1,2);
plot(f2,P12);
xlabel('frequency (Hz)')
ylabel('amplitude');
title('幅度谱');

%% 统计带通滤波后USCBSS的结果，两个指标：MAD和ER
fsampu = 1000;
MADAll = [];
ERAll = [];
[rn, cn] = size(DecompoResults.sources);
for r = 1:rn
    for c = 1:cn
        tmpPulses = DecompoResults.decompo_pulses{r, c};
        tmpSources = DecompoResults.sources{r, c};
        tmpCoV = DecompoResults.CoV{r, c};

        for mu = 1:length(tmpPulses)
            % 计算MAD，单位为ms
            MAD = mad(diff(tmpPulses{mu}/fsampu*1000));
            % 计算能量占比
            [psd, freq] = pwelch(tmpSources(:, mu), [], [], [], fsampu);
            % frequency band of interest
            freqbandOI = [6, 14];
            idxBand = freq >= freqbandOI(1) & freq <= freqbandOI(2);
            idxTotal = freq > 0;
            energyInBand = trapz(freq(idxBand), psd(idxBand));
            energyTotal = trapz(freq(idxTotal), psd(idxTotal));
            energyRatio = energyInBand / energyTotal * 100;
            MADAll(end+1) = MAD;
            ERAll(end+1) = energyRatio;
        end
    end
end

figure;
yyaxis left;
boxplot(MADAll, 'Positions', 1);
hold on;
yyaxis right;
boxplot(ERAll, 'Positions', 2);

figure;
histogram(MADAll);
figure;
histogram(ERAll);

%% 是否滤波对结果的影响，效果绘制
fsampu = 1000;
% 未滤波的结果绘制
load('./Data/experiment/ICdata/R16/USCBSS_DecompResult.mat');
s1 = decompoMUFiltered.Source(:, 1);
t1 = decompoMUFiltered.Twitch(:, 1);
p1 = decompoMUFiltered.Pulse{1};
figure;
subplot(2,1,1);
plot(t1);
xlim([4000,8000]);
xticks(4000:1000:8000);
xticklabels(4:1:8);
xlabel('t/s'); ylabel('amplitude');
title('MU 1 twitch curve');
subplot(2,1,2);
plot(s1);
hold on;
plot(p1, s1(p1), 'ro');
xlim([4000,8000]);
xticks(4000:1000:8000);
xticklabels(4:1:8);
xlabel('t/s'); ylabel('amplitude');
title('MU 1 source');

% 滤波结果绘制
load('./Data/experiment/ICdata/R16/USCBSS_compo25_filter.mat');
tmpPulses = DecompoResults.decompo_pulses{3, 11};
tmpSources = DecompoResults.sources{3, 11};
tmpSourcesRaw = DecompoResults.sourceFirst{3, 11};
s2 = tmpSources(:, 1);
t2 = tmpSourcesRaw(:, 1);
p2 = tmpPulses{1};
figure;
subplot(2,1,1);
plot(t2);
xlim([4000,8000]);
xticks(4000:1000:8000);
xticklabels(4:1:8);
xlabel('t/s'); ylabel('amplitude');
title('twitch curve');
subplot(2,1,2);
plot(s2);
hold on;
plot(p2, s2(p2), 'ro');
xlim([4000,8000]);
xticks(4000:1000:8000);
xticklabels(4:1:8);
xlabel('t/s'); ylabel('amplitude');
title('source');

%%
plotDecomps(decompoMUFiltered.Pulse, [], 1000, 0, 0, []);
title('MU from US')

plotDecomps(pulsesRef(1:8), [], 1000, 0, 0, []);
title('MU from sEMG')
