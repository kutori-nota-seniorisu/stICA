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
for ni = 1:2
    IPTs = decomps{ni}.IPTs;
    PNRs = decomps{ni}.PNRs;
    for i = 1:size(IPTs, 1)
        figure;
        plot(IPTs(i, :));
        hold on;
        plot(edges(2), IPTs(i, edges(2)), 'ro');
        plot(edges(1), IPTs(i, edges(1)), 'ro');
        set(gcf,'unit','normalized','position',[0.05,0.1,0.9,0.6]);
        title(num2str(PNRs(i)));
        % xlim([edges(2), edges(1)])
        % ylim([0, 1e-3])
    end
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
