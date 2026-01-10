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
