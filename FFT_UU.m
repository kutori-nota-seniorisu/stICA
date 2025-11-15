% 傅里叶变换绘图
% clear; close all; clc;
%% 生成待处理信号
fs = 2000; %Hz 采样频率
Ts = 1/fs; %采样周期 1ms
t = (0:L-1)*Ts; %时间分辨率用于生成待处理信号,总共2s的信号数据
% 待处理信号的频率
f1 = 50;  f2 = 100;
f3 = 200; f4 = 400;
x1 = 2*0.5*sin(2*pi*f1*t);
x2 = 2*0.2*sin(2*pi*f2*t);
x3 = 2*0.3*sin(2*pi*f3*t);
x4 = 2*0.6*sin(2*pi*f4*t);
x_source = x1 + x2 + x3 + x4; %待处理信号由四个分量组成
% x_source = T_filter;
L  = length(x_source); %序列长度
%% 平时使用的傅里叶变换---1
Y_1 = fft(x_source); % 傅里叶变换后为复数
% 计算双边谱P2;然后根据P2和偶值信号长度N计算单侧频谱P1。
P2_1 = abs(Y_1/L);
P1_1 = P2_1(1:L/2+1);
P1_1(2:end-1) = 2*P1_1(2:end-1); %！！双边变单边,所以乘以2！！！！！！！！！
% P1 = 2*abs(Y(1:N/2+1))/N; %！！双边变单边,所以乘以2！！！！！！！！！
Y_angle_1 = angle(Y_1(1:L/2+1));
delta_f_1 = fs/L; %频率分辨率 0.5
f_1 = (0:1:L/2)*(delta_f_1); %将横坐标转化，显示为频率最大为500 (0, 0.5, ..., 500)共1001个数
figure;
subplot(3,1,1);plot(t,x_source);title('原信号');grid on;
subplot(3,1,2);plot(f_1,P1_1);grid on;title('Single-Sided Amplitude Spectrum of S(t)');
xlabel('f (Hz)');ylabel('|P1(f)|');
subplot(3,1,3);plot(f_1,Y_angle_1);title('原信号频谱相位特性');grid on;

%%
% 演示频率混叠
Fs = 1000;
t = 0:1/Fs:1-1/Fs;

% 两个不同频率但采样后看起来相同的信号
f1 = 100;    % 100 Hz
f2 = 900;    % 900 Hz

x1 = sin(2*pi*f1*t);
x2 = sin(2*pi*f2*t);

figure;
subplot(2,1,1);
plot(t(1:100), x1(1:100), 'b-', 'LineWidth', 2);
hold on;
plot(t(1:100), x2(1:100), 'r--', 'LineWidth', 2);
title('时域信号：100Hz和900Hz在1000Hz采样下看起来相同');
legend('100 Hz', '900 Hz');
xlabel('时间 (s)');
ylabel('幅度');

% 频率分析
Y1 = fft(x1);
Y2 = fft(x2);
f = Fs*(0:length(t)/2)/length(t);

P1 = abs(Y1/length(t));
P1 = P1(1:length(t)/2+1);
P1(2:end-1) = 2*P1(2:end-1);

P2 = abs(Y2/length(t));
P2 = P2(1:length(t)/2+1);
P2(2:end-1) = 2*P2(2:end-1);

subplot(2,1,2);
plot(f, P1, 'b-', 'LineWidth', 2);
hold on;
plot(f, P2, 'r--', 'LineWidth', 2);
title('频谱显示两者在100Hz处有相同的峰值');
xlabel('频率 (Hz)');
ylabel('幅度');
legend('100 Hz信号', '900 Hz信号');

%%
% 展示FFT的完整对称性
Fs = 1000;
N = 64;  % 使用较小点数便于观察
t = (0:N-1)/Fs;
x = sin(2*pi*50*t) + 0.5*sin(2*pi*150*t);

Y = fft(x);
f_full = Fs*(0:N-1)/N;      % 完整频率向量 (-Fs/2 到 Fs/2)
f_shifted = f_full - Fs/2;  % 将0频率移到中心

% 重新排列FFT结果（fftshift）
Y_shifted = fftshift(Y);

figure;
subplot(3,1,1);
plot(t, x);
title('时域信号');
xlabel('时间 (s)');
ylabel('幅度');

subplot(3,1,2);
stem(f_full, abs(Y), 'filled');
title('FFT结果 - 原始顺序 (0 到 Fs)');
xlabel('频率 (Hz)');
ylabel('|Y(f)|');

subplot(3,1,3);
stem(f_shifted, abs(Y_shifted), 'filled');
title('FFT结果 - 零频居中 (-Fs/2 到 Fs/2)');
xlabel('频率 (Hz)');
ylabel('|Y(f)|');
xlim([-Fs/2, Fs/2]);