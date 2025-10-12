addpath('Func');
clear all;
%% 导入MUAPT数据
iPulses_num = 1;% 选择第几个MU的MUAPT
iEMG_filename = 'IEMG_UU_8000_Trial2_1114_iPulses';
Time_filename = 'TimeResponse2';% 用于卷积的单周期响应曲线
Time_Compo_filename = 'Time_component2';% 最终保存的MU时间分量的名称

fsampu = 2000; % 采样率
data = importdata(['./Data/iEMG/' iEMG_filename '.eaf']); % 读取eaf文件
muNum = max(data.data(:,2)); % MU的个数
iPulses = {};
for mu = 1:muNum
    % iPulses就是这个eaf文件里分解得到的spike train，每个cell表示一个MU，里面的数字是该MU每次放电的时刻
    iPulses{mu} = round(data.data(find(data.data(:,2)==mu),1)'*fsampu); 
end

plotDecomps(iPulses,[],2000,1,1,[]);%~ 绘制MUAPT的图

pulses = zeros(iPulses{1,iPulses_num}(end),1);
pulses(iPulses{1,iPulses_num}) = 1;

figure
stem(pulses);
ylim([0 2]);


%% 绘制手动标定的单个周期运动单位机械响应（速度）

TimeResponse = xlsread(['./Data/simulation/MU_time_response/' Time_filename '.xlsx']);

x = TimeResponse(:,1); 
y = TimeResponse(:,2);
y_norm = normalize(y);

% figure;
% scatter(x,y_norm,'filled')

%~ 插值
% 原始数据的采样率
originalSamplingRate = 1;
originalTimePoints = x;
originalData = y_norm; 

% 目标采样率（例如：重采样为原来的targetSamplingRate倍）
targetSamplingRate = 4;

% 计算重采样后的时间点或数据点位置
originalTimePoints = (0:length(originalData)-1) * originalSamplingRate;
targetTimePoints = (0:1/targetSamplingRate:(length(originalData)-1)/originalSamplingRate);

targetTimePoints_norm = normalization(targetTimePoints,0,200);

% 使用interp1函数进行插值
resampledData = interp1(originalTimePoints, originalData, targetTimePoints, 'linear');
resampledData_norm = normalization(resampledData,-1,1);

figure;
scatter(targetTimePoints_norm,resampledData_norm,'filled')

%% 与MU的spike train做卷积
MU_pulses = pulses(1:2*fsampu);
MU_conv = conv(MU_pulses,resampledData);

% 添加高斯噪声
snr = 20; % 信噪比为~dB
MU_noisy_conv = awgn(MU_conv, snr, 'measured');

figure
plot(MU_noisy_conv)

%% 存储卷积得到的时间分量
% save(['./Data/simulation/MU_time_response/' Time_Compo_filename '.mat'],'MU_noisy_conv')

%% 拟合采样点，这种方法不好，还是对称的曲线
% 定义要进行拟合的多项式阶数
polyOrder = 5;

% 使用polyfit函数拟合多项式曲线
coefficients = polyfit(targetTimePoints_norm, resampledData_norm, polyOrder);

% 创建拟合曲线的x值范围
xFit = 1:0.1:200;

% 使用polyval函数计算拟合曲线的y值
yFit = polyval(coefficients, xFit);

% 绘制原始数据散点图和拟合曲线
figure;
scatter(targetTimePoints_norm, resampledData_norm, 'filled', 'DisplayName', 'Data Points');
hold on;
plot(xFit, yFit, 'r', 'LineWidth', 2, 'DisplayName', 'Fitted Curve');
hold off;
