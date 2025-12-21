% STA算法实现
clear; clc; close all;
addpath('./Func');

% EMG采样率
fsamp = 2048;
% US采样率
fsampu = 2000;
% 提取MUTA的阈值
thresh = 0.65;
% US-STA所用窗口100ms
M = 100/1000*fsampu;

motion = 2;
emgFile = ['./Data/experiment/25-07-04/M' num2str(motion) 'L1T1_decompsRef.mat'];
tviFile = ['./Data/experiment/25-07-04/TVIData_15000_S_wrl_M' num2str(motion) '_level1_trial1_Single_25-07-04.mat'];

% emgFile = './Data/experiment/ICdata/R17/R17_decompsRef.mat';
% tviFile = './Data/experiment/ICdata/R17/v_2d_all.mat';


%% 导入US数据
load(tviFile);
% 数据预处理
disp('开始数据预处理');
tic;

% TVIData = cat(3, zeros(119, 128, 2), v_2d_all);
TVIData = cat(3, zeros(395, 128, 20), TVIData);

% filter the TVI data
% TVIDataFilter = TVIData;
TVIDataFilter = TVIData(:, :, 2001:end);

% 轴向0.5MHz低通滤波
[Be1, Ae1] = butter(4, 0.5/(7.7*4)*2, 'low');
for i = 1:size(TVIDataFilter, 3)
    tmp = TVIDataFilter(:, :, i);
    tmp = filtfilt(Be1, Ae1, tmp);
    TVIDataFilter(:, :, i) = tmp;
end

% 时间5-100Hz带通滤波
[Be2, Ae2] = butter(4, [5, 100]/fsampu*2);
for r = 1:size(TVIDataFilter, 1)
    for c = 1:size(TVIDataFilter, 2)
        tmp = squeeze(TVIDataFilter(r, c, :));
        tmp = filtfilt(Be2, Ae2, tmp);
        TVIDataFilter(r, c, :) = tmp;
    end
end

% 对每一列降采样
for i = 1:size(TVIDataFilter, 3)
    tmp = TVIDataFilter(:, :, i);
    tmp = resample(tmp, 128, 395);
    TVITmp(:, :, i) = tmp;
end
TVIDataFilter = TVITmp;
clear TVITmp;

toc;
disp(['数据预处理用时' num2str(toc)]);

%% 导入EMG数据
load(emgFile);
pulsesAll = decompsRef.Pulses;
muapsAll = decompsRef.MUAPs;
numMU = length(pulsesAll);
% MU spike train的偏移
shiftAP = zeros(1, numMU);
for mu = 1:numMU
    tmpArray = plotArrayPotential(muapsAll{mu}, 1, 0);
    tmpArrayDiff2 = cell(0);
    for nr = 1:size(tmpArray, 1)
        for nc = 1:size(tmpArray, 2)
            if ~isempty(tmpArray{nr, nc})
                % 对每个通道上的MUAP进行二次差分
                tmpArrayDiff2{nr, nc} = diff(diff(tmpArray{nr, nc}));
            end
        end
    end
    [~, ~, ~, pos] = plotArrayPotential(tmpArrayDiff2, 1, 0);
    tmptmp = tmpArrayDiff2{pos(1), pos(2)};
    tmpInd = find(abs(tmptmp)>5*std(tmptmp));
    if ~isempty(tmpInd)
        shiftAP(mu) = tmpInd(1) - 64;
    end
end
pulsesRef = cell(1, numMU);
for mu = 1:numMU
    % 转换成超声样本点时刻
    pulsesRef{mu} = round((pulsesAll{mu}+shiftAP(mu))/fsamp*fsampu);
    % 如果TVI去除了一部分，那就需要执行下面这一步。
    pulsesRef{mu} = pulsesRef{mu}-fsampu;
    % US-STA窗口为100帧，故去除无法提取到完整100帧图像的pulse
    pulsesRef{mu}(pulsesRef{mu}<=(0+M/2)) = [];
    pulsesRef{mu}(pulsesRef{mu}>=(size(TVIDataFilter, 3)-M/2)) = [];
end

%% STA处理
tic;
twitchFrames = cell(0);
varFrames = cell(0);
activation = cell(0);
twitchCurves = cell(0);
for mu = 1:numMU
    tmpPulse = pulsesRef{mu};

    frameSTA = zeros(0);
    varSTA = zeros(0);
    for n = -M/2+1:M/2
        tmpInd = tmpPulse + n;
        tmpTVI = TVIDataFilter(:, :, tmpInd);
        frameSTA(:, :, n+M/2) = mean(tmpTVI, 3);
        varSTA(:, :, n+M/2) = var(tmpTVI, 0, 3);
    end
    
    twitchFrames{mu} = frameSTA;
    varFrames{mu} = varSTA;
    directions = sum(frameSTA(:,:,1:M/2)./varSTA(:,:,1:M/2),3)-sum(frameSTA(:,:,M/2+1:end)./varSTA(:,:,M/2+1:end),3);
    directionsSign = sign(directions);
    MUIntensity = sum(frameSTA.^2./varSTA, 3).*-directionsSign;
    activation{mu} = MUIntensity;

    ax=figure;
    % actArea{1}是正激活度大于阈值，actArea{2}是负激活度小于阈值
    [~, actArea] = plotActivationZone(muapsAll{mu}, activation{mu}, [], [], thresh);
    twitchCurves{mu} = twitchCurve(twitchFrames{mu}, actArea);
    saveas(ax, num2str(mu), 'png');
end
toc;
disp(['STA处理用时' num2str(toc)]);

% save('./Data/experiment/ICdata/R17/STA_DecompResult.mat', 'twitchCurves', 'varFrames', 'activation', 'twitchFrames');
save(['./Data/experiment/25-07-04/M' num2str(motion) 'L1T1_STA_DecompResult.mat'], 'twitchCurves', 'varFrames', 'activation', 'twitchFrames');