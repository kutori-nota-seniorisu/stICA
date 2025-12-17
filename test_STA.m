clear; clc; close all;
addpath('./Func');

% EMG采样率
fsamp = 2048;
% US采样率
fsampu = 1000;
% 提取MUTA的阈值
thresh = 0.65;
% US-STA所用窗口100帧
M = 100;
% EMG文件路径
emgFile = './Data/experiment/ICdata/R17/R17_decompsRef.mat';
% US文件路径
tviFile = './Data/experiment/ICdata/R17/v_2d_all.mat';

%% 导入US数据
load(tviFile);
% 数据预处理
disp('开始数据预处理');
tic;
TVIData = cat(3, zeros(119, 128, 2), v_2d_all);
% filter the TVI data
TVIDataFilter = TVIData;
% 时间5-100Hz带通滤波
[Be2, Ae2] = butter(4, [5, 100]/fsampu*2);
for r = 1:size(TVIDataFilter, 1)
    for c = 1:size(TVIDataFilter, 2)
        tmp = squeeze(TVIDataFilter(r, c, :));
        tmp = filtfilt(Be2, Ae2, tmp);
        TVIDataFilter(r, c, :) = tmp;
    end
end
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
    % US-STA窗口为100帧，故去除无法提取到完整100帧图像的pulse
    pulsesRef{mu}(pulsesRef{mu}<=(0+50)) = [];
    pulsesRef{mu}(pulsesRef{mu}>=(size(TVIDataFilter, 3)-50)) = [];
end

%% STA处理
tic;
twitchFrames = cell(0);
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

save('./Data/experiment/ICdata/R17/STA_DecompResult.mat', 'twitchCurves', 'activation', 'twitchFrames');