%% 用US结果进行STA
clear; clc; close all;
addpath('./Func');

% for sub = [16, 17]
sub = 16;
disp(['Sub=' num2str(sub)]);
load(['./Data/experiment/ICdata/R' num2str(sub) '/USCBSS_DecompResult_MPD50.mat']);

numMU = length(decompoMUFiltered.MU);

%% twitch区域，STA
win = [-100, 100];
fsampu = 2000;

% 导入超声速度图像
% tviFile = ['./Data/experiment/ICdata/R' num2str(sub) '/v_2d_all.mat'];
tviFile = './Data/experiment/25-07-04/TVIData_15000_S_wrl_M1_level1_trial1_Single_25-07-04.mat';
load(tviFile);

TVIData = cat(3, zeros(395, 128, 20), TVIData);

% filter the TVI data
TVIDataFilter = TVIData;

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

[M, N, L] = size(TVIDataFilter);
% 窗口大小
% Row = 10; Col = 10;
Row = 8; Col = 8;
% 窗口移动距离
% dRow = 5; dCol = 5;
dRow = 4; dCol = 4;

for mu = 1:numMU
    tmpPulses = decompoMUFiltered.Pulse{mu};
    if isempty(tmpPulses)
        continue;
    end
    frameSTA = zeros(0);
    varSTA = zeros(0);
    for n = win(1)+1:win(2)
        tmpInd = tmpPulses + n;
        tmpInd(tmpInd <= 0) = [];
        tmpInd(tmpInd >= 30000) = [];
        tmpTVI = TVIDataFilter(:, :, tmpInd);
        % 存储STA图像
        frameSTA(:, :, n-win(1)) = mean(tmpTVI, 3);
        % 存储方差图像
        varSTA(:, :, n-win(1)) = var(tmpTVI, 0, 3);
    end

    r = decompoMUFiltered.Row(mu);
    c = decompoMUFiltered.Col(mu);
    winRow = (1:Row)+(r-1)*dRow;
    winCol = (1:Col)+(c-1)*dCol;
    winRow(winRow>M) = [];
    winCol(winCol>N) = [];

    for r = 1:length(winRow)
        for c = 1:length(winCol)
            frameSTAROI{r, c} = squeeze(frameSTA(winRow(r), winCol(c), :));
        end
    end

    figure;
    ax = subplot('Position', [0.05, 0.05, 0.9, 0.9]);
    [~,maxAmp,maxPP,pos]=plotArrayPotential(frameSTAROI, 1, 1, ax);
    title(['MU ' num2str(mu) ' twitch curve']);
end
% end


%% EMG
% EMG采样频率
fsamp = 2048;
winMUAP = [-50, 50];

emgFile = ['./Data/experiment/ICdata/R' num2str(sub) '/R' num2str(sub) '.mat'];
try
    load(emgFile);
catch ME
    rethrow(ME);
end

% ICData
newdata = Data{1, 2};
newdata(newdata > 32768) = newdata(newdata > 32768) - 2^16;
trigger = newdata(end, :);
[~, edges] = maxk(trigger, 2);
edges = sort(edges);
sEMG = newdata(1:end-2, edges(1):edges(2));

% newdata = Data';
% trigger = newdata(end-1, :);
% [~, edges] = findpeaks(trigger, 'MinPeakDistance', 10000, 'MinPeakProminence', 0.3);
% sEMG = newdata(1:end-2, edges(1):edges(2));

sEMGArray = [];
for ni = 1:2
    data = sEMG(64*(ni-1)+1:64*ni,:);

    % 解码参数设置
    decoderParameters.fsamp = fsamp;
    decoderParameters.TimeDifference = 0;
    decoderParameters.SpatialDifference = 0;
    decoderParameters.ElectrodeType = 18;%13-5*13; 18-mouvi8*8%%%需要注意电极片位置
    decoderParameters.BandpassFilter = 1;%10-500Hz带通滤波
    decoderParameters.LineFilter = 1;%50Hz梳状滤波
    decoderParameters.ChannelFilter = 1;%去除不好的电极channel
    decoderParameters.extendingFactor = 10;%论文里是10c
    decoderParameters.costFcn = 3;
    decoderParameters.iterationNumW = 45;%45
    decoderParameters.iterationNumMU = 30;

    [decompData,dataarray,datafilt,prohibitInd,decompChannelInd] = PreProcess4GUI_v2(data,decoderParameters);

    [rn, cn] = size(dataarray);
    for i = 1:rn
        for j = 1:cn
            if ~isempty(dataarray{i, j})
                sEMGArray(i, j+(cn+2)*(ni-1), :) = dataarray{i, j}(edges(1):edges(2));
            else
                sEMGArray(i, j+(cn+2)*(ni-1), :) = NaN;
            end
        end
    end
end

%% 绘制MUAP
winSize = 128;
lenEMG = diff(edges)+1;

% 存储所有的MUAP阵列
arrayMUAP = cell(1, numMU);

for mu = 1:numMU
    tmpPulses = decompoMUFiltered.Pulse{mu};
    tmpPulses = round(tmpPulses/1000*2048);
    MUAP = [];
    for n = -winSize/2+1:winSize/2
        tmpInd = tmpPulses + n;
        tmpInd(tmpInd <= 0) = [];
        tmpInd(tmpInd >= lenEMG) = [];
        tmpMUAP = sEMGArray(:, :, tmpInd);
        % 存储STA图像
        MUAP(:, :, n+winSize/2) = mean(tmpMUAP, 3);
    end

    % 对MUAP阵列中的NaN进行插值
    MUAP = interpolateMUAP(MUAP);

    % 将数组转成元胞
    for r = 1:size(MUAP, 1)
        for c = 1:size(MUAP, 2)
            tmptmp = squeeze(MUAP(r, c, :));
            if isempty(find(tmptmp))
                arrayMUAP{mu}{r, c} = [];
            else
                arrayMUAP{mu}{r, c} = tmptmp';
            end
        end
    end

    figure;
    ax = subplot('Position', [0.05, 0.05, 0.9, 0.9]);
    plotArrayPotential(arrayMUAP{mu}, 1, 1, ax);
    title(['MU ' num2str(mu) ' MUAP array ']);
    set(gcf,'unit','normalized','position',[0.1,0.5,0.35,0.3]);
end
