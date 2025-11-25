%% 用EMG结果进行STA
clear; clc; close all;
addpath('./Func');

% for sub = [16, 17]
sub = 16;
disp(['Sub=' num2str(sub)]);
load(['./Data/experiment/ICdata/R' num2str(sub) '/pulsesRef.mat']);
numMU = length(pulsesRef);

% EMG采样频率
fsamp = 2048;
winMUAP = [-50, 50];

emgFile = ['./Data/experiment/ICdata/R' num2str(sub) '/R' num2str(sub) '.mat'];
try
    load(emgFile);
catch ME
    rethrow(ME);
end

newdata = Data{1, 2};
newdata(newdata > 32768) = newdata(newdata > 32768) - 2^16;
trigger = newdata(end, :);
[~, edges] = maxk(trigger, 2);
edges = sort(edges);
sEMG = newdata(1:end-2, edges(1):edges(2));
lenEMG = length(sEMG);

for ni = 1:2
    data = sEMG(64*(ni-1)+1:64*ni,:);

    % 解码参数设置
    decoderParameters.fsamp = fsamp;
    decoderParameters.TimeDifference = 0;
    decoderParameters.SpatialDifference = 0;
    decoderParameters.ElectrodeType = 13;%13-5*13; 18-mouvi8*8%%%需要注意电极片位置
    decoderParameters.BandpassFilter = 1;%10-500Hz带通滤波
    decoderParameters.LineFilter = 1;%50Hz梳状滤波
    decoderParameters.ChannelFilter = 1;%去除不好的电极channel
    decoderParameters.extendingFactor = 10;%论文里是10c
    decoderParameters.costFcn = 3;
    decoderParameters.iterationNumW = 45;%45
    decoderParameters.iterationNumMU = 30;

    [decompData,dataarray,datafilt,prohibitInd,decompChannelInd] = PreProcess4GUI_v2(data,decoderParameters);

    for i = 1:5
        for j = 1:13
            if ~isempty(dataarray{i, j})
                sEMGArray(i, j, :) = dataarray{i, j};
            else
                sEMGArray(i, j, :) = NaN;
            end
        end
    end

    for mu = 1:numMU
        tmpPulses = pulsesRef{mu};
        tmpPulses = round(tmpPulses/1000*2048);
        if isempty(tmpPulses)
            continue;
        end
        MUAP = zeros(0);
        for n = winMUAP(1)+1:winMUAP(2)
            tmpInd = tmpPulses + n;
            tmpInd(tmpInd <= 0) = [];
            tmpInd(tmpInd >= lenEMG) = [];
            tmpMUAP = sEMGArray(:, :, tmpInd);
            % 存储STA图像
            MUAP(:, :, n-winMUAP(1)) = mean(tmpMUAP, 3);
        end

        arrayMUAP = {};
        for r = 1:5
            for c = 1:13
                arrayMUAP{r, c} = squeeze(MUAP(r, c, :));
            end
        end

        if ni == 1
            arrayMUAP = rot90(arrayMUAP, 3);
        elseif ni == 2
            arrayMUAP = rot90(arrayMUAP);
        end


        figure;
        ax = subplot('Position', [0.05, 0.05, 0.9, 0.9]);
        [~,maxAmp,maxPP,pos]=plotArrayPotential(arrayMUAP, 1, 1, ax);
        title(['MU ' num2str(mu) ' MUAP array ' num2str(ni)]);
        set(gcf,'unit','normalized','position',[0.4,0.5,0.1,0.3]);
    end
end