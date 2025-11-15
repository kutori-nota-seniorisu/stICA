%% 用CKC反解EMG，两块5*13的电极片
clear; clc; close all;
addpath('./Func');

% EMG采样频率
fsamp = 2048;
% US采样频率
fsampu = 1000;

for sub = 3%[3,4,5,7,10,11,12,14,15,16,17,18]
    
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
    sEMG = newdata(1:end-2, :);

    % 绘制trigger图像
    figure;
    plot(trigger);
    hold on;
    plot(edges(1), trigger(edges(1)), 'ro');
    plot(edges(2), trigger(edges(2)), 'ro');
    
    pulsesRef = {};

    % CKC反解EMG
    % 两块电极片
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
        decomps{ni} = IPTExtraction_gCKC4GUI_v3(decompData,decoderParameters); % classic gradient CKC
        decomps{ni}.decompChannelInd=decompChannelInd;

        pulsesAll = decomps{ni}.MUPulses;
        for mu = 1:length(pulsesAll)
            tmp = pulsesAll{mu};
            tmp = tmp(tmp >= edges(1) & tmp <= edges(2));
            if length(tmp) < 10
                continue;
            end
            tmp = tmp - edges(1);
            tmp = round(tmp/fsamp*fsampu);
            pulsesRef{end+1} = tmp;
        end

        plotDecomps(decomps{ni}.MUPulses, [], fsamp, 0, 0, []);
        % xlim([edges(2), edges(2)+20480]/fsamp);
    end
    % savepath = ['./Data/experiment/ICdata/R' num2str(sub) '/'];
    % save([savepath 'R' num2str(sub) '_decomps.mat'], 'decomps');
    % save([savepath 'pulsesRef.mat'], 'pulsesRef');
end
