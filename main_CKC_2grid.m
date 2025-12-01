%% 用CKC反解EMG
clear; clc; close all;
addpath('./Func');

% EMG采样频率
fsamp = 2048;
% US采样频率
fsampu = 2000;

emgFile = './Data/EMG/25-07-04/M1L1T1.mat';
try
    load(emgFile);
catch ME
    rethrow(ME);
end

%% 根据不同数据类型进行修改
newdata = Data';
% newdata(newdata > 32768) = newdata(newdata > 32768) - 2^16;
trigger = newdata(end-1, :);
[~, edges] = maxk(diff(trigger), 2);
edges = sort(edges);
sEMG = newdata(1:end-2, :);

% 绘制trigger图像
figure;
plot(trigger);
hold on;
plot(edges(1), trigger(edges(1)), 'ro');
plot(edges(2), trigger(edges(2)), 'ro');

%% CKC反解EMG
pulsesRef = {};

for ni = 1:2
    data = sEMG(64*(ni-1)+1:64*ni,:);

    % 解码参数设置
    decoderParameters.fsamp = fsamp;
    decoderParameters.TimeDifference = 0;
    decoderParameters.SpatialDifference = 0;
    decoderParameters.ElectrodeType = 18; %13-5*13; 18-mouvi8*8%%%需要注意电极片位置
    decoderParameters.BandpassFilter = 1; %10-500Hz带通滤波
    decoderParameters.LineFilter = 1; %50Hz梳状滤波
    decoderParameters.ChannelFilter = 1; %去除不好的电极channel
    decoderParameters.extendingFactor = 10; %论文里是10c
    decoderParameters.costFcn = 3;
    decoderParameters.iterationNumW = 45; %45
    decoderParameters.iterationNumMU = 30;

    [decompData,dataarray,datafilt,prohibitInd,decompChannelInd] = PreProcess4GUI_v2(data,decoderParameters);
    decomps{ni} = IPTExtraction_gCKC4GUI_v3(decompData,decoderParameters); % classic gradient CKC
    decomps{ni}.decompChannelInd = decompChannelInd;

    for mu = 1:length(decomps{ni}.MUPulses)
        tmp = decomps{ni}.MUPulses{mu};
        tmp = tmp(tmp >= edges(1) & tmp <= edges(2));
        tmp = tmp - edges(1);
        tmp = round(tmp/fsamp*fsampu);
        pulsesRef{end+1} = tmp;
    end

    plotDecomps(decomps{ni}.MUPulses, [], fsamp, 0, 0, []);
end
save('./Data/experiment/25-07-04/M1L1T1_decomps.mat', 'decomps', 'pulsesRef');

%% 绘制IPT
for ni = 1:2
    IPTs = decomps{ni}.IPTs;
    PNRs = decomps{ni}.PNRs;
    MUPulses = decomps{ni}.MUPulses;
    for i = 1:size(IPTs, 1)
        figure;
        plot(IPTs(i, :));
        hold on;
        plot(MUPulses{i}, IPTs(i, MUPulses{i}), 'ro');
        xline([edges(1), edges(2)], '--');
        set(gcf,'unit','normalized','position',[0.05,0.1,0.9,0.6]);
        title(num2str(PNRs(i)));
    end
end
