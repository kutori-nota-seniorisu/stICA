%% 用CKC反解EMG
clear; clc; close all;
addpath('./Func');

% EMG采样频率
fsamp = 2048;
% US采样频率
fsampu = 2000;

% for sub = [3,4,5,7,10,11,12,14,15,16,17,18]
% emgFile = ['./Data/experiment/ICdata/R' num2str(sub) '/R' num2str(sub) '.mat'];

motion = 2; trial = 2;
emgFile = ['./Data/EMG/25-07-04/M' num2str(motion) 'L1T' num2str(trial) '.mat'];
try
    load(emgFile);
catch ME
    rethrow(ME);
end

%% 根据不同数据类型进行修改
% ICData处理
% newdata = Data{1, 2};
% newdata(newdata > 32768) = newdata(newdata > 32768) - 2^16;
% trigger = newdata(end, :);
% [~, edges] = maxk(trigger, 2);
% edges = sort(edges);
% sEMG = newdata(1:end-2, :);

newdata = Data';
trigger = newdata(end-1, :);
[~, edges] = findpeaks(trigger, 'MinPeakDistance', 10000, 'MinPeakProminence', 0.3);
sEMG = newdata(1:end-2, :);

% 绘制trigger图像
figure;
plot(trigger);
hold on;
plot(edges(1), trigger(edges(1)), 'ro');
plot(edges(2), trigger(edges(2)), 'ro');
drawnow;

%% CKC反解EMG
for ni = 1:2
    data = sEMG(64*(ni-1)+1:64*ni,:);

    % 解码参数设置
    decoderParameters.fsamp = fsamp;
    decoderParameters.TimeDifference = 0;
    decoderParameters.SpatialDifference = 0;
    decoderParameters.ElectrodeType = 18; %13-5*13; 18-mouvi8*8%%%需要注意电极片位置
    decoderParameters.BandpassFilter = 1; %20-500Hz带通滤波
    decoderParameters.LineFilter = 1; %50Hz梳状滤波
    decoderParameters.ChannelFilter = 1; %去除不好的电极channel
    decoderParameters.extendingFactor = 10; %论文里是10c
    decoderParameters.costFcn = 3;
    decoderParameters.iterationNumW = 45; %45
    decoderParameters.iterationNumMU = 30;

    [decompData,dataarray,datafilt,prohibitInd,decompChannelInd] = PreProcess4GUI_v2(data,decoderParameters);
    decomps{ni} = IPTExtraction_gCKC4GUI_v3(decompData,decoderParameters); % classic gradient CKC
    decomps{ni}.decompChannelInd = decompChannelInd;
end

% savepath = ['./Data/experiment/ICdata/R' num2str(sub) '/'];
% save([savepath 'R' num2str(sub) '_decomps.mat'], 'decomps');

save(['./Data/experiment/25-07-04/M' num2str(motion) 'L1T' num2str(trial) '_decomps.mat'], 'decomps');

%% 绘制IPT
for ni = 1:2
    plotDecomps(decomps{ni}.MUPulses, [], fsamp, 0, 0, []);
    IPTs = decomps{ni}.IPTs;
    PNRs = decomps{ni}.PNRs;
    MUPulses = decomps{ni}.MUPulses;

    figure;
    t = tiledlayout('vertical', 'TileSpacing', 'none', 'Padding', 'compact');
    for mu = 1:size(IPTs, 1)
        nexttile;
        plot(IPTs(mu, :));
        hold on;
        plot(MUPulses{mu}, IPTs(mu, MUPulses{mu}), 'ro');
        xline(edges);
        % set(gcf,'unit','normalized','position',[0.05,0.1,0.6,0.3]);
        title(['ni#' num2str(ni) ' MU#' num2str(mu) ' PNR=' num2str(PNRs(mu))]);
    end
end

%% 保存参考脉冲串
pulsesRef = {};
for ni = 1:2
    PNRs = decomps{ni}.PNRs;
    MUPulses = decomps{ni}.MUPulses;
    for mu = 1:length(MUPulses)
        if PNRs(mu) <= 20
            disp(['ni#' num2str(ni) ' MU#' num2str(mu) ' PNR过小，不保留']);
            continue;
        end

        tmp = MUPulses{mu};
        % 截取15s
        tmp = tmp(tmp >= edges(1) & tmp <= edges(2));
        % 对齐零时刻
        tmp = tmp - edges(1);
        % 转换为超声采样率的时刻
        tmp = round(tmp/fsamp*fsampu);

        if length(tmp) <= 15
            disp(['ni#' num2str(ni) ' MU#' num2str(mu) ' 放电时刻过少，不保留']);
            continue;
        end
        if length(tmp) >= 15*35
            disp(['ni#' num2str(ni) ' MU#' num2str(mu) ' 放电时刻过多，不保留']);
            continue;
        end

        pulsesRef{end+1} = tmp;
    end
end
save(['./Data/experiment/25-07-04/M' num2str(motion) 'L1T' num2str(trial) '_pulsesRef.mat'], 'pulsesRef');

