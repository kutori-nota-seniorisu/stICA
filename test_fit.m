clear; clc; close all;
addpath('./Func');

% 肌电采样率2048Hz
fsamp = 2048;
% 超声采样率
fsampu = 1000;

%% 脉冲串时移对齐
% load('./Data/experiment/25-07-04/M1L1T1_decompsRef.mat');
load('./Data/experiment/ICdata/R17/R17_decompsRef.mat')
pulsesAll = decompsRef.Pulses;
muapsAll = decompsRef.MUAPs;
numMU = length(pulsesAll);
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
    pulsesRef{mu}(pulsesRef{mu}<=0) = [];
    pulsesRef{mu}(pulsesRef{mu}>=30*fsampu) = [];
end
% 
% for mu = 1:numMU
%     pulsesRef{mu}(pulsesRef{mu}<=2000) = [];
%     pulsesRef{mu} = pulsesRef{mu} - 2000;
% end

%% STA
load('./Data/experiment/ICdata/R17/USCBSS_DecompResult_MPD50.mat');
twitches = decompoMURaw.Twitch;
% STA窗口100ms，并转换成样本点
winSize = 100/1000*fsampu;
% 数据长度
L = size(twitches, 1);

twitchSTA = cell(0);
varSTA = cell(0);
for i = 1:numMU
    tmpPulse = pulsesRef{i};
    for n = -winSize/2+1:winSize/2
        tmpInd = tmpPulse + n;
        tmpInd(tmpInd <= 0) = [];
        tmpInd(tmpInd >= L) = [];
        tmpTwitch = twitches(tmpInd, :);
        % 存储STA后的曲线
        twitchSTA{i}(n+winSize/2, :) = mean(tmpTwitch);
        % 存储STA后的方差曲线
        varSTA{i}(n+winSize/2, :) = var(tmpTwitch);
    end
    twitchPP{i} = max(twitchSTA{i}) - min(twitchSTA{i});
    figure;
    plot(twitchPP{i}, 'ro');
end
%%
for i = 1:numMU
    tmpTTT = twitchSTA{i};
    figure;
    tiledlayout('vertical', 'TileSpacing', 'tight', 'Padding', 'tight');
    for j = 1:size(tmpTTT, 2)
        nexttile;
        plot(tmpTTT(:, j));
    end
end