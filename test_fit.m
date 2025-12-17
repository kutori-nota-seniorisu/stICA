clear; clc; close all;
addpath('./Func');

% EMG采样率
fsamp = 2048;
% US采样率
fsampu = 1000;
% STA结果路径
STAFile = './Data/experiment/ICdata/R17/STA_DecompResult.mat';
% EMG分解结果文件路径
emgFile = './Data/experiment/ICdata/R17/R17_decompsRef.mat';
% US文件路径
% 导入EMG分解结果
load(emgFile);
% 导入STA结果
load(STAFile);

%% 脉冲串时移对齐
% load('./Data/experiment/25-07-04/M1L1T1_decompsRef.mat');
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

% for mu = 1:numMU
%     pulsesRef{mu}(pulsesRef{mu}<=2000) = [];
%     pulsesRef{mu} = pulsesRef{mu}-2000;
% end

%% twitch curve from BSS
% 用EMG-MU的spike对US-MU的twitch train进行STA，得到twitch curve
% load('./Data/experiment/ICdata/R17/USCBSS_DecompResult.mat')
load('./Data/experiment/ICdata/R17/USCBSS_DecompResult_C15.mat');
% load('./Results/251208/M1L1T1_USCBSS_DecompResultR20C15.mat');
% load('./Data/experiment/ICdata/R17/USCBSS_DecompResult_MPD50.mat');
twitches = decompoMUFiltered.Twitch;
% twitches = decompoMURaw.Twitch;
% STA窗口100帧
winSize = 100;
% 数据长度
L = size(twitches, 1);
% 用EMG对twitch train进行STA的结果存储
twitchCurveBSS = cell(0);
varBSS = cell(0);
twitchPP = cell(0);
for mu = 1:numMU
    tmpPulse = pulsesRef{mu};
    for n = -winSize/2+1:winSize/2
        tmpInd = tmpPulse + n;
        tmpInd(tmpInd <= 0) = [];
        tmpInd(tmpInd >= L) = [];
        tmpTwitch = twitches(tmpInd, :);
        % 存储STA后的曲线
        twitchCurveBSS{mu}(n+winSize/2, :) = mean(tmpTwitch);
        % 存储STA后的方差曲线
        varBSS{mu}(n+winSize/2, :) = var(tmpTwitch);
    end
    directions = sum(twitchCurveBSS{mu}(1:winSize/2,:)./varBSS{mu}(1:winSize/2,:))-sum(twitchCurveBSS{mu}(winSize/2+1:end,:)./varBSS{mu}(winSize/2+1:end,:));
    directionsSign = sign(directions);
    twitchPP{mu} = (max(twitchCurveBSS{mu})-min(twitchCurveBSS{mu})).*-directionsSign;
end

%% twitch area from BSS
% 存储BSS区域
twitchAreaBSS = cell(0);
[M, N] = size(activation{1});
% 窗口大小
Row = 10; Col = 10;
% 窗口移动距离
dRow = 5; dCol = 5;
for mu = 1:length(decompoMUFiltered.MU)
    tmpImg = zeros([M, N]);
    r = decompoMUFiltered.Row(mu);
    c = decompoMUFiltered.Col(mu);
    winRow = (1:Row)+(r-1)*dRow;
    winCol = (1:Col)+(c-1)*dCol;
    winRow(winRow>M) = [];
    winCol(winCol>N) = [];
    tmpImg(winRow, winCol) = 1;
    twitchAreaBSS{mu} = tmpImg;
end

%% twitch area from STA
twitchAreaSTA = cell(0);
for mu = 1:numMU
    tmpImg = cell(0);
    uuimage = activation{mu};
    threshVal = 0.65*max(abs(uuimage(:)));
    tmp = uuimage;
    for na = 1:2
        if na == 1
            tmpInd = find(tmp>threshVal);
        else
            tmpInd = find(tmp<-threshVal);
        end
        tmptmp = zeros([M, N]);
        tmptmp(tmpInd) = 2;
        tmpImg{na} = tmptmp;
    end
    twitchAreaSTA{mu} = tmpImg;
end

%% 将US-MU与EMG-MU进行候选匹配
% 存储匹配结果
matchResult = [];
for mu = 1:numMU
    figure;
    plot(twitchPP{mu}, 'ro');
    [ppMax, idxMax] = max(twitchPP{mu});
    [ppMin, idxMin] = min(twitchPP{mu});
    matchResult(end+1, :) = [mu, idxMax, ppMax, idxMin, ppMin];
end

%% twitch curve相关系数计算
rho = zeros(0);
for i = 1:size(matchResult, 1)
    % 参考MU序号
    refMU = matchResult(i, 1);
    % 正PP值MU序号
    decompMUP = matchResult(i, 2);
    % 负PP值MU序号
    decompMUN = matchResult(i, 4);
    % 参考MU对应的twitch curve
    twitchSTA = twitchCurves{refMU};
    % US-MU对应的twitch curve
    twitchBSS = {twitchCurveBSS{refMU}(:, decompMUP), twitchCurveBSS{refMU}(:, decompMUN)};
    % 参考MU对应的twitch area
    areaSTA = twitchAreaSTA{refMU};
    % US-MU对应的twitch area
    areaBSS = {twitchAreaBSS{decompMUP}, twitchAreaBSS{decompMUN}};

    figure;
    tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    for na = 1:2
        % twitch curve对比图
        nexttile;
        yyaxis left;
        plot(twitchBSS{na});
        yyaxis right;
        plot(twitchSTA{na});
        legend('BSS', 'STA');
        if ~isempty(twitchSTA{na})
            rho(i,na) = corr(twitchSTA{na}', twitchBSS{na});
        else
            rho(i,na)=0;
        end
        title(['rho=' num2str(rho(i,na)) ',na=' num2str(na)]);

        % twitch area对比图
        nexttile;
        imagesc(areaSTA{na}+areaBSS{na});
        title(['na=' num2str(na)]);
        colorbar;
    end
end


%%
for i = 1:size(matchResult, 1)
    % 参考MU序号
    refMU = matchResult(i, 1);
    % 正PP值MU序号
    decompMUP = matchResult(i, 2);
    % 负PP值MU序号
    decompMUN = matchResult(i, 4);


    figure;
    tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    for na = 1:2

        nexttile;
        imagesc(areaSTA{na}+areaBSS{na});

        colorbar;
    end
end
%%
for mu = 1:numMU
    tmpTTT = twitchCurveBSS{mu};
    figure;
    tiledlayout('flow', 'TileSpacing', 'tight', 'Padding', 'tight');
    for j = 1:size(tmpTTT, 2)
        nexttile;
        plot(tmpTTT(:, j));
        title(['MU #' num2str(j)])
    end
end

plotDecomps(decompoMUFiltered.Pulse, [], fsampu, 0, 0, []);

%%
load('./Data/experiment/ICdata/R17/v_2d_all.mat');
TVIData = cat(3, zeros(119, 128, 2), v_2d_all);
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

[M, N, L] = size(TVIDataFilter);
% 窗口大小
Row = 10; Col = 10;
% Row = 8; Col = 8;
% 窗口移动距离
dRow = 5; dCol = 5;
% dRow = 4; dCol = 4;

for mu = 1:numMU
    tmpPulse = decompoMUFiltered.Pulse{mu};
    frameSTA = zeros(0);
    varBSS = zeros(0);
    for n = -winSize/2+1:winSize/2
        tmpInd = tmpPulse + n;
        tmpInd(tmpInd <= 0) = [];
        tmpInd(tmpInd >= L) = [];
        tmpTVI = TVIDataFilter(:, :, tmpInd);
        % 存储STA图像
        frameSTA(:, :, n+winSize/2) = mean(tmpTVI, 3);
        % 存储方差图像
        varBSS(:, :, n+winSize/2) = var(tmpTVI, 0, 3);
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
    set(gcf,'unit','normalized','position',[0.1,0.5,0.35,0.3]);
end