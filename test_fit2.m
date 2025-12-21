clear; clc; close all;
addpath('./Func');

% EMG采样率
fsamp = 2048;
% US采样率
fsampu = 1000;
% STA窗口100ms
winSize = 100/1000*fsampu;

STAFile = './Data/experiment/ICdata/R3/output_info.mat';
BSSFile = './Data/experiment/ICdata/R3/USCBSS_DecompResult.mat';
emgFile = './Data/experiment/ICdata/R3/R3.mat';
% tviFile = './Data/experiment/ICdata/R17/v_2d_all.mat';


load(emgFile);
% load(tviFile);
load(STAFile);
load(BSSFile);

newdata = Data{1, 2};
newdata(newdata > 32768) = newdata(newdata > 32768) - 2^16;
trigger = newdata(end, :);
[~, edges] = maxk(trigger, 2);
edges = sort(edges);


%% twitch curve from BSS
% 用EMG-MU的spike对US-MU的twitch train进行STA，得到twitch curve
twitches = decompoMUFiltered.Twitch(:,1:4);
% twitches = decompoMURaw.Twitch;
% 用EMG对twitch train进行STA的结果存储
twitchCurveBSS = cell(0);
varBSS = cell(0);
twitchPP = cell(0);
numMU = length(output_info);
for mu = 1:numMU
    output = output_info(mu);
    if isempty(output.Number_active_MUs)
        continue;
    end  
    tmpPulse = output.MUFiring_in_region;
    tmpPulse = tmpPulse-edges(1);
    tmpPulse(tmpPulse<=50) = [];
    tmpPulse(tmpPulse>=(30000-50)) = [];
    for n = -winSize/2+1:winSize/2
        tmpInd = tmpPulse + n;
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
% [M, N] = size(activation{1});
M = 39; N = 128;
% 窗口大小
Row = 10; Col = 10;
% 窗口移动距离
dRow = 5; dCol = 5;
for mu = 1:4
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
twitchCurves = cell(0);
for mu = 1:numMU
    output = output_info(mu);
    if isempty(output.Number_active_MUs)
        continue;
    end
    tmpImg = cell(0);

    tmpInd = output.positive_indicies_lin;
    tmptmp = zeros([39, 128]);
    tmptmp(tmpInd) = 2;
    tmpImg{1} = tmptmp;

    tmpInd = output.negative_indicies_lin;
    tmptmp = zeros([39, 128]);
    tmptmp(tmpInd) = 2;
    tmpImg{2} = tmptmp;

    twitchAreaSTA{mu} = tmpImg;

    tmpCurve = cell(0);
    tmpCurve{1} = output.positive_STA;
    tmpCurve{2} = output.negative_STA;
    twitchCurves{mu} = tmpCurve;
end

%% 将US-MU与EMG-MU进行候选匹配
% 存储匹配结果
matchResult = [];
for mu = 1:numMU
    figure;
    plot(twitchPP{mu}, 'ro');
    [ppMax, idxMax] = max(twitchPP{mu});
    [ppMin, idxMin] = min(twitchPP{mu});
    if isempty(ppMax)
        matchResult(end+1, :) = [mu,mu,mu,mu,mu];
    else
        matchResult(end+1, :) = [mu, idxMax, ppMax, idxMin, ppMin];
    end
end

%% twitch curve相关系数计算
rho = zeros(0);
for i = 1:size(matchResult, 1)
    output = output_info(i);
    if isempty(output.Number_active_MUs)
        continue;
    end
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
for mu = 1:numMU
    tmpTTT = twitchCurveBSS{mu};
    if isempty(tmpTTT)
        continue;
    end
    figure;
    tiledlayout('flow', 'TileSpacing', 'tight', 'Padding', 'tight');
    for j = 1:size(tmpTTT, 2)
        nexttile;
        plot(tmpTTT(:, j));
        title(['MU #' num2str(j)])
    end
end

plotDecomps(decompoMUFiltered.Pulse, [], fsampu, 0, 0, []);
