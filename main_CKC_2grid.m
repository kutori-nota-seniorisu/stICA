%% 用CKC反解EMG
clear; clc; close all;
addpath('./Func');

% EMG采样频率
fsamp = 2048;
% 是否导入反解结果，1为导入已反解的结果，0为用CKC反解
isLoad = 0;

% for sub = [3,4,5,7,10,11,12,14,15,16,17,18]
% sub = 17;
% emgFile = ['./Data/experiment/ICdata/R' num2str(sub) '/R' num2str(sub) '.mat'];

motion = 2;
trial = 1;
disp(['M' num2str(motion) 'L1T' num2str(trial)]);
emgFile = ['./Data/EMG/25-07-04/M' num2str(motion) 'L1T' num2str(trial) '.mat'];

try
    load(emgFile);
catch ME
    rethrow(ME);
end

%% step1 EMG处理，根据不同数据类型进行修改
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

%% step2 CKC反解EMG
sEMGArray = [];
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

    if ~isLoad
        decomps{ni} = IPTExtraction_gCKC4GUI_v3(decompData,decoderParameters); % classic gradient CKC
        decomps{ni}.decompChannelInd = decompChannelInd;
    end
    [rn, cn] = size(dataarray);

    % 由元胞转成数组，方便后续计算MUAP
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

if isLoad
    load('./Data/experiment/25-07-04/M1L1T1_decompsRaw.mat');
end

%% step3 绘制IPT，查看分解结果
for ni = 1:2
    plotDecomps(decomps{ni}.MUPulses, [], fsamp, 0, 0, []);
    IPTs = decomps{ni}.IPTs;
    PNRs = decomps{ni}.PNRs;
    MUPulses = decomps{ni}.MUPulses;

    figure;
    t = tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'compact');
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

%% step4 按照筛选要求，保存参考脉冲串
pulsesRef = {};
CoV = [];
for ni = 1:2
    PNRs = decomps{ni}.PNRs;
    MUPulses = decomps{ni}.MUPulses;
    for mu = 1:length(MUPulses)
        if PNRs(mu) <= 20
            disp(['ni#' num2str(ni) ' MU#' num2str(mu) ' PNR过小，不保留']);
            continue;
        end

        tmpPulse = MUPulses{mu};
        % 截取15s
        tmpPulse = tmpPulse(tmpPulse >= edges(1) & tmpPulse <= edges(2));
        % 对齐零时刻
        tmpPulse = tmpPulse - edges(1);
        % 放电变异率
        tmpCoV = std(diff(tmpPulse))/mean(diff(tmpPulse));

        if tmpCoV >= 0.5
            disp(['ni#' num2str(ni) ' MU#' num2str(mu) ' 放电变异率过大，不保留']);
            continue;
        end
        if length(tmpPulse) <= 15
            disp(['ni#' num2str(ni) ' MU#' num2str(mu) ' 放电时刻过少，不保留']);
            continue;
        end
        if length(tmpPulse) >= 15*35
            disp(['ni#' num2str(ni) ' MU#' num2str(mu) ' 放电时刻过多，不保留']);
            continue;
        end

        pulsesRef{end+1} = tmpPulse;
        CoV(end+1) = tmpCoV;
    end
end

%% step5 参考脉冲串去重
numMU = length(pulsesRef);
dIPI = round(0.0005*fsamp);
spikeRoAMatrix = zeros(numMU, numMU);
% mmm=[];
for i = 1:numMU
    for j = i+1:numMU
        % 计算两个脉冲串之间的RoA
        [rr, lag] = RoA(pulsesRef{j}, pulsesRef{i}, 100, dIPI);
        spikeRoAMatrix(i, j) = rr;
        % mmm(end+1,:) = [i,j,rr,lag];
    end
end
% mmm = array2table(mmm,'VariableNames',{'ref','decomp','RoA','Lag'});

% 标记哪些元素被保留
selectedIndices = true(1, numMU);
% 标记哪些对已经处理过
processedPairs = false(numMU, numMU);

for i = 1:numMU
    % 如果已经被标记删除，跳过
    if ~selectedIndices(i)
        continue;
    end

    for j = i+1:numMU
        % 跳过已删除或已处理的配对
        if ~selectedIndices(j) || processedPairs(i, j)
            continue;
        end

        % 如果RoA大于0.9
        if spikeRoAMatrix(i, j) >= 0.9
            % 比较CoV，选择较小的一个
            if CoV(i) <= CoV(j)
                % 删除第j个
                selectedIndices(j) = false;
            else
                % 删除第i个
                selectedIndices(i) = false;
                % 如果第i个被删除，跳出内层循环
                break;
            end
        end

        processedPairs(i, j) = true;
    end
end

pulsesRef = pulsesRef(selectedIndices);
CoV = CoV(selectedIndices);

% 输出结果信息
fprintf('原始元素数量: %d\n', numMU);
fprintf('保留元素数量: %d\n', length(pulsesRef));
fprintf('删除元素数量: %d\n', numMU - length(pulsesRef));

% 显示被保留的索引
fprintf('被保留的索引: ');
fprintf('%d ', find(selectedIndices));
fprintf('\n');

% 显示被删除的索引
deletedIndices = find(~selectedIndices);
fprintf('被删除的索引: ');
if ~isempty(deletedIndices)
    fprintf('%d ', deletedIndices);
else
    fprintf('无');
end
fprintf('\n');

%% step6 绘制MUAP
numMU = length(pulsesRef);
winSize = 128;
lenEMG = diff(edges)+1;

% 存储所有的MUAP阵列
arrayMUAP = cell(1, numMU);

for mu = 1:numMU
    tmp = pulsesRef{mu};
    MUAP = [];
    for n = -winSize/2+1:winSize/2
        tmpInd = tmp + n;
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

%% step7 保存结果
% 分解的原始结果
% save(['./Data/experiment/ICdata/R' num2str(sub) '/R' num2str(sub) '_decompsRaw.mat'], 'decomps');
% save(['./Data/experiment/25-07-04/M' num2str(motion) 'L1T' num2str(trial) '_decompsRaw.mat'], 'decomps');

% 用于充当参考的MU结果
decompsRef.MUAPs = arrayMUAP;
decompsRef.Pulses = pulsesRef;
decompsRef.CoVs = CoV;

% save(['./Data/experiment/ICdata/R' num2str(sub) '/R' num2str(sub) '_decompsRef.mat'], 'decompsRef');
% save(['./Data/experiment/25-07-04/M' num2str(motion) 'L1T' num2str(trial) '_decompsRef.mat'], 'decompsRef');
% end