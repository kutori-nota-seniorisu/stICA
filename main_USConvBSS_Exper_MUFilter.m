% 对US分解得到的MU进行筛选
clear; clc; close all;
addpath('./Func');
% 导入数据
load(['./Data/experiment/24-06-21/UUS-iEMG/S1M1L1T1_USCBSS_compo25_2.mat']);

%% step1 以MAD和能量占比筛选
saveMUs = [];
saveRows = [];
saveCols = [];
savePulses = {};
saveSources = [];
saveTwitches = [];
saveEnergyRatio = [];
saveCoV = [];
saveMAD = [];

fsampu = 2000;
% 计算放电串的MAD，大于25ms则去除
% 计算估计源在6-14Hz内的能量占比，小于20%则去除
[rn, cn] = size(DecompoResults.sources);
for r = 1:rn
    for c = 1:cn
        tmpPulses = DecompoResults.decompo_pulses{r, c};
        tmpSources = DecompoResults.sources{r, c};
        tmpTwitches = DecompoResults.twitches{r, c};
        tmpCoV = DecompoResults.CoV{r, c};

        for mu = 1:length(tmpPulses)
            % 计算MAD，单位为ms
            MAD = mad(diff(tmpPulses{mu}/fsampu*1000));
            % 计算能量占比
            [psd, freq] = pwelch(tmpSources(:, mu), [], [], [], fsampu);
            % frequency band of interest
            freqbandOI = [6, 14];
            idxBand = freq >= freqbandOI(1) & freq <= freqbandOI(2);
            idxTotal = freq > 0;
            energyInBand = trapz(freq(idxBand), psd(idxBand));
            energyTotal = trapz(freq(idxTotal), psd(idxTotal));
            energyRatio = energyInBand / energyTotal * 100;

            if energyRatio >= 20 && MAD <= 100
                disp(['r=' num2str(r) ',c=' num2str(c) ',mu=' num2str(mu) '保留！']);
                saveMUs(end+1) = mu;
                saveRows(end+1) = r;
                saveCols(end+1) = c;
                savePulses{end+1} = tmpPulses{mu};
                saveSources(:, end+1) = tmpSources(:, mu);
                saveTwitches(:, end+1) = tmpTwitches(:, mu);
                saveEnergyRatio(end+1) = energyRatio;
                saveMAD(end+1) = MAD;
                saveCoV(end+1) = tmpCoV(mu);
            end
        end
    end
end
decompoMURaw.MU = saveMUs;
decompoMURaw.Row = saveRows;
decompoMURaw.Col = saveCols;
decompoMURaw.Pulse = savePulses;
decompoMURaw.Source = saveSources;
decompoMURaw.Twitch = saveTwitches;
decompoMURaw.ER = saveEnergyRatio;
decompoMURaw.MAD = saveMAD;
decompoMURaw.CoV = saveCoV;

clear saveMUs saveRows saveCols savePulses saveSources saveTwitches saveEnergyRatio saveCoV;

%% step2.1 以互相关系数为标准去重
% 第一步：提取所有source信号
sources = decompoMURaw.Source;
numSources = length(decompoMURaw.MU);

% 第二步：计算带有时移的互相关系数矩阵
corrThreshold = 0.3;  % 相关系数阈值

% 初始化相关系数矩阵
crossCorrMatrix = zeros(numSources, numSources);

% sources每一列的均值都近似于零，因而corr的值与xcorr在lag=0处的值近似相等
% 使用corr计算，会忽视具有时延的相似性，因而导致冗余保留。（是否存在时延相似？）
% 那就筛两遍。第一次用corr初步筛选，第二次用xcorr精细筛选。
crossCorrMatrix = abs(corr(sources));

% 第三步：筛选处理
selectedIndices = true(1, numSources);  % 标记哪些元素被保留
processedPairs = false(numSources, numSources);  % 标记哪些对已经处理过

for i = 1:numSources
    if ~selectedIndices(i)  % 如果已经被标记删除，跳过
        continue;
    end

    for j = i+1:numSources
        if ~selectedIndices(j) || processedPairs(i,j)  % 跳过已删除或已处理的配对
            continue;
        end

        % 如果互相关系数大于阈值
        if crossCorrMatrix(i, j) > corrThreshold
            % 比较CoV，选择较小的一个
            if decompoMURaw.CoV(i) <= decompoMURaw.CoV(j)
                selectedIndices(j) = false;  % 删除第j个
            else
                selectedIndices(i) = false;  % 删除第i个
                break;  % 如果第i个被删除，跳出内层循环
            end
        end

        processedPairs(i,j) = true;
        processedPairs(j,i) = true;
    end
end

% 第四步：创建筛选后的结构体数组
decompoMUFiltered.MU = decompoMURaw.MU(selectedIndices);
decompoMUFiltered.Row = decompoMURaw.Row(selectedIndices);
decompoMUFiltered.Col = decompoMURaw.Col(selectedIndices);
decompoMUFiltered.Pulse = decompoMURaw.Pulse(selectedIndices);
decompoMUFiltered.Source = decompoMURaw.Source(:, selectedIndices);
decompoMUFiltered.Twitch = decompoMURaw.Twitch(:, selectedIndices);
decompoMUFiltered.ER = decompoMURaw.ER(selectedIndices);
decompoMUFiltered.MAD = decompoMURaw.MAD(selectedIndices);
decompoMUFiltered.CoV = decompoMURaw.CoV(selectedIndices);

% 输出结果信息
fprintf('原始元素数量: %d\n', numSources);
fprintf('保留元素数量: %d\n', length(decompoMUFiltered.MU));
fprintf('删除元素数量: %d\n', numSources - length(decompoMUFiltered.MU));

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

%% step2.2 以互相关系数为标准去重
% 第一步：提取所有source信号
sources = decompoMUFiltered.Source;
numSources = length(decompoMUFiltered.MU);

% 第二步：计算带有时移的互相关系数矩阵
corrThreshold = 0.3;  % 相关系数阈值

% 初始化相关系数矩阵
crossCorrMatrix = zeros(numSources, numSources);

% sources每一列的均值都近似于零，因而corr的值与xcorr在lag=0处的值近似相等
% 使用corr计算，会忽视具有时延的相似性，因而导致冗余保留。（是否存在时延相似？）
% 那就筛两遍。第一次用corr初步筛选，第二次用xcorr精细筛选。
% crossCorrMatrix = abs(corr(sources));

% 计算每对信号的最大互相关系数（考虑时移）
for i = 1:numSources
    for j = i+1:numSources  % 只计算上三角部分，避免重复
        % 计算互相关系数，考虑时移
        [corrVals, ~] = xcorr(sources(:, i), sources(:, j), 'coeff');
        % 取绝对值最大值（考虑正负相关）
        maxCorr = max(abs(corrVals));
        crossCorrMatrix(i,j) = maxCorr;
        crossCorrMatrix(j,i) = maxCorr;  % 对称矩阵
    end
end

% 第三步：筛选处理
selectedIndices = true(1, numSources);  % 标记哪些元素被保留
processedPairs = false(numSources, numSources);  % 标记哪些对已经处理过

for i = 1:numSources
    if ~selectedIndices(i)  % 如果已经被标记删除，跳过
        continue;
    end

    for j = i+1:numSources
        if ~selectedIndices(j) || processedPairs(i,j)  % 跳过已删除或已处理的配对
            continue;
        end

        % 如果互相关系数大于阈值
        if crossCorrMatrix(i, j) > corrThreshold
            % 比较CoV，选择较小的一个
            if decompoMUFiltered.CoV(i) <= decompoMUFiltered.CoV(j)
                selectedIndices(j) = false;  % 删除第j个
            else
                selectedIndices(i) = false;  % 删除第i个
                break;  % 如果第i个被删除，跳出内层循环
            end
        end

        processedPairs(i,j) = true;
        processedPairs(j,i) = true;
    end
end

% 第四步：创建筛选后的结构体数组
decompoMUFiltered.MU = decompoMUFiltered.MU(selectedIndices);
decompoMUFiltered.Row = decompoMUFiltered.Row(selectedIndices);
decompoMUFiltered.Col = decompoMUFiltered.Col(selectedIndices);
decompoMUFiltered.Pulse = decompoMUFiltered.Pulse(selectedIndices);
decompoMUFiltered.Source = decompoMUFiltered.Source(:, selectedIndices);
decompoMUFiltered.Twitch = decompoMUFiltered.Twitch(:, selectedIndices);
decompoMUFiltered.ER = decompoMUFiltered.ER(selectedIndices);
decompoMUFiltered.MAD = decompoMUFiltered.MAD(selectedIndices);
decompoMUFiltered.CoV = decompoMUFiltered.CoV(selectedIndices);

% 输出结果信息
fprintf('原始元素数量: %d\n', numSources);
fprintf('保留元素数量: %d\n', length(decompoMUFiltered.MU));
fprintf('删除元素数量: %d\n', numSources - length(decompoMUFiltered.MU));

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

%% 更加紧凑的子图绘制
figNum = 30;
numMU = length(decompoMUFiltered.MU);
for axes = 1:ceil(numMU/figNum)
    figure;
    t = tiledlayout(6, 5, "TileSpacing", 'tight', 'Padding', 'compact');
    for i = (axes-1)*figNum+1:min(axes*figNum, numMU)
        source = decompoMUFiltered.Source(:, i);
        pulse = decompoMUFiltered.Pulse{i};
        ER = decompoMUFiltered.ER(i);
        nexttile;
        plot(source);
        hold on;
        plot(pulse, source(pulse), 'ro');
        title(['MU # ' num2str(i) ',ER=' num2str(ER)]);
        xticks(0:1e4:3e4);
        xticklabels(0:5:15);
    end
    xlabel(t, 'time (s)');
    ylabel(t, 'amplitude');
end

numMU = length(decompoMUFiltered.MU);
for axes = 1:ceil(numMU/figNum)
    figure;
    t = tiledlayout(6, 5, "TileSpacing", 'tight', 'Padding', 'compact');
    for i = (axes-1)*figNum+1:min(axes*figNum, numMU)
        source = decompoMUFiltered.Source(:, i);
        [pxx, f] = pwelch(source, [], [], [], fsampu);
        nexttile;
        % semilogy(f, pxx);
        plot(f, pxx);
        title(['MU # ' num2str(i)]);
    end
    xlabel(t, 'freq (s)');
    ylabel(t, '功率谱密度');
end

%% 脉冲串匹配
load(['./Data/experiment/ICdata/R' sub '/pulsesRef.mat']);
% 匹配容差为0~5ms
winSize = [0, 5]/1000*fsampu;
% lim=100ms，转换成样本点作为输入参数
lim = 100/1000*fsampu;

fsampu = 1000;
dIPI = round(0.010*fsampu);

matchResultRaw = [];
for i = 1:length(pulsesRef)
    for j = 1:length(decompoMUFiltered.MU)
        [PulseStat,SourceID,Lag,Sens,Miss,FalseAlarms,Specificity] = testSinResults(pulsesRef{i},decompoMUFiltered.Pulse{j},dIPI,0);
        [Sen,FA,Pre,Spe,Acc] = accEvaluation(decompoMUFiltered.Pulse{j},pulsesRef{i},dIPI,100);
        [rr, ~] = RoA(decompoMUFiltered.Pulse{j},pulsesRef{i},100, dIPI);
        matchResultRaw(end+1,:) = [i,j,rr,Lag,Sens,Sen,Miss,FalseAlarms,FA,Specificity,Spe,Pre,Acc];
    end
end
matchResultRaw = array2table(matchResultRaw,'VariableNames',{'ref','decomp','RoA','Lag','Sens1','Sens2', 'Miss', 'FA1', 'FA2', 'Spe1','Spe2', 'Pre', 'Acc'});

plotDecomps(decompoMUFiltered.Pulse, [], fsampu, 0, 0, []);
% plotDecomps(decompoMURaw.Pulse, [], fsampu, 0, 0, []);
plotDecomps(pulsesRef, [], fsampu, 0, 0, []);

% matchResultRaw = [];
% for i = 1:length(pulsesRef)
%     if length(pulsesRef{i}) < 150 || length(pulsesRef{i}) > 450
%         disp(['i=' num2str(i) '参考脉冲序列无效']);
%         continue;
%     end
%     for j = 1:length(decompoMUFiltered.MU)
%         xcor = fxcorr(pulsesRef{i}, decompoMUFiltered.Pulse{j}, lim);
%         winSums = zeros(1, 2*lim+1-(winSize(2)-winSize(1)));
%         for k = 1:2*lim+1-(winSize(2)-winSize(1))
%             winSums(k) = sum(xcor(k:k+(winSize(2)-winSize(1))));
%         end
%         [maxMatchSum, maxIdx] = max(winSums);
%         % 最佳匹配的时移
%         Lag = maxIdx - winSize(1) - lim - 1;
%         % 匹配上的脉冲个数，真阳性
%         TP = maxMatchSum;
%         % 假阳性
%         FP = length(decompoMUFiltered.Pulse{j}) - maxMatchSum;
%         % 假阴性
%         FN = length(pulsesRef{i}) - maxMatchSum;
%         % 匹配率
%         matchRatio = TP / (TP + FN) * 100;
%         % RoA
%         rr = TP / (TP + FP + FN);
%         % rr = RoA(decompoMUFiltered(j).pulse, pulsesRef{i}, 100, 5);
% 
%         matchResultRaw(end+1, :) = [i, j, Lag, TP, FP, FN, matchRatio, rr];
%     end
% end
% 
% matchResultRaw = array2table(matchResultRaw, 'VariableNames', {'ref', 'decomp', 'Lag', 'TP', 'FP', 'FN', 'match ratio', 'RoA'});

%% 参考脉冲
data = importdata('./Data/EMG/24-06-21/iEMG_S1_M1_level1_trial1_24-06-21_UUS.eaf'); % 读取eaf文件
fsampu = 1000; % 采样率
muNum = max(data.data(:,2)); % MU的个数
pulsesRef = {};
for mu = 1:muNum
    % iPulses就是这个eaf文件里分解得到的spike train，每个cell表示一个MU，里面的数字是该MU每次放电的时刻
    tmp = round(data.data(find(data.data(:,2)==mu),1)'*fsampu);
    tmp = tmp(tmp >= 3*fsampu & tmp <= 13*fsampu);
    tmp = tmp - 3*fsampu;
    pulsesRef{mu} = tmp;
    CoV(mu) = std(diff(pulsesRef{mu}))/mean(diff(pulsesRef{mu}));
    MAD(mu) = mad(diff(pulsesRef{mu}/fsampu*1000));
end
plotDecomps(pulsesRef, [], fsampu, 0, 0, []);

%% 计算筛选后MU的PNR
for sub = 16%[3,4,5,7,10,11,12,14,15,16,17,18]
    % 导入数据
    load(['./Data/experiment/ICdata/R' num2str(sub) '/USCBSS_DecompResult_MPD60.mat']);
    PNRs = [];
    Rows = [];
    for i = 1:length(decompoMUFiltered.MU)
        IPT = decompoMUFiltered.Source(:, i);
        pulse = decompoMUFiltered.Pulse{i};
        Rows(i) = decompoMUFiltered.Row(i);
        PNRs(i) = 10*log10( mean( IPT(pulse).^2 ) / mean( IPT(setdiff(1:length(IPT), pulse)).^2 ) );
    end
    % save(['./Data/experiment/ICdata/R' num2str(sub) '/PNRs.mat'], 'PNRs', 'Rows');
end

PNRsAll = [];
RowsAll = [];
for sub = [3,4,5,7,10,11,12,14,15,16,17,18]
    % 导入数据
    load(['./Data/experiment/ICdata/R' num2str(sub) '/PNRs.mat']);
    PNRsAll = [PNRsAll, PNRs];
    RowsAll = [RowsAll, Rows];
end
mean(PNRsAll)
std(PNRsAll)
median(PNRsAll)

% plotDecomps({pulsesRef{5},decompoMUFiltered(2).pulse}, [], 1000, 0, 0, []);

%% 绘制23*25个区域的估计源信号
fsampu = 1000;
[rn, cn] = size(DecompoResults.sources);
for r = 1:rn
    for c = 1:cn
        tmpPulses = DecompoResults.decompo_pulses{r, c};
        tmpSources = DecompoResults.sources{r, c};
        tmpTwitches = DecompoResults.twitches{r, c};
        tmpCoV = DecompoResults.CoV{r, c};

        ax = figure;
        for mu = 1:length(tmpPulses)
            % 计算MAD，单位为ms
            MAD = mad(diff(tmpPulses{mu}/fsampu*1000));
            % 计算能量占比
            [psd, freq] = pwelch(tmpSources(:, mu), [], [], [], fsampu);
            % frequency band of interest
            freqbandOI = [6, 14];
            idxBand = freq >= freqbandOI(1) & freq <= freqbandOI(2);
            idxTotal = freq > 0;
            energyInBand = trapz(freq(idxBand), psd(idxBand));
            energyTotal = trapz(freq(idxTotal), psd(idxTotal));
            energyRatio = energyInBand / energyTotal * 100;

            subplot(10,5,mu+floor((mu-1)/5)*5);
            plot(tmpTwitches(:,mu));
            xlim([4000,8000]);
            xticks(4000:1000:8000);
            xticklabels(4:1:8);
            xlabel('t/s'); ylabel('amplitude');
            title(['twitch mu=' num2str(mu)])

            subplot(10,5,mu+ceil(mu/5)*5);
            plot(tmpSources(:,mu));
            xlim([4000,8000]);
            xticks(4000:1000:8000);
            xticklabels(4:1:8);
            xlabel('t/s'); ylabel('amplitude');
            title(['source mu=' num2str(mu) '，MAD=' num2str(MAD) ',ER=' num2str(energyRatio) '%']);
        end

        set(gcf,'unit','normalized','position',[0,0,1,1]);
        saveas(ax, ['./Results/L1T1P1/r' num2str(r) 'c' num2str(c)], 'png');
        % saveas(ax, ['./Results/R' sub '_LowPass/r' num2str(r) 'c' num2str(c)], 'png');
        close;
    end
end

%% 2s结果
fsampu = 1000;
[rn, cn] = size(DecompoResults.sources);
for r = 1:rn
    for c = 1:cn
        tmpPulses = DecompoResults.decompo_pulses{r, c};
        tmpSources = DecompoResults.sources{r, c};
        tmpTwitches = DecompoResults.twitches{r, c};
        tmpCoV = DecompoResults.CoV{r, c};

        ax = figure;
        for mu = 1:length(tmpPulses)
            % 计算MAD，单位为ms
            MAD = mad(diff(tmpPulses{mu}/fsampu*1000));
            % 计算能量占比
            [psd, freq] = pwelch(tmpSources(:, mu), [], [], [], fsampu);
            % frequency band of interest
            freqbandOI = [6, 14];
            idxBand = freq >= freqbandOI(1) & freq <= freqbandOI(2);
            idxTotal = freq > 0;
            energyInBand = trapz(freq(idxBand), psd(idxBand));
            energyTotal = trapz(freq(idxTotal), psd(idxTotal));
            energyRatio = energyInBand / energyTotal * 100;

            subplot(10,5,mu+floor((mu-1)/5)*5);
            plot(tmpTwitches(:,mu));
            % xlim([4000,8000]);
            % xticks(4000:1000:8000);
            xticklabels(0:1:2);
            xlabel('t/s'); ylabel('amplitude');
            title(['twitch mu=' num2str(mu)])

            subplot(10,5,mu+ceil(mu/5)*5);
            plot(tmpSources(:,mu));
            % xlim([4000,8000]);
            % xticks(4000:1000:8000);
            xticklabels(0:1:2);
            xlabel('t/s'); ylabel('amplitude');
            title(['source mu=' num2str(mu) '，MAD=' num2str(MAD) ',ER=' num2str(energyRatio) '%']);
        end

        set(gcf,'unit','normalized','position',[0,0,1,1]);
        saveas(ax, ['./Results/L1T1P1_3~5s/r' num2str(r) 'c' num2str(c)], 'png');
        % saveas(ax, ['./Results/R' sub '_LowPass/r' num2str(r) 'c' num2str(c)], 'png');
        close;
    end
end

%% 统计三元组出现情况
% 假设 ccMax 是互相关系数方阵，对角线为1
% 设置相关系数阈值
threshold = 0.3;

% 获取矩阵大小
n = size(ccMax, 1);

% 存储所有符合条件的三元组
triplets = [];

% 遍历所有可能的三元组组合
for a = 1:n
    for b = (a+1):n
        % 检查a和b的相关系数是否小于阈值
        if abs(ccMax(a,b)) < threshold
            for c = 1:n
                % c不能等于a或b
                if c == a || c == b
                    continue;
                end

                % 检查a和c，b和c的相关系数是否都大于阈值
                if abs(ccMax(a,c)) > threshold && abs(ccMax(b,c)) > threshold
                    % 找到符合条件的三元组
                    triplets = [triplets; a, b, c];

                    % 输出信息
                    fprintf('发现三元组: a=%d, b=%d, c=%d\n', a, b, c);
                    fprintf('  ccMax(a,b)=%.3f < %.3f\n', ccMax(a,b), threshold);
                    fprintf('  ccMax(a,c)=%.3f > %.3f\n', ccMax(a,c), threshold);
                    fprintf('  ccMax(b,c)=%.3f > %.3f\n', ccMax(b,c), threshold);
                    fprintf('---\n');
                end
            end
        end
    end
end

% 输出统计结果
if isempty(triplets)
    fprintf('未找到符合条件的三元组\n');
else
    fprintf('\n=== 统计结果 ===\n');
    fprintf('共找到 %d 个符合条件的三元组\n', size(triplets, 1));

    % 统计每个节点出现的频率
    allNodes = unique([triplets(:,1); triplets(:,2); triplets(:,3)]);
    nodeCount = zeros(length(allNodes), 2);

    for i = 1:length(allNodes)
        node = allNodes(i);
        count = sum(triplets(:) == node);
        nodeCount(i, :) = [node, count];
    end

    % 按出现频率排序
    [~, idx] = sort(nodeCount(:,2), 'descend');
    nodeCount = nodeCount(idx, :);

    fprintf('\n节点出现频率:\n');
    for i = 1:size(nodeCount, 1)
        fprintf('  节点 %d: 出现 %d 次\n', nodeCount(i,1), nodeCount(i,2));
    end
end

%% 二次筛选
% 计算带有时延的互相关系数
[cc, ~] = xcorr(decompoMUFiltered.Source, 'coeff');
numMU = length(decompoMUFiltered.MU);
% 转成三维数组
cc = reshape(cc', numMU, numMU, []);
% 计算互相关系数的最大值
ccMax = max(abs(cc), [], 3);
% % 生成掩膜。只关注矩阵的上三角部分，不包含对角线，避免重复比较。
% % mask是一个逻辑矩阵，下三角及主对角线部分均为零，上三角部分只有大于阈值的位置才为一。
% mask = triu(ccMax > 0.7, 1);
% % 初始化保留索引
% toKeep = true(size(decompoSourceAll, 2), 1);
% % 遍历mask矩阵，找到所有要删除的j
% for i = 1:size(decompoSourceAll, 2)
%     j_toRemove = find(mask(i, :));
%     toKeep(j_toRemove) = false;
% end
%
% decompoSourceAll = decompoSourceAll(:, toKeep);
% decompoCoVAll = decompoCoVAll(toKeep);
% decompoPulseAll = decompoPulseAll(toKeep);

%%
figure;
subplot(2,1,1);
plot(decompoMUFiltered(1).sourceRaw);
xlim([0,4000]);
xticks(0:1000:4000);
xticklabels(0:1:4);
xlabel('t/s'); ylabel('amplitude');
title('一阶段迭代结果')
subplot(2,1,2);
plot(decompoMUFiltered(1).source);
hold on;
plot(decompoMUFiltered(1).pulse, decompoMUFiltered(1).source(decompoMUFiltered(1).pulse), 'ro');
xlim([0,4000]);
xticks(0:1000:4000);
xticklabels(0:1:4);
xlabel('t/s'); ylabel('amplitude');
title('二阶段迭代结果')


%%
numMU = length(decompoMURaw.MU);
for i = 1:numMU
    sss = decompoMURaw.Source(:, i);
    ttt = decompoMURaw.Twitch(:, i);
    ppp = decompoMURaw.Pulse{i};
    figure;
    subplot(2,1,1);
    plot(ttt);
    xticklabels(0:1:10);
    xlabel('t (s)')
    title('twitch')
    subplot(2,1,2);
    plot(sss);
    hold on;
    plot(ppp, sss(ppp), 'ro');
    xticklabels(0:1:10);
    xlabel('t (s)')
    title('estimated source')
    set(gcf,'unit','normalized','position',[0.3,0.4,0.4,0.3]);
end

%%
numMU = length(decompoMUFiltered.MU);
for i = 1:numMU
    sss = decompoMUFiltered.Source(:, i);
    ttt = decompoMUFiltered.Twitch(:, i);
    ppp = decompoMUFiltered.Pulse{i};
    figure;
    subplot(2,1,1);
    plot(ttt);
    xticklabels(0:1:10);
    xlabel('t (s)')
    title('twitch')
    subplot(2,1,2);
    plot(sss);
    hold on;
    plot(ppp, sss(ppp), 'ro');
    xticklabels(0:1:10);
    xlabel('t (s)')
    title('estimated source')
    set(gcf,'unit','normalized','position',[0.3,0.4,0.4,0.3]);
end