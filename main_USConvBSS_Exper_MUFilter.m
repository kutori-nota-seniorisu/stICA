%% 仿真数据Pulse筛选
clear; clc; close all;
addpath('./Func');
Sub = '16';
% 导入数据
load(['./Data/experiment/ICdata/R' Sub '/USCBSS_compo25.mat']);
%% step1 以MAD和能量占比筛选
decompoMURaw = struct('row', {}, 'col', {}, 'pulse', {}, 'source', {}, 'CoV', {});
fsampu = 1000;
% decompoSources = [];
% 计算放电串的MAD，大于25ms则去除
% 计算估计源在6-14Hz内的能量占比，小于20%则去除
for r = 16:23
    for c = 1:25
        tmpPulses = DecompoResults.decompo_pulses{r, c};
        tmpSources = DecompoResults.sources{r, c};
        tmpSourcesRaw = DecompoResults.sourceFirst{r, c};
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
            plot(tmpSourcesRaw(:,mu));
            xlim([0,4000]);
            xticks(0:1000:4000);
            xticklabels(0:1:4);
            xlabel('t/s'); ylabel('amplitude');
            title(['mu=' num2str(mu) '一阶段迭代'])

            subplot(10,5,mu+ceil(mu/5)*5);
            plot(tmpSources(:,mu));
            xlim([0,4000]);
            xticks(0:1000:4000);
            xticklabels(0:1:4);
            xlabel('t/s'); ylabel('amplitude');
            title(['mu=' num2str(mu) '二阶段迭代，MAD=' num2str(MAD) ',ER=' num2str(energyRatio) '%']);

            % L = 30000;
            % tmp = abs(fft(tmpSources(:, mu))/L);
            % sFFT = tmp(1:L/2+1);
            % sFFT(2:end-1) = 2*sFFT(2:end-1);
            % E = sFFT.^2;
            % freq = (0:1:L/2)*fsampu/L;
            % freqbandOI = [6, 14];
            % idxBand = freq >= freqbandOI(1) & freq <= freqbandOI(2);
            % idxTotal = freq > 0;
            % energyInBand = trapz(freq(idxBand), E(idxBand));
            % energyTotal = trapz(freq(idxTotal), E(idxTotal));
            % energyRatio = energyInBand / energyTotal * 100;

            % if MAD > 25
            %     disp(['r=' num2str(r) ',c=' num2str(c) ',mu=' num2str(mu) 'MAD过大']);
            % end
            % if energyRatio < 20
            %     disp(['r=' num2str(r) ',c=' num2str(c) ',mu=' num2str(mu) '能量占比过小']);
            % end

            if MAD <= 25 && energyRatio >= 20
                disp(['r=' num2str(r) ',c=' num2str(c) ',mu=' num2str(mu) '保留！']);
                decompoMURaw(end+1).MU = mu;
                decompoMURaw(end).row = r;
                decompoMURaw(end).col = c;
                decompoMURaw(end).pulse = tmpPulses{mu};
                decompoMURaw(end).source = tmpSources(:, mu);
                decompoMURaw(end).sourceRaw = tmpSourcesRaw(:, mu);
                decompoMURaw(end).energyRatio = energyRatio;
                decompoMURaw(end).CoV = tmpCoV(mu);
                % decompoSources(:, end+1) = tmpSources(:, mu);
            end
        end

        set(gcf,'unit','normalized','position',[0,0,1,1]);
        saveas(ax, ['./Results/R' Sub '/r' num2str(r) 'c' num2str(c)], 'png');
        close;
    end
end
%% step2 以互相关系数为标准去重
% decompoMU = struct('row', {}, 'col', {}, 'pulse', {}, 'source', {}, 'CoV', {});
% decompo
% 假设 decompoMURaw 是您的结构体数组，包含 row, col, source, CoV 四个字段

% 第一步：提取所有source信号
sources = [decompoMURaw.source];
numSources = length(decompoMURaw);

% 第二步：计算带有时移的互相关系数矩阵
corrThreshold = 0.3;  % 相关系数阈值

% 初始化相关系数矩阵
crossCorrMatrix = zeros(numSources, numSources);

% 计算每对信号的最大互相关系数（考虑时移）
for i = 1:numSources
    for j = i+1:numSources  % 只计算上三角部分，避免重复
        % 计算互相关系数，考虑时移
        % [corrVals, ~] = xcorr(sources(:, i), sources(:, j), 'coeff');
        % 取绝对值最大值（考虑正负相关）
        % maxCorr = max(abs(corrVals));
        % crossCorrMatrix(i,j) = maxCorr;
        % crossCorrMatrix(j,i) = maxCorr;  % 对称矩阵

        [corrVals, ~] = corr(sources(:, i), sources(:, j));
        crossCorrMatrix(i,j) = abs(corrVals);
        crossCorrMatrix(j,i) = abs(corrVals);
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
            if decompoMURaw(i).CoV <= decompoMURaw(j).CoV
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
decompoMUFiltered = decompoMURaw(selectedIndices);

% 输出结果信息
fprintf('原始元素数量: %d\n', numSources);
fprintf('筛选后元素数量: %d\n', length(decompoMUFiltered));
fprintf('删除元素数量: %d\n', numSources - length(decompoMUFiltered));

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

%% 脉冲串匹配
load(['./Data/experiment/ICdata/R' Sub '/pulsesRef.mat']);
% 匹配容差为0~5ms
winSize = [0, 5]/1000*fsampu;
% lim=100ms，转换成样本点作为输入参数
lim = 100/1000*fsampu;

% decompoPulses = {};
% for i = 1:length(decompoMUFiltered)
%     decompoPulses{i} = decompoMUFiltered(i).pulse;
% end
% [matchResult, matchResultRaw] = PulseMatch(decompoPulses, pulsesRef, 0.3, 1000);

matchResultRaw = [];
for i = 1:length(pulsesRef)
    if length(pulsesRef{i}) < 150 || length(pulsesRef{i}) > 450
        disp(['i=' num2str(i) '参考脉冲序列无效']);
        continue;
    end
    for j = 1:length(decompoMUFiltered)
        xcor = fxcorr(pulsesRef{i}, decompoMUFiltered(j).pulse, lim);
        winSums = zeros(1, 2*lim+1-(winSize(2)-winSize(1)));
        for k = 1:2*lim+1-(winSize(2)-winSize(1))
            winSums(k) = sum(xcor(k:k+(winSize(2)-winSize(1))));
        end
        [maxMatchSum, maxIdx] = max(winSums);
        % 最佳匹配的时移
        Lag = maxIdx - winSize(1) - lim - 1;
        % 匹配上的脉冲个数，真阳性
        TP = maxMatchSum;
        % 假阳性
        FP = length(decompoMUFiltered(j).pulse) - maxMatchSum;
        % 假阴性
        FN = length(pulsesRef{i}) - maxMatchSum;
        % 匹配率
        matchRatio = TP / (TP + FN) * 100;
        % RoA
        rr = TP / (TP + FP + FN);
        % rr = RoA(decompoMUFiltered(j).pulse, pulsesRef{i}, 100, 5);

        matchResultRaw(end+1, :) = [i, j, Lag, TP, FP, FN, matchRatio, rr];
    end
end

matchResultRaw = array2table(matchResultRaw, 'VariableNames', {'ref', 'decomp', 'Lag', 'TP', 'FP', 'FN', 'match ratio', 'RoA'});

%% 计算筛选后MU的PNR
for sub = [3,4,5,7,10,11,12,14,15,16,17,18]
    % 导入数据
    load(['./Data/experiment/ICdata/R' num2str(sub) '/USCBSS_decomp_result.mat']);
    PNRs = [];
    Rows = [];
    for i = 1:length(decompoMUFiltered)
        IPT = decompoMUFiltered(i).source;
        pulse = decompoMUFiltered(i).pulse;
        % 10*log10(mean(tT(compInd).^2)/mean(tT(setdiff([1:length(tT)],compInd)).^2))
        Rows(i) = decompoMUFiltered(i).row;
        PNRs(i) = 10*log10( mean( IPT(pulse).^2 ) / mean( IPT(setdiff(1:length(IPT), pulse)).^2 ) );
    end
    save(['./Data/experiment/ICdata/R' num2str(sub) '/PNRs.mat'], 'PNRs', 'Rows');
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

%% 保存结果
save(['./Data/experiment/ICdata/R' Sub '/USCBSS_decomp_result.mat'], 'decompoMURaw', 'decompoMUFiltered', 'matchResultRaw');

% plotDecomps({pulsesRef{5},decompoMUFiltered(2).pulse}, [], 1000, 0, 0, []);

%% 绘制23*25个区域的估计源信号
for r = 1:23
    for c = 1:25
        tmpPulses = DecompoResults.decompo_pulses{r, c};
        tmpSources = DecompoResults.sources{r, c};
        tmpSourcesRaw = DecompoResults.sourceFirst{r, c};
        tmpCoV = DecompoResults.CoV{r, c};
        ax=figure;
        for mu = 1:length(tmpPulses)
            subplot(10,5,mu+floor((mu-1)/5)*5);
            plot(tmpSourcesRaw(:,mu));
            xlim([0,4000]);
            xticks(0:1000:4000);
            xticklabels(0:1:4);
            xlabel('t/s'); ylabel('amplitude');
            title(['mu=' num2str(mu) '一阶段迭代结果'])

            subplot(10,5,mu+ceil(mu/5)*5);
            plot(tmpSources(:,mu));
            xlim([0,4000]);
            xticks(0:1000:4000);
            xticklabels(0:1:4);
            xlabel('t/s'); ylabel('amplitude');
            title(['mu=' num2str(mu) '二阶段迭代结果'])
        end
        set(gcf,'unit','normalized','position',[0,0,1,1]);
        saveas(ax, ['./Results/R' Sub '/r' num2str(r) 'c' num2str(c)], 'png');
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
[cc, ~] = xcorr(decompoSources, 'coeff');
% 转成三维数组
cc = reshape(cc', size(decompoSources, 2), size(decompoSources, 2), []);
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

%% 绘制筛选后的MU
close all;
fsampu = 1000;
numMU = length(decompoMURaw);
L = 30000;
f = (0:1:L/2) * fsampu / L;
for muii = 1:5%numMU
    % 傅里叶变换
    tmp = abs(fft(decompoMURaw(muii).source)/L);
    sFFT = tmp(1:L/2+1);
    sFFT(2:end-1) = 2*sFFT(2:end-1);
    energyRatio = sum(sFFT(1:100)) / sum(sFFT);

    figure;
    subplot(4,1,1);
    scatter(decompoMURaw(muii).pulse, ones(length(decompoMURaw(muii).pulse)), 1000, "black", '|');
    xlim([0, L])
    title('pulse train')

    subplot(4,1,2);
    plot(decompoMURaw(muii).sourceRaw);
    % xlim([10000, 12000])
    % ylim([-5, 5])
    title('estimated source step1')

    subplot(4,1,3);
    plot(decompoMURaw(muii).source);
    % xlim([10000, 12000])
    % ylim([-2, 2])
    title('estimated source step2')

    subplot(4,1,4);
    plot(f, sFFT);
    title(['ratio='  num2str(energyRatio*100) '%']);
    sgtitle(['MU #' num2str(muii)]);
end

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