function [PulseStat,SourceID,Lag,Sensitivity,Misses,FalseAlarms,Specificity] = testSinResults(sInd,tInd,d,displayMode,cost)
% [PulseStat,SourceID,Lag,Sensitivity,Misses,FalseAlarms,Specificity] = testSinResutls2(sInd,tInd,d,displayMode,cost)
% Calculates different matching statistics between the tested pulse sequence tInd and reference pulse sequences sInd.
% The method returns the ID of the sInd sequence exhibiting the best match with tInd sequence.
% 计算两个放电串之间的统计学指标
% - 灵敏度Sensitivity
% - 假阴性率Misses
% - 误报率FalseAlarms
% - 特异性Specificity
% - 真假阴阳性，存储在PulseStat内
% - 时延Lag

% INPUTS:
%   sInd - cell structure containing the discharge times [in samples] of reference pulse sequences
%   tInd - ischarge times [in samples] of tested pulse sequence
%   d - allowed pulse-to-pulse missmatch [in samples] between two compared pulse sequences; set d=0 for exact matching of pulses;
%   displayMode - use displayMode>0 for debug mode (displays the results of comparisson with all tInd pulse seguences)
%   cost - (optional) - heuristic score of tested pulse sequence (used only in debug mode)
% OUTPUTS:
%   SourceID  - best matched source number (ID in sInd)
%   PulseStat  - cell structure containing all the statistics of tInd  w.r.t. the best mached delayed repetition of sInd
%   Lag  - delay of best matched source sInd{SourceID} w.r.t. tInd
%   TP - true  positive statistics (w.r.t. the best matched pulse sequence)
%   FP - false positive statistics (w.r.t. the best matched pulse sequence)
%   FN - false negative statistics (w.r.t. the best matched pulse sequence)
%   TN - true  negative statistics (w.r.t. the best matched pulse sequence)



% sInd tends to be a cell with length = 1
% 直接设置偏移极限为100个样本，也许可以设置为可以更改的参数。
lim = 100;

if isempty(sInd) || isempty(tInd)
    PulseStat.TP = 0;
    PulseStat.FP = 0;
    PulseStat.FN = 1;
    SourceID = 0;
    Lag = 0;
    Sensitivity = 0;
    Misses = 1;
    FalseAlarms = 1;
    Specificity = 0;
    return ;
end

% for k = 1:length(sInd)
%     a = fxcorr(tInd,sInd{k},lim);
%     [~,Lag(k)] = max(a);
%     m(k) = min([sum(a(max(1,Lag(k)-d):min(length(a),Lag(k)+d))),length(sInd{k}),length(tInd)]);
% end

% 统计不同时延情况下，两个放电串具有相同脉冲的个数
% 输入的放电串是01序列（？）
a = fxcorr(tInd,sInd,lim);
% 最佳匹配对应的平移量为Lag
[~,Lag] = max(a);
% missmatch范围内的相同脉冲总和，并且与放电串长度对比。
% MaxPeak为两个放电串重合脉冲的总个数。
MaxPeak = min([sum(a(max(1,Lag-d):min(length(a),Lag+d))), length(sInd), length(tInd)]);

% 完全没有重合
if MaxPeak == 0
    PulseStat.TP = 0;
    PulseStat.FP = 0;
    PulseStat.FN = 1;
    SourceID = 0;
    Lag = 0;
    Sensitivity = 0;
    Misses = 1;
    FalseAlarms = 1;
    Specificity = 0;
    return ;
end

% [MaxPeak,~] = max(m);
Lag = Lag-lim-1;    % calculate the actual lag, + means PT1 proceeds PT2, - means PT2 proceeds PT1

% Ind = find(m==MaxPeak);
% 找到平移距离最小的平移量
[~,SourceID] = min(abs(Lag));
% SourceID = Ind(SourceID);
% Lag = Lag(SourceID);

% sInd means the standard pulse train while tInd means the test pulse train
% TP - true positive 真阳性，实际为阳性，预测结果为阳性。
% 这里定义为 解码与参考重合的脉冲个数，也即 MaxPeak
PulseStat.TP = MaxPeak; % the number of correctly identified discharges
% FP - false positive 假阳性，实际为阴性，预测结果为阳性。
% 这里定义为 在解码但不在参考的脉冲个数，也即 length(tInd)-MaxPeak
PulseStat.FP = length(tInd)-MaxPeak; % the number of misplaced discharges
% FN - false negative 假阴性，实际为阳性，预测结果为阴性
% 这里定义为 在参考但不在解码的脉冲个数，也即 length(sInd)-MaxPeak
PulseStat.FN = length(sInd)-MaxPeak; % the number of unidentified discharges
% TN - true negative 真阴性，实际为阴性，预测结果为阴性
% 对这个公式存疑
PulseStat.TN = max([sInd,tInd])-length(sInd)-PulseStat.FN;

% 灵敏度定义为 实际阳性中正确预测出阳性的概率
Sensitivity = PulseStat.TP/(PulseStat.TP+PulseStat.FN);
% 假阴性率定义为 实际阳性中预测出阴性的概率
Misses = PulseStat.FN/(PulseStat.TP+PulseStat.FN);
FalseAlarms = PulseStat.FP/(PulseStat.FP+PulseStat.FN);
% FalseAlarms formula has been corrected depending on the [Ref: Accurate identification of motor unit discharge patterns
% from high-density surface EMG and validation with a novel signal-based performance metric]
% 特异性定义为 实际阴性中正确预测出阴性的概率
Specificity = PulseStat.TN/(PulseStat.FP+PulseStat.TN);

if displayMode > 0
    if nargin < 5
        disp(['   SourceID = ',num2str(SourceID),'   lag = ',num2str(Lag),'   Sensitivity = ',num2str(Sensitivity),' ( ',num2str(PulseStat.TP),' / ',num2str(PulseStat.FP),' )   Misses = ',num2str(Misses),' ( ',num2str(PulseStat.TN),' / ',num2str(PulseStat.FN),' )']);
    else
        disp(['   SourceID = ',num2str(SourceID),'   lag = ',num2str(Lag),'   Sensitivity = ',num2str(Sensitivity),' ( ',num2str(PulseStat.TP),' / ',num2str(PulseStat.FP),' )   Misses = ',num2str(Misses),' ( ',num2str(PulseStat.TN),' / ',num2str(PulseStat.FN),' );  IPI cost = ',num2str(cost)]);
    end
end


