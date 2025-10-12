function [dS,dSens,dLag] = sortIPTcorr(sInd,tInd,dIPI,Ind,MatchTrsh)
% sorts the pulse sequences in sInd according to their match with pulse
% sequence in tInd. Parameters Ind and MatchTrsh are for debugging purposes
% only and are optional.
% 寻找sInd中与tInd匹配的脉冲串。
% OUTPUTS:
% dS - array of binary values (0 and 1); dS(k)==1 indicates the match
% between the pulse sequence in tInd and the pulse sequence in sInd{k}
% dSens - array of Sensitivity values; dSens(k) contains the sensitivity of
% match between the pulse sequence in tInd and the pulse sequence in sInd{k}
% dLag - array of Lags; dLag(k) contains the ralative delay of pulse
% sequence sInd{k} w.r.t. the pulse sequence tInd

if (nargin<4)
    Ind = 1:length(sInd);
end
if (nargin <5)
    MatchTrsh = 0.5;
end

dS = zeros(length(sInd),1);
dSens = zeros(length(sInd),1);
dLag = zeros(length(sInd),1);
for k = 1:length(sInd)
    % 正算+反算，计算两种情况下的参数
    [PulseStat,SourceID,Lag,Sens,Miss,~,~] = testSinResults(sInd{k},tInd,dIPI,0);
    [PulseStat2,SourceID2,Lag2,Sens2,Miss2,~,~] = testSinResults(tInd,sInd{k},dIPI,0);
    % [max(灵敏度)>匹配阈值 且 min(预测结果为阳性的个数)>50] 或者 max(灵敏度)>0.5
    % 满足这个条件的，才认为两个放电串匹配
    % 前者的灵敏度阈值低，但是要求预测阳性结果的个数满足条件，后者的灵敏度要求高。为什么是50和0.5？经验参数？
    if ((max(Sens,Sens2)>MatchTrsh && min(PulseStat.TP+PulseStat.FP,PulseStat2.TP+PulseStat2.FP)>50) || (max(Sens,Sens2)>0.5))
        if Sens > Sens2
            % disp(['',num2str(Ind(k)),':   i = ',num2str(SourceID),'   lag = ',num2str(Lag),'   Sensitivity = ',num2str(Sens),' ( ',num2str(PulseStat.TP),' / ',num2str(PulseStat.TP+PulseStat.FN), ' )   Misses = ',num2str( Miss ),' ( ',num2str(PulseStat.FN),' / ',num2str(PulseStat.TP+PulseStat.FP),' )']);
            dS(k,1) = 1;
            dSens(k,1) = Sens;
            dLag(k) = Lag;
        else
            % disp(['',num2str(Ind(k)),':   i = ',num2str(SourceID2),'   lag = ',num2str(Lag2),'   Sensitivity = ',num2str(Sens2),' ( ',num2str(PulseStat2.TP),' / ',num2str(PulseStat2.TP+PulseStat2.FN),' )   Misses = ',num2str(Miss2),' ( ',num2str(PulseStat2.FN),' / ',num2str(PulseStat2.TP+PulseStat2.FP),' )']);
            dS(k,1) = 1;
            dSens(k,1) = Sens2;
            dLag(k) = Lag2;
        end
    end
end

