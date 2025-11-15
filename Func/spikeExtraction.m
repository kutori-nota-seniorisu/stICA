function [spikeInd,SIL,C] = spikeExtraction(IPT,fsamp)
tmpInd = isnan(IPT);
if sum(tmpInd)==length(IPT)
    spikeInd = [];
    SIL = 0;
    C = [];
    return;
end
PulseDistance = round(fsamp/100);
% PulseDistance = 40;

[idx,rawC,sumd,D] = kmeans(IPT,2,'start',[mean(IPT);0]);
% [idx,rawC,sumd,D] = kmeans(IPT,2);

% D = sqrt(D);

Ind1 = find(idx==1);
Ind2 = find(idx==2);
% sum_intraCluster = sum(sumd);
% sum_interCluster = sum(D(Ind2,1))+sum(D(Ind1,2));
% SIL = abs(sum_intraCluster-sum_interCluster)/max(sum_intraCluster,sum_interCluster);


%移除距离过近的尖峰，以避免重复检测
if rawC(1)<rawC(2)
    spikeInd = remRepeatedInd(IPT,Ind2,PulseDistance);
%     silInd = setdiff([1:length(IPT)],spikeInd);
else
    spikeInd = remRepeatedInd(IPT,Ind1,PulseDistance);
%     silInd = spikeInd;
end
% C(1) = max(rawC);
% C(2) = min(rawC);

%计算两个簇的均值 C，分别对应于尖峰和非尖峰信号。
C(1) = mean(IPT(spikeInd));
C(2) = mean(IPT(setdiff([1:length(IPT)],spikeInd)));

% mean(IPT(Ind1))

%计算距离矩阵 D，其中每一行代表信号中的一个点与两个簇中心的距离
D(:,1) = abs(IPT-C(1));
D(:,2) = abs(IPT-C(2));

flag = ones(length(idx),1);
flag(spikeInd) = -1;

%计算分离系数（轮廓指数） SIL，这是通过比较每个点到两个簇中心的距离，并取平均值得到的。
SIL = mean(((D(:,1)-D(:,2)).*flag)./max(D,[],2));

% mean(IPT(Ind2))