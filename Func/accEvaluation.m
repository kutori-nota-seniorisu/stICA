function [Sen,FA,Pre,Spe,Acc] = accEvaluation(ST_decomp,ST_true,dIPI,lim)
% -- Created by CC -- %
% 计算灵敏度Sen、误报率FA、精确率Pre、特异性Spe、准确率Acc
xcor = zeros(2*lim+1,1);
for k = -lim:lim
    xcor(k+lim+1) = length(intersect(ST_decomp,ST_true+k));
end 

[~,Lag] = max(xcor);
TP = min([sum(xcor(max(1,Lag-dIPI):min(length(xcor),Lag+dIPI))),length(ST_decomp),length(ST_true)]);
FP = length(ST_decomp)-TP;
FN = length(ST_true)-TP;
TN = max([ST_decomp,ST_true])-length(ST_decomp)-FN;

Sen = TP/(TP+FN);
Pre = TP/(TP+FP);
FA = FP/(FP+FN);
Spe = TN/(TN+FP);
Acc = (TP+TN)/(TP+TN+FP+FN);
