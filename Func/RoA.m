function [r,Lag] = RoA(Ind1,Ind2,Lim,dIPI)
% -- Created by CC -- Email: cedric_c@126.com %
% 计算两个放电串之间的一致率RoA
% INPUTS:
%   Ind1:第一个放电串，通常为解码得到的结果。
%   Ind2:第二个放电串，通常为参考数据。
%   注：Ind1与Ind2应该为数组而非元胞，且均为以样本格式而非时间格式表示的放电串。具体而言，时间格式表示的放电时刻*采样率=样本格式表示的放电时刻。
%   Lim:考虑到放电时刻会出现整体的偏移，让放电串在-Lim到Lim范围先平移再寻找最大的匹配情况。以样本格式表示。
%   dIPI:RoA计算时的容忍程度
% OUTPUTS:
%   r:RoA
%   Lag:平移量，或者称为时延。

xcor = zeros(2*Lim+1,1);
for k = -Lim:Lim
    % intersect(A,B)返回A与B的交集，去除重复项，并且排序
    % step1:将一个放电串在容差范围内整体平移。
    % step2:计算两者的交集，并且计算交集的大小。这个值表示在当前平移条件下，两个放电串之间相同脉冲的个数。
    % step3:容差的范围为±Lim，令k遍历这个范围。xcor中的元素表示当平移量为k时，两个放电串具有相同脉冲的个数。
    xcor(k+Lim+1) = length(intersect(Ind1,Ind2+k));
end

% 找到xcor中最大值所在的索引。也许可以认为，这个a指向了一个最佳匹配序列。
a = find(xcor == max(xcor));
% 如果有不止一个最大值，那么选择距离Lim最近的那个索引（为什么要这么选？）
if length(a)>1
    tmp = abs(a-Lim);
    [~,tmpind] = min(tmp);
    a = a(tmpind);
end

% 对a-dIPI到a+dIPI范围内相同脉冲的个数求和。max与min函数用于防止a-dIPI和a+dIPI超出数组索引。
A = sum( xcor( max(1,a-dIPI) : min(a+dIPI,length(xcor)) ) );
% A表示两个放电串共有脉冲的个数，这个值不可能超过放电串的长度。
if A>min(length(Ind1),length(Ind2))
    A = min(length(Ind1),length(Ind2));
end

% 平移量，小于零表示向左平移，大于零表示向右平移，平移量为Lag的绝对值。
% 我觉得应该为Lag=a-Lim-1。这个想法和testSinResults中是一样的，所以改了一下。
% Lag = a-Lim;
Lag = a-Lim-1;

% RoA计算公式为：C/(C+A+B)
% 其中C表示两个放电串共有脉冲的个数，A表示只在放电串1中的脉冲个数，B表示只在放电串2中的脉冲个数。
r = A/(length(Ind1)+length(Ind2)-A);
end