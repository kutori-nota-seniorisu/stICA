function [ y ] = normalization( x,ymin,ymax )
%NORMALIZATION 将数据x归一化到任意区间[ymin,ymax]范围的方法
%   输入参数x：需要被归一化的数据
%   输入参数ymin：归一化的区间[ymin,ymax]下限
%   输入参数ymax：归一化的区间[ymin,ymax]上限
%   输出参数y：归一化到区间[ymin,ymax]的数据
xmax=max(x);%计算最大值
xmin=min(x);%计算最小值
y = (ymax-ymin)*(x-xmin)/(xmax-xmin) + ymin;
end