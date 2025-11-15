function [minInd,newx,maxInd] = findHillLim(y,x,Lim)
% find the first peak point near x. newx is the peak point, minInd is the left
% valley, maxInd is the right valley. 

if nargin < 3;
    Lim = 100;
end 

tmp = diff(y);

if x > length(tmp) || abs(tmp(x)) == 0
    minInd = x;
    maxInd = x;
    newx = x;
    return ;
end 

dx = tmp(x)/abs(tmp(x));
% if temp(x)<0, searching peak towards left, otherwise, searching towards
% right 

newx = x;
while newx > 1 && newx < length(y) && y(newx+dx) > y(newx)
    newx = newx+dx;
end 

minInd = newx;
maxInd = newx;

while minInd > 1 && y(minInd-1) < y(minInd) && newx-minInd < Lim
    minInd = minInd-1;
end 
while maxInd < length(y)-1 && y(maxInd+1) < y(maxInd) && maxInd-newx < Lim
    maxInd = maxInd+1;
end 