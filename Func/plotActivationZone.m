function [uufilter,actArea] = plotActivationZone(muap,uuimage,pos,IED,thresh,flag_plot)

if nargin<6
    flag_plot = 1;
end

% threshVal = 0;
% uufilter = cell(0);
if isempty(uuimage)
    threshVal = 0;
else
    threshVal = thresh*max(abs(uuimage(:)));
end

if flag_plot
    % res_width = 0.9/size(muap,2)/IED;
    % res_height = 0.9/(size(muap,1)*IED+60);
    % ax = subplot('Position', [0.05,50*res_height+0.1,0.9,res_height*size(muapsAll,1)*IED]);
    ax = subplot('Position', [0.05,0.5,0.9,0.45]);
    plotArrayPotential(muap,1,1,ax);
    subplot('Position', [0.05,0.05,0.9,0.4]);
    hold on;
end

tmp = uuimage;
for na = 1:2
    if na == 1
        tmpInd = find(tmp>threshVal);
    else
        tmpInd = find(tmp<-threshVal);
    end
    actArea{na} = [rem(tmpInd,size(tmp,1)),ceil(tmpInd/size(tmp,1))];
    actArea{na}(find(actArea{na}==0)) = size(tmp,1);
end

tmp(abs(tmp)<=threshVal) = 0;
uufilter = tmp;
if flag_plot
    imagesc([0,size(uufilter,2)],[size(uufilter,1),0],abs(uufilter));
    colormap("jet");
    colorbar;
    axis('tight','off');
end