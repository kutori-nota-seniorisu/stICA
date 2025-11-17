% 对超声信号进行STA
clear; clc; close all;
addpath('./Func');

win = [-50, 50];
fsampu = 1000;

for sub = [16]
    disp(['Sub=' num2str(sub)]);
    % 导入超声速度图像
    tviFile = ['./Data/experiment/ICdata/R' num2str(sub) '/v_2d_all.mat'];
    load(tviFile);
    TVIData = cat(3, zeros(119, 128, 2), v_2d_all);
    % filter the TVI data
    TVIDataFilter = TVIData;
    % 时间5-100Hz带通滤波
    [Be2, Ae2] = butter(4, [5, 100]/fsampu*2);
    for r = 1:size(TVIDataFilter, 1)
        for c = 1:size(TVIDataFilter, 2)
            tmp = squeeze(TVIDataFilter(r, c, :));
            tmp = filtfilt(Be2, Ae2, tmp);
            TVIDataFilter(r, c, :) = tmp;
        end
    end

    [M, N, L] = size(TVIDataFilter);
    % 窗口大小
    Row = 10; Col = 10;
    % 窗口移动距离
    dRow = 5; dCol = 5;
    
    load(['./Data/experiment/ICdata/R' num2str(sub) '/USCBSS_DecompResult_MPD60.mat']);

    for mu = 1:length(decompoMUFiltered.MU)
        tmpPulses = decompoMUFiltered.Pulse{mu};
        if isempty(tmpPulses)
            continue;
        end
        frameSTA = zeros(0);
        varSTA = zeros(0);
        for n = win(1)+1:win(2)
            tmpInd = tmpPulses + n;
            tmpInd(tmpInd <= 0) = [];
            tmpInd(tmpInd >= 30000) = [];
            tmpTVI = TVIDataFilter(:, :, tmpInd);
            % 存储STA图像
            frameSTA(:, :, n-win(1)) = mean(tmpTVI, 3);
            % 存储方差图像
            varSTA(:, :, n-win(1)) = var(tmpTVI, 0, 3);
        end

        r = decompoMUFiltered.Row(mu);
        c = decompoMUFiltered.Col(mu);
        winRow = (1:Row)+(r-1)*dRow;
        winCol = (1:Col)+(c-1)*dCol;
        winRow(winRow>M) = [];
        winCol(winCol>N) = [];

        for r = 1:length(winRow)
            for c = 1:length(winCol)
                frameSTAROI{r, c} = squeeze(frameSTA(winRow(r), winCol(c), :));
            end
        end

        % [rn, cn, ~] = size(frameSTA);
        % tmpSTA = cell(0);
        % for r = 1:rn
        %     for c = 1:cn
        %         tmpSTA{r, c} = squeeze(frameSTA(r, c, :));
        %     end
        % end
        
        figure;
        ax = subplot('Position', [0.05, 0.05, 0.9, 0.9]);
        [~,maxAmp,maxPP,pos]=plotArrayPotential(frameSTAROI, 1, 1, ax);
    end
end