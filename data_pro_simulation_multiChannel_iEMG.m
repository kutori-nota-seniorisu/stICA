%% 仿真数据处理
clear; clc; close all;
addpath('./Func')
weight = 0.9;
model = 'NC_NMFm';
% 导入数据
for set_i = 1:10
    load(['F:/EEEMG/stICA/Data/simulation/datasets' num2str(set_i) '/UV_compo12_' model num2str(weight) '.mat']);

    % 分解后的空间成分图像
    figure;
    for i = 1:size(S,2)
        subplot(2,6,i);
        imagesc(reshape(S(:,i),400,128));
        colorbar;
        title(['Space #' num2str(i)]);
    end
    sgtitle('分解后的空间成分');

    % 源空间成分图像
    figure;
    for i = 1:size(Xs,2)
        subplot(2,6,i);
        imagesc(reshape(Xs(:,i),400,128));
        colorbar;
        title(['Space #' num2str(i)]);
    end
    sgtitle('源空间成分');

    % spike train extraction
    flag_extract = 'CBSS';%'smooth','BPM','HWM','CBSS'
    decompo_pulses = {};

    if strcmp(flag_extract,'smooth')
        disp('平滑处理');
        % move median
        T_smooth = smoothdata(T, "movmedian", 8);

        % z-score标准化
        T_mean = mean(T_smooth);
        T_std = std(T_smooth);
        T_norm = (T_smooth-T_mean)./T_std;

        % 采样率2000Hz，故1ms=2样本点。设置MinPeakDistance时需要注意ms与样本点的转换
        for i = 1:size(T_norm,2)
            [~, locs] = findpeaks(T_norm(:,i),'MinPeakDistance',90);%,'MinPeakHeight',0.8);
            decompo_pulses{end+1} = locs';
        end

    elseif strcmp(flag_extract,'BPM')
        disp('BPM带通滤波');
        % 5-25Hz带通滤波
        fsampu = 2000;
        [B,A] = butter(4,[5,25]/fsampu*2);
        T_filter = filtfilt(B,A,T);

        % z-score标准化
        T_mean = mean(T_filter);
        T_std = std(T_filter);
        T_norm = (T_filter-T_mean)./T_std;

        for i = 1:size(T_norm,2)
            [~, locs] = findpeaks(-T_norm(:,i),'MinPeakDistance',50,'MinPeakHeight',0.5);
            decompo_pulses{end+1} = locs';
        end

    elseif strcmp(flag_extract,'HWM')
        disp('HWM哈尔小波变换');
        % HWM哈尔小波变换
        fsampu = 2000;
        scale = fsampu/45;
        for i = 1:size(T,2)
            tmp = cwt(T(:,i),scale,'haar');
            T_HWM(:,i) = tmp';
        end

        % z-score标准化
        T_mean = mean(T_HWM);
        T_std = std(T_HWM);
        T_norm = (T_HWM-T_mean)./T_std;

        for i = 1:size(T_norm,2)
            [~, locs] = findpeaks(T_norm(:,i),'MinPeakDistance',50*2,'MinPeakHeight',0.9);
            decompo_pulses{end+1} = locs';
        end
    elseif strcmp(flag_extract,'CBSS')
        disp('卷积盲源分离');
        % 拓展因子
        exFactor = 20;
        % nMAD
        nMAD = 4;
        % MPD(ms)
        MPD = 20;

        for i = 1:12
            % 是否要在T的前面加上负号？我认为不需要，迭代结果不受符号的影响
            [source, ~, PT, ~, w_new] = blindDeconvPeakFinding(T(:, i)', exFactor, nMAD, MPD*2);
            decompo_pulses{i} = PT;
            T_new(:, i) = source';
            W(:, i) = w_new;
        end
    end

    % 绘制某一成分
    % figure;
    % subplot(2,1,1);
    % plot(-T(:,1));
    % title('原始分量')
    % subplot(2,1,2);
    % plot(-T_norm(:,1));
    % title('滤波&标准化处理');

    % 导入simulation iEMG
    load(['F:/EEEMG/stICA/Data/simulation/MU_time_response/TimeCompoDatasets' num2str(set_i) '/ipulses.mat']);

    % 空间成分匹配
    % 计算相关系数并选出最大匹配
    corr_matrix_space = corr(S,Xs);
    corr_matrix_space(:,end-1:end) = [];
    [max_values_space, indices_space] = max(abs(corr_matrix_space),[],2);

    % 在重复匹配中找最大匹配
    matchresult_space_raw = [(1:12)',indices_space, max_values_space];
    matchresult_space = matchresult_space_raw;
    for mu = 1:size(S,2)
        tmpInd = find(matchresult_space(:,2) == mu);
        if length(tmpInd) > 1
            [~,tmpInd2] = max(matchresult_space(tmpInd,3));
            deleteInd = setdiff(1:length(tmpInd),tmpInd2);
            matchresult_space(tmpInd(deleteInd),:) = [];
        end
    end
    matchresult_space_raw = array2table(matchresult_space_raw, 'VariableNames', {'decomp', 'ref', 'space'});
    matchresult_space = array2table(matchresult_space, 'VariableNames', {'decomp', 'ref', 'space'});

    % 时间成分匹配
    matchresult_time_raw = [];
    % 计算RoA
    for i = 1:length(decompo_pulses)
        for j = 1:length(ipulses)
            [Array1, Array2] = meshgrid(decompo_pulses{i}, ipulses{j});
            diff_values = Array1 - Array2;
            valid_elements = diff_values <= 60 & diff_values >= 0;
            count = sum(valid_elements(:));
            r = count/(length(decompo_pulses{i})+length(ipulses{j})-count);
            if r > 1
                r = 1;
            end
            spike_ROA_matrix(i,j) = r;
            matchresult_time_raw(end+1,:) = [i, j, r];
        end
    end

    % 选出最大的匹配
    matchresult_time = matchresult_time_raw;
    for mu = 1:length(decompo_pulses)
        tmpInd = find(matchresult_time(:,1) == mu);
        if length(tmpInd) > 1
            [~, tmpInd2] = max(matchresult_time(tmpInd,3));
            deleteInd = setdiff(1:length(tmpInd), tmpInd2);
            matchresult_time(tmpInd(deleteInd), :) = [];
        end
    end
    for mu = 1:length(ipulses)
        tmpInd = find(matchresult_time(:,2) == mu);
        if length(tmpInd) > 1
            [~, tmpInd2] = max(matchresult_time(tmpInd,3));
            deleteInd = setdiff(1:length(tmpInd), tmpInd2);
            matchresult_time(tmpInd(deleteInd), :) = [];
        end
    end
    matchresult_time_raw = array2table(matchresult_time_raw, 'VariableNames', {'decomp', 'ref', 'time'});
    matchresult_time = array2table(matchresult_time, 'VariableNames', {'decomp', 'ref', 'time'});

    matchresult_final = innerjoin(matchresult_space, matchresult_time, 'Keys', {'decomp', 'ref'});

    matchresult_space_all{set_i} = matchresult_space;
    matchresult_time_all{set_i} = matchresult_time;
    matchresult_final_all{set_i} = matchresult_final;
    close all;
end
% 保存匹配结果
save(['F:/EEEMG/stICA/Results/compo12_' model num2str(weight) '_' flag_extract '_10sets.mat'],'matchresult_space_all','matchresult_time_all','matchresult_final_all');

