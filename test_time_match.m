% 测试不同的时间成分匹配方法。目前采用的是方法3，另外CBSS的方法也实现了，参数还需要优化。

%% 时间成分匹配1。用T和Xt的相关系数做比较，只能作为仿真数据的辅助佐证。
load('Data/simulation/datasets1/UV_compo12_NC_NMFm.mat');
% 绘制分解后时间成分图像
figure;
for i = 1:size(T,2)
    subplot(2,6,i);
    plot(T(:,i));
    title(['Time #' num2str(i)]);
end
sgtitle('分解后的时间成分')
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

% 绘制源时间成分图像
figure;
for i = 1:size(Xt,2)
    subplot(2,6,i);
    plot(Xt(:,i));
    title(['Time #' num2str(i)]);
end
sgtitle('源时间成分')
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

corr_matrix_time = corr(T,Xt);
[max_values_time, indices_time] = max(abs(corr_matrix_time),[],2);

matchresult_time_raw = [(1:12)',indices_time, max_values_time];
matchresult_time = matchresult_time_raw;
for mu = 1:size(T,2)
    tmpInd = find(matchresult_time(:,2) == mu);
    if length(tmpInd) > 1
        [~,tmpInd2] = max(matchresult_time(tmpInd,3));
        deleteInd = setdiff(1:length(tmpInd),tmpInd2);
        matchresult_time(tmpInd(deleteInd),:) = [];
    end
end

figure;
imagesc(abs(corr_matrix_time));
colorbar;
colormap gray;

%% 时间成分匹配2。提取脉冲串并计算RoA，调用了RoA函数。
load('Data/simulation/datasets1/UV_compo12_NC_NMFm.mat');
load('Data/simulation/MU_time_response/TimeCompoDatasets1/ipulses.mat');
fsampu = 2000;
[B,A] = butter(4,[5,25]/fsampu*2);
T_filter = filtfilt(B,A,T);

% z-score标准化
T_mean = mean(T_filter);
T_std = std(T_filter);
T_norm = (T_filter-T_mean)./T_std;

matchresult_time_raw = [];
decompo_pulses = {};
for i = 1:size(T_norm,2)
    [~, locs] = findpeaks(-T_norm(:,i), 'MinPeakDistance', 50, 'MinPeakHeight', 0.5);
    decompo_pulses{i} = locs';
end
for i = 1:length(decompo_pulses)
    for j = 1:length(ipulses)
        [r, lag] = RoA(decompo_pulses{i}, ipulses{j}, 100, 30);
        matchresult_time_raw(end+1, :) = [i, j, r, lag];
    end
end

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

%% 时间成分匹配3。以识别出的脉冲串直接进行匹配，由于空间匹配的相关系数是负值，因而在提取脉冲串时需要添加负号。
load('Data/simulation/datasets1/UV_compo12_NC_NMFm0.9.mat');
load('Data/simulation/MU_time_response/TimeCompoDatasets1/ipulses.mat');
fsampu = 2000;
[B,A] = butter(4,[5,25]/fsampu*2);
T_filter = filtfilt(B,A,T);

% z-score标准化
T_mean = mean(T_filter);
T_std = std(T_filter);
T_norm = (T_filter-T_mean)./T_std;

figure;
for i = 1:size(Xt,2)
    subplot(2,6,i);
    plot(Xt(:,i));
       title(['Time #' num2str(i)]);
end 
sgtitle('源时间成分');
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

% 把处理后的时间成分画出来
figure;
for i = 1:size(T,2)
    subplot(2,6,i);
    plot(-T_norm(:,i));
    title(['Time #' num2str(i)]);
end
sgtitle('滤波&标准化 时间成分');
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

matchresult_time_raw = [];
decompo_pulses = {};
% 提取脉冲串
for i = 1:size(T_norm,2)
    [~, locs] = findpeaks(-T_norm(:,i),'MinPeakDistance',50,'MinPeakHeight',0.5);
    decompo_pulses{i} = locs';
end
% 计算RoA
for i = 1:length(decompo_pulses)
    for j = 1:length(ipulses)
        [Array1, Array2] = meshgrid(decompo_pulses{i}, ipulses{j});
        diff_values = Array1 - Array2;
        valid_elements = diff_values <= 45 & diff_values >= 5;
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

%% 时间成分匹配4。以识别出的脉冲串直接进行匹配，采用HWM方法，不需要加负号
load('Data/simulation/datasets1/UV_compo12_NC_NMFm.mat');
load('Data/simulation/MU_time_response/TimeCompoDatasets1/ipulses.mat');
fsampu = 2000;
scale = fsampu/45;
for i = 1:size(T,2)
    tmp = cwt(T(:,i), scale, 'haar');
    T_HWM(:,i) = tmp';
end

% z-score标准化
T_mean = mean(T_HWM);
T_std = std(T_HWM);
T_norm = (T_HWM-T_mean)./T_std;

% 把处理后的时间成分画出来
figure;
for i = 1:size(T,2)
    subplot(2,6,i);
    plot(-T_norm(:,i));
    title(['Time #' num2str(i)]);
end
sgtitle('滤波&标准化 时间成分');
set(gcf,'unit','normalized','position',[0.05,0.5,0.9,0.3]);

matchresult_time_raw = [];
decompo_pulses = {};
% 提取脉冲串
for i = 1:size(T_norm,2)
    [~, locs] = findpeaks(T_norm(:,i),'MinPeakDistance',100,'MinPeakHeight',0.9);
    decompo_pulses{i} = locs';
end
% 计算RoA
for i = 1:length(decompo_pulses)
    for j = 1:length(ipulses)
        [Array1, Array2] = meshgrid(decompo_pulses{i}, ipulses{j});
        diff_values = abs(Array1 - Array2);
        valid_elements = diff_values <= 60;
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

