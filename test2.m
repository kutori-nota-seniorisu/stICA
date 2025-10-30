%% 统计在RoA计算中，各脉冲间隔的最小值，然后绘制箱线图进行统计
clear; clc; close all;
for setset = 1:10
    datasets_num = num2str(setset);
    load(['./Results/datasets' datasets_num '_result.mat']);
    load(['Data/simulation/MU_time_response/TimeCompoDatasets' datasets_num '/ipulses.mat']);

    for i = 1:size(matchresult_time, 1)
        decomp = matchresult_time.decomp(i);
        ref = matchresult_time.ref(i);
        [Array1, Array2] = meshgrid(decompoPulseAll{decomp}, ipulses{ref});
        diff_values = Array1 - Array2;
        [minAbsVals, minInd] = min(abs(diff_values), [], 2);
        for j = 1:length(minInd)
            minVals(j+(i-1)*length(minInd)) = diff_values(j, minInd(j));
        end
        % minVals{i} = tmp_minVals;
    end

    save(['./Results/datasets' datasets_num '_RoA_Interval.mat'], 'minVals');
end

%%
figure;
for setset = 1:10
    datasets_num = num2str(setset);
    load(['./Results/datasets' datasets_num '_RoA_Interval.mat']);
    boxplot(minVals, 'Positions', setset, 'Colors', 'r');
    [lower_adj(setset), upper_adj(setset), distance(setset), IQR(setset)] = calculate_adjacent_distance(minVals);
    hold on;
end
plot(1:10, IQR, 'g-*');
plot(1:10, distance, 'b-*')
xticks(1:10)
xticklabels(1:10);
ylim([0,40])
ylabel('RoA Interval');
legend('IQR', '上下邻距离')
set(gcf,'unit','normalized','position',[0.05,0.1,0.9,0.3]);

function [lower_adj, upper_adj, distance, IQR] = calculate_adjacent_distance(data)
% 计算箱线图上下邻距离
% 输入: data - 数据向量
% 输出: lower_adj - 下邻, upper_adj - 上邻, distance - 距离

Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
IQR = Q3 - Q1;

% 理论边界
theoretical_lower = Q1 - 1.5 * IQR;
theoretical_upper = Q3 + 1.5 * IQR;

data(data<theoretical_lower) = [];
data(data>theoretical_upper) = [];

% 实际边界（考虑数据范围）
valid_data = data(isfinite(data)); % 去除无穷大和NaN
lower_adj = max(min(valid_data), theoretical_lower);
upper_adj = min(max(valid_data), theoretical_upper);

distance = upper_adj - lower_adj;
end

