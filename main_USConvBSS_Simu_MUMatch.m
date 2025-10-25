%% 适用于USConvBSS的RoA计算与MU匹配
load(['Data/simulation/MU_time_response/TimeCompoDatasets' datasets_num '/ipulses.mat']);

matchresult_time_raw = [];
% 计算RoA
for i = 1:length(decompoPulseAll)
    for j = 1:length(ipulses)
        [Array1, Array2] = meshgrid(decompoPulseAll{i}, ipulses{j});
        diff_values = Array1 - Array2;
        valid_elements = diff_values <= 15*2 & diff_values >= 0;
        count = sum(valid_elements(:));
        r = count/(length(decompoPulseAll{i})+length(ipulses{j})-count);
        if r > 1
            r = 1;
        end
        spike_ROA_matrix(i, j) = r;
        matchresult_time_raw(end+1,:) = [i, j, r];
    end
end
matchresult_time = matchresult_time_raw;
for mu = 1:length(decompoPulseAll)
    tmpInd = find(matchresult_time(:,1) == mu);
    if length(tmpInd) > 1
        [~, tmpInd2] = max(matchresult_time(tmpInd,3));
        deleteInd = setdiff(1:length(tmpInd), tmpInd2);
        matchresult_time(tmpInd(deleteInd), :) = [];
    end
end
% matchresult_time1 = matchresult_time;
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

plotDecomps(ipulses, [], 2000, 0, 0, []);

save(['./Results/datasets' datasets_num '_resultnew.mat'], 'decompoCoVAll', 'decompoPulseAll', 'decompoSourceFirstAll', 'decompoSourceAll', 'matchresult_time', 'matchresult_time_raw', 'spike_ROA_matrix');