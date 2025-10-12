% 对时间成分与空间成分进行匹配，并保留同时在时域与空域上匹配的结果
%% 时间成分匹配
load(['Data/simulation/MU_time_response/TimeCompoDatasets' datasets_num '/ipulses.mat']);

matchresult_time_raw = [];
% 计算RoA
for i = 1:length(decompo_pulses)
    for j = 1:length(ipulses)
        [Array1, Array2] = meshgrid(decompo_pulses{i}, ipulses{j});
        diff_values = Array1 - Array2;
        valid_elements = diff_values <= (exFactor+10)*2 & diff_values >= 0;
        count = sum(valid_elements(:));
        r = count/(length(decompo_pulses{i})+length(ipulses{j})-count);
        if r > 1
            r = 1;
        end
        spike_ROA_matrix(i, j) = r;
        matchresult_time_raw(end+1,:) = [i, j, r];
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
matchresult_time_raw = array2table(matchresult_time_raw, 'VariableNames', {'decomp', 'ref', 'time'});
matchresult_time = array2table(matchresult_time, 'VariableNames', {'decomp', 'ref', 'time'});

%% 空间成分匹配
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

%% 时空结果匹配
matchresult_final = innerjoin(matchresult_space, matchresult_time, 'Keys', {'decomp', 'ref'});

%% 空间图像二值化处理
for i = 1:size(S, 2)
    % thresh = 0.5 * max(abs(S(:, i)));
    tmpS = abs(S(:, i));
    thresh = 0.5 * (max(tmpS) - min(tmpS)) + min(tmpS);
    tmpInd = find(tmpS > thresh);

    % 没有转换成二维坐标
    actArea{i} = tmpInd;

    % 空间图像二值化
    tmpS(find(tmpS <= thresh)) = 0;
    tmpS(find(tmpS > thresh)) = 1;
    S_bin(:, i) = tmpS;
end

%% 绘制空间成分与源成分的对比图
for mi = 1:size(matchresult_space_raw, 1)
    decomp = matchresult_space_raw.decomp(mi);
    ref = matchresult_space_raw.ref(mi);

    figure;
    % 绘制空间成分
    subplot(1,5,1);
    imagesc(reshape(S(:,decomp),400,128)); title(['S' num2str(decomp)]);
    colorbar;
    % 绘制二值化空间成分
    subplot(1,5,2);
    imagesc(reshape(S_bin(:,decomp),400,128)); title(['S' num2str(decomp) ' (0-1)']);
    % 绘制空间源成分
    subplot(1,5,3);
    imagesc(reshape(Xs(:,ref),400,128)); title(['Xs' num2str(ref)]);
    colorbar;
    % 绘制空间成分的频数统计
    subplot(1,5,4);
    histogram(S(:, decomp)); title(['S' num2str(decomp)]);
    % 绘制空间源成分的频数统计
    subplot(1,5,5);
    histogram(Xs(:, ref)); title(['Xs' num2str(ref)]);
    set(gcf,'unit','normalized','position',[0.1,0.3,0.8,0.4]);
end

%% 绘制时空匹配成分图
for mi = 1:size(matchresult_final, 1)
    decomp = matchresult_final.decomp(mi);
    ref = matchresult_final.ref(mi);

    figure;
    % 空间匹配
    subplot(4,2,[1,3]);
    imagesc(reshape(S(:,decomp),400,128));
    title(['S' num2str(decomp) ',Corr = ' num2str(matchresult_final.space(mi))]);
    subplot(4,2,[2,4]);
    imagesc(reshape(Xs(:,ref),400,128));
    title(['Xs' num2str(ref)]);

    % 时间匹配
    subplot(4,2,[5,6]);
    plot(T_norm(:,decomp))
    title(['T' num2str(decomp)]);
    subplot(4,2,[7,8]);
    plot(Xt(:,ref))
    title(['Xt' num2str(ref)]);

    set(gcf,'unit','normalized','position',[0.35,0.1,0.3,0.7]);

    plotDecomps({decompo_pulses{decomp}, ipulses{ref}}, [], 2000, 0, 0, []);
    yticks(1:2);
    yticklabels({['T' num2str(decomp)], ['Xt' num2str(ref)]});
    title(['RoA=' num2str(matchresult_final.time(mi))]);
end


