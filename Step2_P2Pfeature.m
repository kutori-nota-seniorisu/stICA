% 峰峰值 & Map & 聚类平均曲线（未进行refine）

clearvars -except twitchFrames_all varFrames_all

% 参数配置
removeNonBody = 0;    % 是否屏蔽非人体组织区域
ilim = 132; jlim = 80;    % 若屏蔽，组织区域的边界
threshold = 0.7;         % 二值化阈值比例默认0.7
silthreshold = 0.8;      % 得到的区域聚类的KMeans 最小Silhouette阈值，默认0.8

savefolder='Z:\Result\25-07-04前臂屈肌肌电超声\feature\median';
% savefolder='D:\ICdata';
motionnum = 4;
savename='p2p';

datalength=100;%一共多少帧的数据（还是默认1-100帧）


for M=1:motionnum

% load(['Z:\Result\25-07-16前臂屈肌肌电超声\USSTA_5-35Hz_M' num2str(M) '_0s_15s_-49_to_100_Frame_1000Hz.mat'])
% load(['Z:\Result\25-07-17前臂屈肌肌电超声\USSTA_5-35Hz_M' num2str(M) '_30ms_15s_1_to_100_Frame_1000Hz_median_CoV=0.7.mat'])
load(['Z:\Result\25-07-04前臂屈肌肌电超声\USSTA_5-35Hz_M' num2str(M) '_0s_15s_-49_to_100_Frame_2000Hz.mat'])
[levelnum, trialnum] = size(twitchFrames_all{M});

if removeNonBody==1
saveFile = [savefolder '\' savename '_M' num2str(M) '_ilim' num2str(ilim) '_jlim' num2str(jlim) '_1_100.mat'];
else
  saveFile = [savefolder '\' savename '_M' num2str(M) '_1_100.mat'];
end

% 1. 峰峰值提取 & 二值化 Map 
% for M = motionnum
    for L = 1:levelnum
        for T = 1:trialnum
            dataPack = twitchFrames_all{M}{L, T};
            if isempty(dataPack), continue; end
            munum = size(dataPack, 2);
            map =  {};  % 预分配

            for mu = 1:munum
                traces = dataPack(:, mu);            % 单/多探头数据
                pnum = numel(traces);
                tfpp = cell(pnum,1);

                for probe = 1:pnum
                    tmp3D = traces{probe};           % [ii, jj, time]
                    [ii, jj, tt] = size(tmp3D);
                    % 预分配并填充 tf_array
                    tf_array = cell(ii, jj);
                    for i = 1:ii
                        for j = 1:jj
                            vec = squeeze(tmp3D(i,j,:));
                            if removeNonBody && j>=jlim
                                vec = zeros(size(vec));
                            end
                             if removeNonBody && i>=ilim
                                vec = zeros(size(vec));
                            end
                            tf_array{i,j} = vec;
                        end
                    end
                    % 计算峰峰值
                    ppmap = zeros(ii, jj);
                    for i = 1:ii
                        for j = 1:jj
                            ppmap(i,j) = peak2peak(tf_array{i,j}(end-datalength+1:end));%%%%%%%%%%%%%%%%%
                        end
                    end
                    tfpp{probe} = ppmap;
                end
                property{M}{L,T}.tfpp{:,mu} = tfpp;
                % 二值化 Map
                for probe = 1:pnum
                    ppmap = tfpp{probe};
                    thresh = threshold * max(ppmap, [], 'all');
                    binMap = ppmap;
                    binMap(ppmap < thresh) = 0;
                    map{probe, mu} = binMap;
                end
            end
            map_all{M}{L,T} = map;
        end
    end
% end

% 2. 基于 Map 聚类 & 计算 Mean Twitch
curvesofmap_all = cell(size(map_all));
meantwitch_all  = cell(size(map_all));
% for M = motionnum
    for L = 1:levelnum
        for T = 1:trialnum
            map = map_all{M}{L,T};
            traces = twitchFrames_all{M}{L,T};
            if isempty(map), continue; end
            munum = size(map,2);

            for mu = 1:munum
                for probe = 1:size(map,1)
                    areaMask = map{probe,mu} > 0;
                    curves = get_single_sta(areaMask, traces{probe,mu});
                    nC = numel(curves);
                    if nC < 1, continue; end
                    % 用 vertcat 加速拼接
                    dataMat = vertcat(curves.Data);  % nC×T
                    Kmax = min(10, nC);
                    sils = -inf(Kmax,1);
                    idxAll = cell(Kmax,1);
                    for k = 2:Kmax
                        [idx, ~] = kmeans(dataMat, k, 'Replicates', 10);
                        sils(k) = mean(silhouette(dataMat, idx));
                        idxAll{k} = idx;
                    end
                    [bestSil, bK] = max(sils);
                    if bestSil < silthreshold
                        idx = ones(nC,1);
                        bK = 1;
                    else
                        idx = idxAll{bK};
                    end
                    % 存储分组信息与 Mean Twitch
                    for c = 1:nC
                        curves(c).group = idx(c);
                    end
                    curvesofmap_all{M}{L,T}{probe,mu} = curves;
                    for g = 1:bK
                        meantwitch_all{M}{L,T}{probe,mu}{g} = mean(dataMat(idx==g,:),1);
                    end
                end
            end
        end
    end
% end
%%
% saveFile = ['Z:\Result\25-07-04\p2p_M1_1-14MHz.mat'];
save(saveFile,"property",'meantwitch_all',"curvesofmap_all","map_all");
end
