% clear all
framenum=100;
threshold=0.7;
motionnum=1;
levelnum=1;
trialnum=10;
silthreshold=0.8;
% saveFile = ['Z:\Result\25-07-04\feature\emma_cut1-10.mat'];
saveFile = ['D:\ICdata\emma.mat'];

for M = 1:motionnum
    % clearvars twitchFrames_all varFrames_all
    % datafile=['Z:\Result\25-07-04\USSTA\USSTA_M' num2str(M) '_0s_15s_-49_to_100_Frame_2000Hz.mat'];
    % load(datafile);
    % [levelnum, trialnum] = size(twitchFrames_all{M});
    for L = 1:levelnum
        for T = 1:trialnum
            twitch_trial = twitchFrames_all{M}{L, T};
            var_trial=varFrames_all{M}{L,T};
            munum = size( twitch_trial, 2);
            map =  {};  % 预分配
            for mu = 1:munum
                traces =  twitch_trial(:, mu);            % 单/多探头数据
                pnum = numel(traces);
                emma = cell(pnum,1);
                for probe = 1:pnum
                    frames_sta=twitch_trial{probe,mu};
                    var_sta=var_trial{probe,mu};
                    directions = sum(frames_sta(:,:,1:framenum/2)./var_sta(:,:,1:framenum/2),3)-sum(frames_sta(:,:,framenum/2+1:framenum)./var_sta(:,:,framenum/2+1:framenum),3);
                    directions_sign = sign(directions);
                    MUintensity = sum(frames_sta.^2./var_sta,3).*-directions_sign;
                    % MUintensity(1:10,:)=0;
                    % MUintensity_resample = downsample(MUintensity,3);%对轴向降采样3倍，使每个像素约为0.3*0.3mm
                    emma{probe}=MUintensity;
                    thresh = threshold * max(abs(MUintensity(:)));
                    binMap=MUintensity;
                    binMap(abs(MUintensity)<thresh)=0;
                    map{probe,mu}=binMap;
                    
                end
                property{M}{L,T}.emma{:,mu} = emma;
            end
            map_all{M}{L,T}=map;
        end
    end
end

curvesofmap_all = cell(size(map_all));
meantwitch_all  = cell(size(map_all));
for M = 1:motionnum
    % clearvars twitchFrames_all varFrames_all
    % datafile=['Z:\Result\25-07-04\USSTA\USSTA_M' num2str(M) '_0s_15s_-49_to_100_Frame_2000Hz.mat'];
    % load(datafile);
    for L = 1:levelnum
        for T = 1:trialnum
            
            map = map_all{M}{L,T};
            traces = twitchFrames_all{M}{L,T};
            if isempty(map), continue; end
            munum = size(map,2);

            for mu = 1:munum
                for probe = 1:size(map,1)
                    areaMask = abs(map{probe,mu}) > 0;
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
end

%%

save(saveFile,"property",'meantwitch_all',"curvesofmap_all","map_all");
