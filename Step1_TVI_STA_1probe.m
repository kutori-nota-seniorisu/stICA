% 使用TVI和sEMG的放电串得到TVI-STA和VAR数据
% by KYM 25-7


%% 参数
clearvars -except CKCdecomp trackingdecomp

fsampemg=2048;%肌电采集频率
fsampu = 2000;%探头帧率
framenum = 100;%总帧数
ignoreN   = 30;     % 想忽略的帧数
starttime = ignoreN / fsampu;   % 原来是 0，现在改成 30 帧对应的时间
starttime_ms = round(starttime *fsampu);
endtime = 15;
mode = 2;%STA的方法:1-前后framenum帧 2-后framenum帧 3- -50到100帧
activationDir = 'Z:\Result\25-07-04前臂屈肌肌电超声';  % 可更改为任意目标路径
saveflag = 1; % 是否保存

% 导入CKC
decompfile=['Z:\Result\25-07-04\MUST\含反解\5-35Hz\CKCdecomp_M' num2str(M) '.mat'];
decompuse=CKCdecomp;

gn=1;%单探头
filenameTemplate = 'USSTA_5-35Hz_M%d_%dms_%ds_%d_to_%d_Frame_%dHz_median.mat';  % 参数化命名模板
[start_frame, end_frame] = getSTAFrameRange(framenum, mode);
for M=1:4
    % decompfile=['Z:\Result\25-07-04\MUST\含反解\5-35Hz\CKCdecomp_M' num2str(M) '.mat'];
    % decompuse=importdata(decompfile,"CKCdecomp");
    for L=1:4
        for T=1:2
            %% 导入数据
            decomp=decompuse{M}{L,T};
            MUpt=decomp.MUPulses;
            mucov=decomp.CoV;
            tvifile = ['D:\BaiduNetdiskDownload\25-07-04\downsampleTVIData_30000_S_wrl_M' num2str(M) '_level' num2str(L) '_trial' num2str(T) '_Single_25-07-04'];
            % tvifile = ['Z:\data\25-07-17单指EMG+US\TVI_axial_downsample\downsampleTVIData_15000_S_hzy_finger_M' num2str(M) '_level' num2str(L) '_trial' num2str(T) '_Single_25-07-17'];    
            load([tvifile,'.mat']);
            
            %% 对齐放电串
            emgdata=decomp.datafilt;
            if size(emgdata,1)==128
                celltype=7;
            else
                celltype=6;
            end
            emgcell=sig2cell(emgdata,celltype);
            muaplength=100;
            muaps_all=muapExtraction(emgcell,MUpt,muaplength,'STA');% 中心点在N/2+1;
            num_mu=length(MUpt);
            shift_ap =  zeros(1,num_mu);
            peakdetect=zeros(1,num_mu);
            for mu = 1:num_mu
                % for mu=3
                tmpArray = plotArrayPotential(muaps_all{mu},1,0);
                tmpArray_diff2 = cell(size(tmpArray));
                tmpArray_diff2 = cellfun(@(x) diff(diff(x)),tmpArray,'UniformOutput', false);
          
                [~,~,~,pos] = plotArrayPotential(tmpArray_diff2,1,0);
                tmptmp = tmpArray_diff2{pos(1),pos(2)};

                z = (tmptmp - mean(tmptmp)) / std(tmptmp);% Zscore一下
                tmpInd = find(abs(z) > 5);%找到大于5倍std的一次
                if isempty(tmpInd)%防止信噪比低，如果没有大于5倍std的就找大于3倍std的。
                    tmpInd=find(abs(z)>3);
                end
                if ~isempty(tmpInd)
                    shift_ap(mu) = tmpInd(1)-(muaplength/2+1);
                    peakdetect(mu)=1;
                end
            end
            pulses_new = cellfun(@(x, s) x + s, MUpt, num2cell(shift_ap), 'UniformOutput', false);

            
            %% 通过sEMG和肌内肌电得到MUAP（认为肌内肌电的时间是准的所以不用调整肌内肌电的时延），计算每一帧的每个点的intensity
            %初始化
            twitchFrames = cell(0);%sta后的曲线
            varFrames = cell(0);%sta后的曲线
            num_mu=size(pulses_new,2);%得到肌内的MU数量
            %对每个MU，对每个探头都做一次STA
            for mu=1:num_mu
                if mucov(mu)>0.7
                    continue
                end
                tmp_pt1=pulses_new{mu};
                tmp_pt=round(tmp_pt1*fsampu/fsampemg);
                [frames_sta,var_sta] = ultrasoundSTA(framenum,tmp_pt,TVIData_filter,gn,fsampu,starttime,endtime,mode);%计算一共M帧的STA和Var
                twitchFrames{gn,mu} = frames_sta;
                % varFrames{gn,mu}=downsample(var_sta,3);
                varFrames{gn,mu}=var_sta;
                disp([ 'L/T/#/MU:' num2str(L) '/' num2str(T) '/' num2str(mu) '/' num2str(num_mu)])
               
            end
            twitchFrames_all{M}{L,T}=twitchFrames;
            varFrames_all{M}{L,T}=varFrames;
            clearvars TVIData_filter
        end
    end
    %% 保存
    if saveflag
        filename = sprintf(filenameTemplate, M,starttime_ms,endtime, start_frame,end_frame, fsampu);
        activationFile = fullfile(activationDir, filename);
        % activationFile = ['Z:\Result\肱二头肌针刺双超声STA\USSTA_M1_0-15s_-50_100Frame_2000Hz'];
        save(activationFile,"twitchFrames_all","varFrames_all");
    end
    clearvars twitchFrames_all varFrames_all

end



function [frames_sta, var_sta] = ultrasoundSTA(framenum, tmp_pt, TVIData_filter, gn, fsampu, starttime, endtime, mode)
    % 1. 根据 mode 确定相对帧的偏移范围 (n_range)
    switch mode
        case 1 % 前后对称
            n_range = (-framenum/2 + 1) : (framenum/2);
        case 2 % 仅放电后
            n_range = 1 : framenum;
        case 3 % 前50后100 (非对称)
            pre_frames = framenum / 2; 
            post_frames = framenum;
            n_range = (-pre_frames + 1) : post_frames;
        otherwise
            error('Mode must be 1, 2, or 3');
    end

    % 2. 预分配输出内存 (获取图像尺寸)
    [imgH, imgW, ~] = size(TVIData_filter{gn});
    num_frames_out = length(n_range);
    
    frames_sta = zeros(imgH, imgW, num_frames_out);
    var_sta    = zeros(imgH, imgW, num_frames_out);
    
    % 预计算边界阈值，避免在循环中重复乘法
    idx_min = starttime * fsampu;
    idx_max = endtime * fsampu;

    % 3. 统一的主循环
    for k = 1:num_frames_out
        n = n_range(k); % 当前相对于放电时刻的偏移帧数
        
        % 计算绝对索引并进行边界过滤
        tmp_ind = tmp_pt + n;
        tmp_ind(tmp_ind <= idx_min | tmp_ind >= idx_max) = [];

        % 如果当前没有任何有效的触发点，保持为0并跳过
        if isempty(tmp_ind)
            continue;
        end

        % 提取当前延迟下的所有图像帧 [H x W x NumSpikes]
        tmp_tvi = TVIData_filter{gn}(:, :, tmp_ind);
        num_spike = length(tmp_ind);

        
        % 1. 计算中心趋势 (Median 或 Mean)
        stam = median(tmp_tvi, 3); 
        
        % 2. 计算离散度 (方差)
        % 优化：移除原本的 for j 循环，使用矩阵减法 (Implicit Expansion) 加速
        % 计算 (x - center)^2
        diff_sq = (tmp_tvi - stam).^2; 
        % 在第三维度求和并除以 N-1
        varm = sum(diff_sq, 3) / (num_spike - 1);
        
        % ---------------------------------------

        % 存入结果矩阵
        frames_sta(:, :, k) = stam;
        var_sta(:, :, k)    = varm;
    end
end



function [start_frame, end_frame] = getSTAFrameRange(framenum, mode)
% 返回相对于放电时刻的STA起止帧编号
% mode = 1：对称前后帧
% mode = 2：仅使用放电后的帧
% mode = 3：前50后100（固定比例）

    switch mode
        case 1 % 使用前后对称 framenum
            start_frame = -framenum/2 + 1;
            end_frame = framenum/2;

        case 2 % 使用放电后的 framenum 帧
            start_frame = 1;
            end_frame = framenum;

        case 3 % 放电前50帧，后100帧，共150帧
            start_frame = -framenum/2 + 1; % -49（如果framenum=100时前半段）
            end_frame = framenum;         % +100（后半段）

        otherwise
            error('模式 mode 必须为 1、2 或 3');
    end
end