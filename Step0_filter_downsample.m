fsampu=1000;
for M=1:3
    for L=1:3
        for T=3
            % 需要修改此处file和下面的file
            tvifile = ['E:\TVIData\25-12-02\25-12-02\TVIData_15000_S_wy_shen_M' num2str(M) '_level' num2str(L) '_trial' num2str(T) '_Single_25-12-02'];
            % intra_pt=eaf2pulse(eaffile,fsampu);%导入肌内肌电的放电时刻
            tmp = load([tvifile,'.mat']);
            tvi_raw{1}=cat(3,zeros(395,128,20),tmp.TVIData);%%为什么这里要cat 20-求TVI的时候去除了20帧的数据，所以补充回来
            clearvars tmp
            %% 对TVI滤波(不变)
            TVIData_filter = tvi_raw;
            lowFre = 0.5/(7.7*4)*2;
            % lowFre = 1/(7.7*4)*2;
            % highFre= 14/(7.7*4)*2;
            % [Be1,Ae1] = butter(4,[lowFre,highFre]);%轴向低通0.5Mz
            [Be1,Ae1] = butter(4,lowFre,"low");%轴向低通0.5Mz
            [Be2,Ae2] = butter(4,[5,100]/fsampu*2);%时向滤波5-100Hz
            % [Be2,Ae2] = butter(4,5/fsampu*2,"low");
            for gn = 1
                for i = 1:size(TVIData_filter{gn},3)
                    %                                             [gn,i]
                    tmp = TVIData_filter{gn}(:,:,i);
                    tmp = filtfilt(Be1,Ae1,tmp);
                    TVIData_filter{gn}(:,:,i) = tmp;
                end
                for r = 1:size(TVIData_filter{gn},1)
                    for c = 1:size(TVIData_filter{gn},2)
                        %                         [gn,r,c]
                        tmp = TVIData_filter{gn}(r,c,:);
                        tmp = reshape(tmp,1,size(TVIData_filter{gn},3));
                        tmp = filtfilt(Be2,Ae2,tmp);
                        TVIData_filter{gn}(r,c,:) = tmp;
                    end
                end
            end
            frames_TVI=TVIData_filter{1};
            clearvars TVIData_filter;
            % frames_TVI=tmp{1}(:,:,1:end);
            clearvars tmp tvi_raw
            %%
            % 原数据尺寸
            [n1, n2, n3] = size(frames_TVI);

            % 计算需要补多少行
            pad_len = ceil(n1/3)*3 - n1;

            % 用 NaN 补到 3 的整数倍
            data_pad = cat(1, frames_TVI, NaN(pad_len, n2, n3));

            % 重排维度，把第一维按 3 分成一组
            data_group = reshape(data_pad, 3, [], n2, n3);

            % 求平均，忽略 NaN
            TVIData_filter{1}= squeeze(mean(data_group, 1, 'omitnan'));

            savename=['Z:\data\25-12-02\TVI_axial_downsample\downsampleTVIData_15000_S_wy_shen_M' num2str(M) '_level' num2str(L) '_trial' num2str(T) '_Single_25-12-02'];
            save(savename,"TVIData_filter")
            clearvars TVIData_filter data_group data_pad frames_TVI tmp
        end
    end
end