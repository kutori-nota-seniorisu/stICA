clear; clc; close all;
addpath('./Func');
for usePSO = 0%1:-1:0
    for setset = 9%1:10
        for weight = 0.9%0.1:0.4:0.9
            for not_Xt = 1%0:1

                %% 参数设置
                % 潜在源的数目/主成分的个数
                numCompo = 12;
                % 数据集编号
                datasets_num = num2str(setset);
                % datasets_num = '2';

                % 是否从本地导入TVIdata
                flag_loadTVIdata = 0;

                % 分解方法
                % flag_decomp = 'SVD';
                % flag_decomp = 'NMF';
                flag_decomp = 'load_UVdata';

                % 模型选择
                model_tail = 'NC_NMFm';
                load_NMF_data_name = ['UV_compo' num2str(numCompo) '_' model_tail];

                flag_PSO = usePSO; % 针对优化算法，1表示使用"粒子群优化算法"，0表示使用"共轭梯度下降法"。

                flag_ICA = 1; %~ 是否使用ICA算法进行独立成分分析，迭代过程
                flag_plotXsXt = 1; %~ 是否绘制源图和源信号
                flag_plotUV = 0; %~ 是否绘制SVD/NMF分解后的信号
                flag_movie = 0; %~ 是否将混合后的TVI图存储成movie

                space_datasets = ['image_mat' datasets_num];% 空间分量数据集所在的文件夹
                time_datasets = ['TimeCompoDatasets' datasets_num];% 时间分量数据集所在的文件夹

                %% 保存迭代记录
                if ~exist('./Log', 'dir')
                    mkdir('./Log');
                end
                if not_Xt
                    if usePSO
                        diary(['./Log/log_set' datasets_num '_NMFm' num2str(weight) '_PSO_' char(datetime('now','Format','yyyy_MM_dd')) '.txt']);
                    else
                        diary(['./Log/log_set' datasets_num '_NMFm' num2str(weight) '_' char(datetime('now','Format','yyyy_MM_dd')) '.txt']);
                    end
                else
                    if usePSO
                        diary(['./Log/log_set' datasets_num '_NMFmXt' num2str(weight) '_PSO_' char(datetime('now','Format','yyyy_MM_dd')) '.txt']);
                    else
                        diary(['./Log/log_set' datasets_num '_NMFmXt' num2str(weight) '_' char(datetime('now','Format','yyyy_MM_dd')) '.txt']);
                    end
                end
                disp(['程序开始于:' char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'))]);

                %% 导入空间源成分
                folderPath = ['F:/EEEMG/stICA/Data/simulation/figure/' space_datasets '/'];
                fileFormat = '*.mat';
                % 使用dir函数获取文件夹中符合文件格式的文件信息
                fileList = dir(fullfile(folderPath, fileFormat));

                matFiles = {};
                for i = 1:numel(fileList)
                    if endsWith(fileList(i).name, '.mat')
                        matFiles{end+1} = fullfile(folderPath, fileList(i).name);
                    end
                end

                image_all = [];
                for i = 1:numel(matFiles)
                    load(matFiles{i});
                    image_all(:,:,i) = image_data;
                end

                %~ 高斯噪声图像
                imageGS1 = randn(400,128); image_all(:,:,end+1) = imageGS1*3;
                imageGS2 = randn(400,128); image_all(:,:,end+1) = imageGS2*3;

                Xs = reshape(image_all,size(image_all,1)*size(image_all,2),size(image_all,3));

                %% 导入时间源成分
                waveGS1 = randn(4000,1);
                waveGS2 = randn(4000,1);
                for i = 1:10
                    load(['./Data/simulation/MU_time_response/' time_datasets '/Time_component' num2str(i) '.mat']);
                    Xt_temp(:,i) = MU_noisy_conv(1:4000);
                    clear MU_noisy_conv;
                end
                Xt = [Xt_temp waveGS1 waveGS2];

                %% 下面是stICA算法的各部分实现
                params.dataSize = [400 128 numCompo];
                Ws0 = eye(numCompo)+randn(numCompo)*0.9; %~ Ws作为变量初始化

                %% 矩阵分解
                tic;
                disp('开始矩阵分解');

                if flag_loadTVIdata
                    disp('从本地导入TVI数据');
                    load(['./Data/simulation/datasets' datasets_num '/TVIdata_compo' num2str(numCompo) '.mat']);
                end

                if strcmp(flag_decomp,'SVD')
                    disp('采用了SVD分解')
                    if ~flag_loadTVIdata
                        TVIdata = Xs*Xt';
                    end

                    [U,D,V] = svd(TVIdata,0); %~ svd奇异值分解，0表示精简分解，减少计算量
                    Lambda = D(1:numCompo, 1:numCompo); %~ 缩放矩阵（大写的Lambda）
                    U_hat = U(:,1:numCompo)*sqrt(Lambda);%~ 取svd分解后的空间分量U的前k个主成分作为U_hat（降维）。疑问：*sqrt(Lambda)的作用？
                    V_hat = V(:,1:numCompo)*sqrt(Lambda);

                elseif strcmp(flag_decomp,'NMF')
                    disp('采用了NMF分解')
                    if ~flag_loadTVIdata
                        TVI_tmp = Xs*Xt';
                        % TVI_tmp_max = max(max(TVI_tmp));
                        % TVI_tmp_min = min(min(TVI_tmp));
                        % TVIdata = (TVI_tmp - TVI_tmp_min)/(TVI_tmp_max - TVI_tmp_min);
                        %~ 归一化
                        TVIdata = normalize(TVI_tmp,'range',[0 2]); %~ 按列，空间分量效果好，时间分量效果差
                        % TVIdata = normalize(TVI_tmp,2,'range',[0 2]);%~ 按行，空间分量效果差，时间分量效果好一些。进行缩放，使其范围在 [0,1] 区间内
                    end

                    %~ NMF分解，并行计算
                    if isempty(gcp('nocreate'))
                        parpool;
                    end
                    opt_nmf = statset('UseParallel',true);
                    [U_hat,V_tmp] = nnmf(TVIdata,numCompo,'options',opt_nmf);%~ NMF分解
                    V_hat = V_tmp';
                    D = eye(numCompo);

                elseif strcmp(flag_decomp,'load_UVdata')
                    disp('导入UV数据')
                    load(['./Data/simulation/datasets' datasets_num '/' load_NMF_data_name '.mat'])
                    clearvars S T

                else
                    disp('选择一种矩阵分解方法');
                end
                toc;
                disp(['矩阵分解用时：' num2str(toc)]);

                %% 参数设置
                % 潜在独立源成分
                params.k = numCompo;
                % 是否对时间成分进行pdf拟合，0表示不拟合
                % params.time = 0;
                % 是否使用带偏度的stICA算法，1表示skew-stICA，0表示stICA(无偏)
                params.skew_stICA = 1;
                % 带偏度的stICA前提下，1表示多概率密度分布函数，0表示单一概率密度分布函数
                params.MultiDistribution = 1;
                % 1表示使用V_hat拟合概率密度函数，0表示使用Xt插值拟合概率密度函数
                params.NMFt = not_Xt;
                % 多概率密度分布前提下，1表示使用NMF之后的空间分量作为标准拟合得到分布函数的参数，0表示使用最原始的空间分量作为标准拟合得到参数
                params.NMF = 1;
                % 概率密度分布函数中峰值对应的横坐标超过params.noise的认为是噪声
                params.noise = 5;
                % 时空ICA的权重参数,alpha=1时表示空间ICA，alpha=0时表示时间ICA
                params.alpha = weight;

                %% 自适应拟合概率密度分布函数，确定参数p、a、b、u
                tic;
                disp('开始概率密度函数参数估计')
                if params.MultiDistribution
                    if params.NMF
                        disp('多概率密度分布拟合 NMF')
                        pts = -10:0.05:50;
                        for i=1:size(U_hat,2)
                            zero_count = find(U_hat(:,i)==0);
                            U_temp = U_hat(:,i);
                            U_temp(zero_count) = [];
                            if ~isempty(U_temp)
                                [fs_temp(:,i), xs(:,i)] = ksdensity(U_temp,pts);
                                if ~isempty(find(fs_temp(:,i) > 1.5))
                                    fs(:,i) = normalize(fs_temp(:,i),'range');
                                else
                                    fs(:,i) = fs_temp(:,i);
                                end
                                [f_max(i), idx_s(i)] = max(fs(:,i));
                            else
                                fs(:,i) = 0;
                                xs(:,i) = 0;
                                f_max(i) = 0;
                                idx_s(i) = 1;
                            end
                        end

                        for i = 1:size(U_hat,2)
                            syms b p
                            %~ 确定参数a
                            a = 2;
                            %~ 求解u
                            u = floor(xs(idx_s(i),i)/4)*4;
                            if u == xs(idx_s(i),i)
                                u = u-1;%~ 防止u的值和最大值的横坐标相等导致a=0，即无偏概率分布函数；
                            end
                            if xs(idx_s(i),i) < params.noise
                                x_eqn = a/(sqrt(b^2-a^2)) + u == xs(idx_s(i),i);
                                %~ 求解b
                                b_temp = double(solve(x_eqn));
                                b = abs(b_temp(1)); %~ 取其中一个解的正数即可
                                %~ 求解p
                                x_eqn_value = a/(sqrt(b^2-a^2)) + u;%~ 概率密度分布函数的峰值对应的横坐标的值
                                p_eqn = p*exp(a*(x_eqn_value-u) - b*sqrt((x_eqn_value-u).^2+1)) == f_max(i);
                                p = double(solve(p_eqn));
                            else
                                a = 1;
                                b = 2;
                                u = 0;
                                p = 1;
                            end

                            %~ 参数整合
                            params.p(i) = p;
                            params.a(i) = a;
                            params.b(i) = b;
                            params.u(i) = u;
                        end
                    else
                        Xs_tmp = reshape(image_all(:,:,1:end-1),size(image_all,1)*size(image_all,2),size(image_all,3)-1);
                        Xs_tmp_max = max(max(Xs_tmp));
                        Xs_tmp_min = min(min(Xs_tmp));
                        Xsdata = (Xs_tmp - Xs_tmp_min)/(Xs_tmp_max - Xs_tmp_min);
                        image_sum = reshape(Xsdata,size(image_all,1),size(image_all,2),size(image_all,3)-1);
                        image_sum(:,:,size(image_all,3)) = image_gs*3;
                        Xs_new = reshape(image_sum,size(image_sum,1)*size(image_sum,2),size(image_sum,3));
                        Xs = normalize(Xs_new,'range',[0 2]);

                        pts = -10:0.05:50;
                        for i=1:size(Xs,2)
                            zero_count = find(Xs(:,i)==0);
                            U_temp = Xs(:,i);
                            U_temp(zero_count) = [];
                            if ~isempty(U_temp)
                                [fs_temp(:,i), xs(:,i)] = ksdensity(U_temp,pts);
                                if ~isempty(find(fs_temp(:,i) > 1.5))
                                    fs(:,i) = normalize(fs_temp(:,i),'range');
                                else
                                    fs(:,i) = fs_temp(:,i);
                                end
                                [f_max(i), idx_s(i)] = max(fs(:,i));
                            else
                                fs(:,i) = 0;
                                xs(:,i) = 0;
                                f_max(i) = 0;
                                idx_s(i) = 1;
                            end
                        end

                        for i = 1:size(Xs,2)
                            syms b p
                            %~ 确定参数a
                            a = 15;
                            %~ 求解u
                            u = floor(xs(idx_s(i),i)/3)*3;
                            if u == xs(idx_s(i),i)
                                u = u-1;%~ 防止u的值和最大值的横坐标相等导致a=0，即无偏概率分布函数；
                            end
                            if xs(idx_s(i),i) < params.noise
                                x_eqn = a/(sqrt(b^2-a^2)) + u == xs(idx_s(i),i);
                                %~ 求解b
                                b_temp = double(solve(x_eqn));
                                b = abs(b_temp(1)); %~ 取其中一个解的正数即可
                                %~ 求解p
                                x_eqn_value = a/(sqrt(b^2-a^2)) + u;%~ 概率密度分布函数的峰值对应的横坐标的值
                                p_eqn = p*exp(a*(x_eqn_value-u) - b*sqrt((x_eqn_value-u).^2+1)) == f_max(i);
                                p = double(solve(p_eqn));
                            else
                                a = 1.5;
                                b = 2.5;
                                u = 0;
                                p = 1;
                            end

                            %~ 参数整合
                            params.p(i) = p;
                            params.a(i) = a;
                            params.b(i) = b;
                            params.u(i) = u;
                        end
                    end
                end
                toc;
                disp(['概率密度函数参数估计用时：' num2str(toc)]);

                %% 时间成分拟合概率密度函数
                if params.NMFt
                    % 存储插值函数结构体以及一阶导数结构体

                else
                    % 存储插值函数结构体以及一阶导数结构体
                    pts = -3 : 0.005 : 3;
                    for i = 1:size(Xt,2)
                        [f, x] = ksdensity(Xt(:,i), pts);
                        pp = spline(x, f);
                        pp_d = fnder(pp, 1);
                        params.pp(i) = pp;
                        params.pp_d(i) = pp_d;
                    end
                end

                %% ICA算法的共轭梯度下降迭代
                if flag_ICA
                    tic;
                    Ws = Ws0;
                    fun = @(Ws) f_function(Ws,U_hat,V_hat,params);

                    if flag_PSO
                        % 参数需要优化调整，暂时不用。
                        % 并行计算能够大幅提高优化速度。
                        %~ 粒子群优化算法
                        % options = optimoptions('particleswarm', ...
                        %     'ObjectiveLimit',1e-3, ...
                        %     'Display', 'iter',...
                        %     'PlotFcn', @pswplotbestf);
                        if isempty(gcp('nocreate'))
                            parpool;
                        end

                        disp('使用粒子群优化算法');
                        options = optimoptions('particleswarm', ...'FunctionTolerance',1e-3,...'SwarmSize',60,...'ObjectiveLimit',1e-3, ...
                            'Display', 'iter',...
                            'MaxStallIterations', 100,...
                            'OutputFcn', @psoOutputFcn,...
                            'UseParallel', true,...
                            'PlotFcn', @pswplotbestf);
                        lb = repmat(-50,numCompo);
                        ub = repmat(50,numCompo);
                        [Ws,fval,exitflag,output] = particleswarm(fun,numCompo.^2,lb,ub,options); %~ 粒子群优化算法，全局
                        disp(output.message);

                        %         [Ws,fval] = particleswarm(fun,num_IC.^2,[],[],options); %~ 粒子群优化算法，全局
                        %     [Ws,fval] = particleswarm(fun,num_IC.^2);

                        %     %~ 遗传算法
                        %     options = optimoptions('ga','ConstraintTolerance',1e-3);
                        %     [Ws,fval] = ga(fun,num_IC.^2,[],[],[],[],[],[],[],options);

                    else
                        % disp('使用粒子群优化算法与共轭梯度法');
                        % options_ps = optimoptions('particleswarm', ...
                        %     'FunctionTolerance',1e-3,...
                        %     'Display', 'iter',...
                        %     'PlotFcn', @pswplotbestf);
                        % lb = repmat(-50,num_IC);
                        % ub = repmat(50,num_IC);
                        % [Ws_ps,fval_ps,exitflag_ps,output_ps] = particleswarm(fun,num_IC.^2,lb,ub,options_ps);
                        %
                        % options_fmin = optimoptions('fminunc',...
                        %     'MaxIterations',28800, ...
                        %     'MaxFunctionEvaluations',28800000, ...
                        %     'OptimalityTolerance',1e-6, ...
                        %     'StepTolerance',1e-8,...
                        %     'Display','iter', ...
                        %     'PlotFcn','optimplotfval');
                        % [Ws,fval_opt,exitflag_opt,output_opt] = fminunc(fun,Ws_ps,options_fmin);

                        % 共轭梯度法
                        % 并行计算似乎无法显著提升优化速度。
                        % options = optimset('Display','iter','PlotFcns',@optimplotfval);
                        disp('使用共轭梯度法');
                        options = optimoptions('fminunc',...
                            'MaxIterations',28800, ...
                            'MaxFunctionEvaluations',500000, ...'OptimalityTolerance',1e-6, ...
                            'StepTolerance',eps,...'SpecifyObjectiveGradient',true,...
                            'OutputFcn',@myOutputFcn,...
                            'Display','iter', ...'PlotFcn','optimplotfirstorderopt'
                            'PlotFcn','optimplotfval');
                        % valid = checkGradients(fun, Ws, 'Display', 'on');
                        [Ws,fval,exitflag,output] = fminunc(fun,Ws,options);
                        disp(output.message);

                        if exitflag == 5
                            disp('进行第二次迭代，采用中心差分');
                            options2 = optimoptions('fminunc',...
                                'MaxIterations',28800, ...
                                'MaxFunctionEvaluations',500000, ...'OptimalityTolerance',1e-6, ...
                                'StepTolerance',eps,...
                                'FiniteDifferenceType', 'central',...
                                'FiniteDifferenceStepSize', 1e-6,...
                                'OutputFcn',@myOutputFcn,...
                                'Display','iter', ...'PlotFcn','optimplotfirstorderopt'
                                'PlotFcn','optimplotfval');
                            [Ws,fval2,exitflag2,output2] = fminunc(fun,Ws,options2);
                        end
                    end
                    toc;
                    disp(['优化算法用时：' num2str(toc)]);

                    Ws = reshape(Ws,numCompo,numCompo);
                    Wt = inv(Ws');
                    S = U_hat*Ws;
                    T = V_hat*Wt;
                end

                %% 绘图分析
                if flag_plotXsXt
                    %~ 绘制源图和源信号
                    figure;
                    for i=1:10
                        subplot(2,10,i)
                        imagesc(reshape(Xs(:,i),[params.dataSize(1) params.dataSize(2)]));
                        title(['Space #' num2str(i)]);
                        subplot(2,10,i+10)
                        plot(-1*Xt(:,i));
                        title(['Temporal #' num2str(i)]);
                    end
                    set(gcf,'unit','normalized','position',[0.1,0.6,0.8,0.32]);
                    sgtitle('源成分图像')
                end

                if flag_plotUV
                    %~ 绘制SVD/NMF分解后的信号
                    Xs_plot = reshape(U_hat,400,128,numCompo);
                    figure;
                    for i = 1:numCompo
                        subplot(2,6,i)
                        imagesc(Xs_plot(:,:,i));
                    end
                    sgtitle('NMF分解后的空间分量')

                    Xt_plot = V_hat;
                    figure
                    plot(Xt_plot(:,1),'b','linewidth',2);
                    hold on
                    plot(Xt_plot(:,2),'k','linewidth',2);
                    hold on
                    plot(Xt_plot(:,3),'g','linewidth',2);
                    hold on
                    plot(Xt_plot(:,4),'r','linewidth',2);
                    if numCompo == 5
                        hold on
                        plot(Xt_plot(:,5),'y','linewidth',1);
                    end
                    sgtitle('分解后的时间分量')
                    legend('正弦波1','正弦波2','锯齿波','三角波','高斯噪声')
                end

                if flag_ICA
                    figure;
                    for i=1:numCompo
                        subplot(2,numCompo,i)
                        imagesc(reshape(S(:,i),[params.dataSize(1) params.dataSize(2)]));
                        title(['Space #' num2str(i)]);
                        % colorbar
                        subplot(2,numCompo,i+numCompo)
                        plot(-1*T(:,i));
                        title(['Temporal #' num2str(i)]);
                    end
                    set(gcf,'unit','normalized','position',[0.1,0.6,0.8,0.32]);
                    sgtitle('处理前图像')
                end


                %% 混叠后的TVIdata绘图
                if flag_movie
                    TVIData = reshape(TVIdata,400,128,[]);

                    mkdir(['./Movie/' date]);%~ 创建视频保存文件夹
                    Movie_file = dir(fullfile(['./Movie/' date],'*.mp4'));
                    if isempty(Movie_file)
                        MoviePath = ['./Movie/' date '/TVI_Movie_1.mp4']; %~ 视频文件的名称
                    else
                        MoviePath = ['./Movie/' date '/TVI_Movie_' num2str(size(Movie_file,1)+1) '.mp4'];
                    end

                    profile = 'MPEG-4'; %~ 视频文件格式
                    writerObj = VideoWriter(MoviePath,profile); %~ 创建视频文件
                    open(writerObj); %~ 打开该视频文件

                    tic
                    figure;
                    for n = 1:size(TVIData,3)
                        TVI_img = TVIData(:,:,n);
                        TVI_img(isnan(TVI_img)) = 0;
                        imagesc(TVI_img)
                        colorbar
                        caxis([-1 1]*max(max(max(TVIData))));
                        title('Tissue Velocity map')
                        ylabel('[cm]')
                        axis equal ij tight
                        set(gca,'XColor','none','box','off')
                        frame = getframe; %~ 获取图像帧
                        writeVideo(writerObj,frame); %~ 将帧写入视频文件中
                    end
                    toc
                    close(writerObj); %~ 关闭视频文件句柄
                end

                %% 关闭日志文件
                disp(['程序结束于:' char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'))]);
                diary off;

                %% 保存变量，在运行完main_simulation_multiChannel_iEMG后执行
                savepath = ['./Data/simulation/datasets' datasets_num ];
                if ~exist(savepath, 'dir')
                    mkdir(savepath);
                end
                % 保存分解结果
                if params.NMFt
                    if usePSO
                        % save([savepath '/' load_NMF_data_name '.mat'], 'S','T','U_hat','V_hat','Xs','Xt','Ws0','Ws','Wt','params');
                        save([savepath '/UV_compo12_NC_NMFm' num2str(params.alpha) '_PSO.mat'], 'S','T','U_hat','V_hat','Xs','Xt','Ws0','Ws','Wt','params');
                    else
                        save([savepath '/UV_compo12_NC_NMFm' num2str(params.alpha) '.mat'], 'S','T','U_hat','V_hat','Xs','Xt','Ws0','Ws','Wt','params');
                    end
                else
                    if usePSO
                        save([savepath '/UV_compo12_NC_NMFmXt' num2str(params.alpha) '_PSO.mat'], 'S','T','U_hat','V_hat','Xs','Xt','Ws0','Ws','Wt','params');
                    else
                        save([savepath '/UV_compo12_NC_NMFmXt' num2str(params.alpha) '.mat'], 'S','T','U_hat','V_hat','Xs','Xt','Ws0','Ws','Wt','params');
                    end
                end
                disp('数据保存完成');
                close all; clc;
            end
        end
    end
end