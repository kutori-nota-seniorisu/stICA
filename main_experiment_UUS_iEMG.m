clear; clc; close all;
addpath('./Func')
%% 修改参数
for level = 1%:4
    for trial = 1%:2
        for probe = 1%:2
            % import_file_name = 'TVI_Data_4000_CPWC1_22-11-14_3Mode_trial1.mat';
            % import_file_name = 'TVIData_StandardTest';
            % import_file_name = 'TVI_Data_10000_CPWC1_22-11-14_3Mode_trial1_new';
            import_file_name = ['TVIData_S1_M1_level' num2str(level) '_trial' num2str(trial) '_Dual_24-06-21_' num2str(probe)];
            % import_file_name = 'TVI_Data_20000_CPWC1_23-01-06_trial3_new';
            % import_file_name = 'TVI_Data_8000_CPWC1_22-11-14_3Mode_trial1_new';

            file_date = '24-06-21/UUS-iEMG/';

            % 导入数据
            load(['./Data/experiment/' file_date import_file_name '.mat'])
            TVIdata = TVIData(:,:,3001:13000);
            params.dataSize = [size(TVIdata,1) size(TVIdata,2) size(TVIdata,3)];

            flag_segment = 0; % 是否分段处理
            datalength = 4000; %~ 每段数据长度

            % 矩阵分解方法
            % flag_decomp = 'SVD';
            flag_decomp = 'NMF';
            % flag_decomp = 'load_NMF_data';

            flag_ICA = 1; %~ 是否使用ICA算法进行独立成分分析，迭代过程

            params.skew_stICA = 1; %~ 是否使用带偏度的stICA算法，1表示skew-stICA，0表示stICA(无偏)
            params.MultiDistribution = 1; %~ 带偏度的stICA前提下，1表示1多概率密度分布函数”，0表示“单一概率密度分布函数”
            params.noise = 10; %~ 概率密度分布函数中峰值对应的横坐标超过params.noise的认为是噪声

            params.k = 12; %~ 潜在源的数目/主成分的个数
            params.alpha = 0.8; %~ 时空ICA的权重参数,alpha=1时表示空间ICA，alpha=0时表示时间ICA
            params.mode = 'st';

            %% 下面是stICA算法的各部分实现
            % for i = 1:floor(params.dataSize(3)/length)
            for seg_i = 1
                %% 矩阵分解
                tic;
                disp('开始矩阵分解')
                if strcmp(flag_decomp,'SVD')
                    if flag_segment
                        TVIdata = TVIData(:,:,(seg_i-1)*datalength+1:seg_i*datalength); %~ 截取某一段数据分析
                    end

                    % 数据预处理
                    TVI_tmp = stripmean(TVIdata,params.mode);%~ 时空组织速度图数据零均值化
                    TVI_data =reshape(TVI_tmp,params.dataSize(1)*params.dataSize(2),[]); %~ 三维数据重排为二维数据

                    Ws0 = eye(params.k)+randn(params.k)*0.9; %~ Ws作为变量初始化

                    %~ SVD分解
                    [U,D,V] = svd(TVI_data,0); %~ svd奇异值分解，0表示精简分解，减少计算量
                    Lambda = D(1:params.k, 1:params.k); %~ 缩放矩阵（大写的Lambda）
                    U_hat = U(:,1:params.k)*sqrt(Lambda);%~ 取svd分解后的空间分量U的前k个主成分作为U_hat（降维）。疑问：*sqrt(Lambda)的作用？
                    V_hat = V(:,1:params.k)*sqrt(Lambda);

                elseif strcmp(flag_decomp,'NMF')
                    if flag_segment
                        TVIdata = TVIData(:,:,(seg_i-1)*datalength+1:seg_i*datalength); %~ 截取某一段数据分析
                    end

                    % 数据预处理
                    TVI_tmp = stripmean(TVIdata,params.mode);%~ 时空组织速度图数据零均值化
                    TVI_data =reshape(TVI_tmp,params.dataSize(1)*params.dataSize(2),[]); %~ 三维数据重排为二维数据

                    Ws0 = eye(params.k)+randn(params.k)*0.9; %~ Ws作为变量初始化

                    % 归一化方式选择
                    % 按列归一化
                    TVI_norm = normalize(TVI_data, 'range', [0 2]);

                    % Min-Max归一化（Min-Max Normalization）
                    % TVI_data_max = max(max(TVI_data));
                    % TVI_data_min = min(min(TVI_data));
                    % TVI_norm = (TVI_data - TVI_data_min)/(TVI_data_max - TVI_data_min);

                    % NMF分解
                    % optNMF = statset('UseParallel', true);
                    [U_hat,V_tmp] = nnmf(TVI_norm, params.k);%, 'options', optNMF);
                    D = eye(params.k);
                    V_hat = V_tmp';

                elseif strcmp(flag_decomp,'load_NMF_data')
                    load('./Data/experiment/U_hat.mat');
                    load('./Data/experiment/V_hat.mat');
                    load('./Data/experiment/Ws0.mat');

                else
                    disp('选择一种矩阵分解方法');
                end
                toc;
                disp(['矩阵分解用时：' num2str(toc)]);

                %% 自适应拟合概率密度分布函数，确定参数p、a、b、u
                tic;
                disp('开始概率密度函数参数估计')
                if params.MultiDistribution
                    % 空间分布拟合
                    pts = -10:0.05:50;
                    for i=1:size(U_hat,2)
                        zero_count = find(U_hat(:,i)==0);
                        U_temp = U_hat(:,i);
                        U_temp(zero_count) = [];
                        if ~isempty(U_temp)
                            % ksdensity核心平滑密度估计
                            [fs_temp(:,i), xs(:,i)] = ksdensity(U_temp,pts);
                            % 如果存在大于1.5的密度值，就对估计的概率密度函数进行归一化。
                            if ~isempty(find(fs_temp(:,i) > 1.5))
                                fs(:,i) = normalize(fs_temp(:,i),'range');
                                disp('归一化')
                            else
                                fs(:,i) = fs_temp(:,i);
                            end
                            % 找到成分i的概率密度函数最大值以及其索引
                            [fs_max(i), idx_s(i)] = max(fs(:,i));
                        else
                            fs(:,i) = 0;
                            xs(:,i) = 0;
                            fs_max(i) = 0;
                            idx_s(i) = 1;
                        end
                    end

                    for i = 1:size(U_hat,2)
                        % 论文中的字母为d,q,μ,r
                        % p对应于d，a对应于q，u对应于μ，b对应于r
                        syms b p
                        %~ 确定参数a
                        a = 15;
                        %~ 求解u
                        % 参数u由PDF顶点的横坐标决定，此处将横坐标分割为长度为4的分段
                        u = floor(xs(idx_s(i),i)/4)*4;
                        if u == xs(idx_s(i),i)
                            u = u-1;%~ 防止u的值和最大值的横坐标相等导致a=0，即无偏概率分布函数；
                        end
                        if xs(idx_s(i),i) < params.noise
                            % 一个等式：a/(sqrt(b^2-a^2)) + u = x_ksdensity(idx(i),i)，其中b是未知数
                            x_eqn = a/(sqrt(b^2-a^2)) + u == xs(idx_s(i),i);
                            %~ 求解b
                            b_temp = double(solve(x_eqn));
                            b = abs(b_temp(1)); %~ 取其中一个解的正数即可
                            %~ 求解p
                            x_eqn_value = a/(sqrt(b^2-a^2)) + u;%~ 概率密度分布函数的峰值对应的横坐标的值
                            p_eqn = p*exp(a*(x_eqn_value-u) - b*sqrt((x_eqn_value-u).^2+1)) == fs_max(i);
                            p = double(solve(p_eqn));
                        else
                            % 噪声的CDF参数
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
                toc;
                disp(['概率密度函数参数估计用时：' num2str(toc)]);

                %% 迭代优化
                tic;
                disp('开始迭代优化')
                if flag_ICA
                    Ws = Ws0;
                    fun = @(Ws) f_function(Ws,U_hat,V_hat,params);

                    %~ 共轭梯度法
                    % options = optimset('MaxIter', 250,'MaxFunEvals',4e4,'Display','iter','PlotFcns',@optimplotfval);
                    disp('使用共轭梯度法');
                    options = optimoptions('fminunc',...
                        'MaxIterations',28800, ...
                        'MaxFunctionEvaluations',500000, ...'OptimalityTolerance',1e-6, ...
                        'Display','iter', ...'PlotFcn','optimplotfirstorderopt'
                        'PlotFcn','optimplotfval');
                    [Ws,fval,exitflag,output] = fminunc(fun,Ws,options);
                    disp(output.message);

                    %         %~ 粒子群优化算法
                    % %         options = optimoptions('particleswarm','ObjectiveLimit',1e-2);
                    % %         [Ws,fval] = particleswarm(fun,params.k.^2,[],[],options); %~ 粒子群优化算法，全局
                    %         options = optimoptions('particleswarm','ObjectiveLimit',1e-3, ...
                    %             'Display', 'iter',...
                    %             'PlotFcn', @pswplotbestf);
                    %         lb = repmat(-50,params.k.^2);
                    %         ub = repmat(50,params.k.^2);
                    %         [Ws,fval] = particleswarm(fun,params.k.^2,lb,ub,options); %~ 粒子群优化算法，全局

                    %~ 遗传算法
                    %         options = optimoptions('ga','ConstraintTolerance',1e-3);
                    %         [Ws,fval] = ga(fun,params.k.^2,[],[],[],[],[],[],[],options);

                    Ws = reshape(Ws,params.k,params.k);
                    Wt = inv(Ws');
                    S = U_hat*Ws;
                    T = V_hat*Wt;
                end
                toc;
                disp(['迭代优化用时：' num2str(toc)]);

                %% 数据后处理
                % T=(T(2:end,:)-mean(T(2:end,:),1))./std(T(2:end,:),1); %~ 幅值标准化，第一行的值过大？
                % T_norm=(T-mean(T))./std(T);

                %% 绘图
                figure;
                for i=1:params.k
                    subplot(2,params.k,2*(i-1)+1)
                    imagesc(reshape(S(:,i),[395 128]));
                    title(['Space #' num2str(i)]);
                    % colorbar
                    subplot(2,params.k,2*(i-1)+2)
                    plot(-1*T(:,i));
                    title(['Temporal #' num2str(i)]);
                end
                set(gcf,'unit','normalized','position',[0.1,0.6,0.8,0.32]);

                %% 保存变量
                savepath = ['./Data/experiment/' file_date 'S1M1L' num2str(level) 'T' num2str(trial) 'P' num2str(probe) '_compo' num2str(params.k) '_' flag_decomp 'm.mat'];
                save(savepath, 'S', 'T', 'Ws0', 'Ws', 'Wt', 'U_hat', 'V_hat', 'params');
            end
        end
    end
end