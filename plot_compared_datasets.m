%% 计算空间的corr和时间分量的ROA（spike匹配度）
clear; clc; 
%% 导入数据
algo_name = 'NMFm';
% algo_name = 'NMFn';
% algo_name = 'SVDn';
accuracy_matrix_space = {};
accuracy_matrix_time = {};
% 10个datasets，6种潜在成分个数设置(7到12)
mean_accuracy_matrix_space = zeros(10,6);
mean_accuracy_matrix_time = zeros(10,6);

compo_num = 12:-1:7;
for m = 1%:length(compo_num)
    for n = 1%:10
        datasets_num = num2str(n);
        file_name = ['UV_compo' num2str(compo_num(m)) '_' algo_name '.mat'];
        time_tolerance = 60;%对应±30ms的误差，都算spike train能对应上
        % file_path = ['D:\Code\Verasonics\运动单位解码\stICA_simple\Data\simulation\datasets' datasets_num '\'];
        file_path = ['F:/EEEMG/stICA_simple/Data/simulation/datasets' datasets_num '/'];
        load([file_path file_name])
        % pulses_path = ['D:\Code\Verasonics\运动单位解码\stICA_simple\Data\simulation\MU_time_response\TimeCompoDatasets' datasets_num '\'];
        pulses_path = ['F:/EEEMG/stICA_simple/Data/simulation/MU_time_response/TimeCompoDatasets' datasets_num '/'];
        load([pulses_path 'ipulses.mat'])
        pulses_2s = {};
        for i = 1:length(ipulses)
            index = ipulses{i} <= 4000;
            pulses_2s{end+1} = ipulses{i}(index);
        end

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 计算空间分量的匹配相关性
        S_reshape = reshape(S,400,128,size(S,2));
        Xs_reshape = reshape(Xs,400,128,size(Xs,2));

        for i = 1:size(S,2)
            for j = 1:size(Xs_reshape,3)
                corr_matrix_space(i,j) = corr2(S_reshape(:,:,i),Xs_reshape(:,:,j));
            end
        end
        [max_values_space, indices_space] = max(abs(corr_matrix_space),[],2);

        if length(max_values_space)==12
            % mink(A,K)返回A中最小的K个元素
            [minValue_space, minIndex_space] = mink(max_values_space,2);
            % 去除相关系数最小的
            max_values_space_withoutMin = max_values_space(~ismember(max_values_space,minValue_space));
        else
            [minValue_space, minIndex_space] = min(max_values_space);
            max_values_space_withoutMin = max_values_space(max_values_space ~= minValue_space);
        end

        mean(max_values_space);

        max_values_space_withoutMin;
        mean(max_values_space_withoutMin);

        accuracy_matrix_space{m,1}(:,n) = max_values_space_withoutMin;
        mean_accuracy_matrix_space(n,m) = mean(max_values_space_withoutMin);


        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 计算时间分量的spike train匹配度
        % 滤波（可以调节参数看效果）
        T_smooth = smoothdata(T,"movmedian",8);
        decompo_pulses = {};
        for i = 1:size(T_smooth,2)
            [~, locs] = findpeaks(T_smooth(:,i),'MinPeakDistance',90); % 最后一个参数视情况调节，对匹配准确度影响较大
            decompo_pulses{end+1} = locs';
        end
        % 绘制spike train
        % plotDecomps(pulses_2s,[],2000,0,0,[]);
        % plotDecomps(decompo_pulses,[],2000,0,0,[]);

        % spike匹配ROA
        for i = 1:length(decompo_pulses)
            for j = 1:length(pulses_2s)
                % 创建一个包含所有组合的矩阵
                [Array1, Array2] = meshgrid(decompo_pulses{i}, pulses_2s{j});

                % 计算差值的绝对值
                diff_values = abs(Array1 - Array2);

                % 找到满足条件的元素
                valid_elements = diff_values < time_tolerance;

                % 计算满足条件的总数目
                count = sum(valid_elements(:));

                spike_ROA_matrix(i,j) = count/size(pulses_2s{j},2);
            end
        end

        % 去除指定索引位置的值
        indices_time = indices_space;
        indices_time(minIndex_space) = [];
        indicesToRemove = find(indices_time > 10);%索引超过10的去除，因为11和12是高斯噪声图像
        indices_time(indicesToRemove) = 1;

        for i = 1:size(indices_time,1)
            spike_ROA_temp = spike_ROA_matrix(i,indices_time(i));
            spike_ROA_temp(spike_ROA_temp>1) = 1;
            spike_ROA(i) = spike_ROA_temp;
        end
        spike_ROA(indicesToRemove) = NaN;

        spike_ROA';
        mean(spike_ROA);

        accuracy_matrix_time{m,1}(:,n) = spike_ROA';

        spike_ROA(indicesToRemove) = [];
        mean_accuracy_matrix_time(n,m) = mean(spike_ROA);

        % clearvars -except algo_name accuracy_matrix_space accuracy_matrix_time mean_accuracy_matrix_space mean_accuracy_matrix_time compo_num m n
    end

end