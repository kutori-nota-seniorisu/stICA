function time_curves=get_single_sta(actarea,tf2)
%% 用于得到muta中每个点的curve值，输出位置和数据
% input: 
%   actarea：逻辑区域
%   tf2：STA后的区域


% 找到逻辑为1的位置的线性索引
[row, col] = find(actarea); % 获取逻辑为1的行和列索引
linear_indices = sub2ind(size(actarea), row, col); % 将行和列索引转换为线性索引
% 初始化存储提取数据的数组
num_points = length(linear_indices); % 逻辑为1的点的数量
% extracted_data = zeros(num_points, length(T)); % 初始化提取数据的存储数组
% 循环提取每个时间点的数据
for i = 1:length(row)
    extracted_data(i,: ) = tf2(row(i),col(i), :); % 提取每个时间点逻辑为1的数据
end
% 创建一个结构数组来存储位置和数据
time_curves = struct('Position', [], 'Data', []);
for i = 1:num_points
    time_curves(i).Position = [row(i), col(i)]; % 存储位置
    time_curves(i).Data = extracted_data(i, :); % 存储数据
end
end