function muap_data = interpolateMUAP(muap_data)
% 对MUAP阵列中的NaN进行插值
% 输入一个MUAP阵列的数组，而不是元胞

% 获取数组大小
[n_rows, n_cols, n_time] = size(muap_data);

% 创建坐标网格
[x_grid, y_grid] = meshgrid(1:n_cols, 1:n_rows);

% 对每个时间点进行插值
for t = 1:n_time
    current_slice = muap_data(:, :, t);

    % 找出NaN的位置
    nan_mask = isnan(current_slice);

    % 如果没有NaN，跳过
    if ~any(nan_mask(:))
        continue;
    end

    % 提取已知点（非NaN）
    known_mask = ~nan_mask;

    % 已知点的坐标和值
    x_known = x_grid(known_mask);
    y_known = y_grid(known_mask);
    z_known = current_slice(known_mask);

    % 需要插值的点
    x_nan = x_grid(nan_mask);
    y_nan = y_grid(nan_mask);

    % 使用griddata进行插值
    if ~isempty(x_nan) && length(z_known) >= 4
        try
            z_interp = griddata(x_known, y_known, z_known, ...
                x_nan, y_nan, 'cubic');

            % 将插值结果放回数组
            for k = 1:length(x_nan)
                muap_data(y_nan(k), x_nan(k), t) = z_interp(k);
            end
        catch
            % 如果插值失败，使用最近邻插值
            z_interp = griddata(x_known, y_known, z_known, ...
                x_nan, y_nan, 'nearest');
            for k = 1:length(x_nan)
                muap_data(y_nan(k), x_nan(k), t) = z_interp(k);
            end
        end
    end
end