function [f, g] = f_function(Ws,U_hat,V_hat,params)

Ws = reshape(Ws,params.k,params.k); %~ 如果采用粒子群优化算法需要用到
Wt = inv(Ws');

% 计算时间熵
if params.NMFt
    %~ 计算时间熵，无分布区别，采用sech2函数，不带偏度的单一分布概率密度函数
    y_t = V_hat*Wt;
    p_t = sech(y_t).^2;
    L_t = sum(mean(log(p_t+eps)));
    DET_t = log(abs(det(Wt)));
    f_t = -(L_t + DET_t);

else
    % 使用Xt的数据插值得到pdf
    y_t = V_hat*Wt;
    for i = 1:12
        p_t(:, i) = ppval(params.pp(i), y_t(:, i));
        idx1 = find(y_t(:, i) < params.pp(i).breaks(1));
        p_t(idx1, i) = 0;
        idx2 = find(y_t(:, i) > params.pp(i).breaks(end));
        p_t(idx2, i) = 0;
    end
    L_t = sum(mean(log(p_t+eps)));
    DET_t = log(abs(det(Wt)));
    f_t = -(L_t + DET_t);
end

% 计算空间熵
if params.skew_stICA
    if params.MultiDistribution
        %~ 带偏度的多分布概率密度函数
        p = params.p;
        a = params.a;
        b = params.b;
        u = params.u;

        y_s = U_hat*Ws-u;
        s_s = sqrt(1+y_s.^2);
        A_s = a.*y_s;
        B_s = b.*s_s;
        logp_s = A_s-B_s;

        L_s = sum(mean(logp_s+log(p+eps)));
        DET_s = log(abs(det(Ws)));
        f_s = -(L_s + DET_s);

    else
        %~ 带偏度的单一分布概率密度函数
        p = 1e3;
        a = 15;
        b = 16.5431;

        y_s = U_hat*Ws;
        s_s = sqrt(1+y_s.^2);
        A_s = a.*y_s;
        B_s = b.*s_s;
        logp_s = A_s-B_s;
        L_s = sum(mean(logp_s+log(p+eps)));
        DET_s = log(abs(det(Ws)));
        f_s = -(L_s + DET_s);
    end
else
    %~ 不带偏度的单一分布概率密度函数
    y_s = U_hat*Ws;
    p_s = sech(y_s).^2;
    L_s = sum(mean(log(p_s+eps)));
    DET_s = log(abs(det(Ws)));
    f_s = -(L_s + DET_s);
end

%～ 计算总熵
f = params.alpha*f_s + (1-params.alpha)*f_t;

% 梯度
if nargout > 1
    % 时间熵的梯度
    % sigma_d = sech(y_t).^2;
    % sigma_dd = -2*sech(y_t).^3.*sinh(y_t);
    tmp = -2*sech(y_t)'.*sinh(y_t)';
    sigma = reshape(tmp, params.k, 1, []);
    x = reshape(V_hat', 1, params.k, []);
    % phi_y = mean(-2*sech(y_t).*sinh(y_t));
    phi_yt = mean(sigma .* x, 3);
    % x = mean(V_hat);
    g_t = -(Ws + phi_yt);

    % 用Xt拟合
    % for i = 1:12
    %     dp_t(:,i) = ppval(params.pp_d(i), y_t(:,i));
    % end
    % phi_y = mean(dp_t ./ p_t);
    % % x = mean(V_hat);
    % g_t = Ws + phi_y' * mean(V_hat);

    % 空间熵的梯度
    if params.skew_stICA
        % 使用具有偏度的多分布概率密度函数
        % sigma_d = p .* exp(a.*y_s - b.*sqrt(y_s.^2+1));
        % sigma_dd = sigma_d .* (a - b.*y_s./sqrt(y_s.^2+1));
        tmp = (a - b.*y_s./sqrt(y_s.^2+1))';
        sigma = reshape(tmp, params.k, 1, []);
        x = reshape(U_hat', 1, params.k, []);
        phi_ys = mean(sigma .* x, 3);
        % x = mean(U_hat);
        g_s = -(Wt + phi_ys);
    else
        % 使用不具有偏度的双曲正割函数
        % sigma_d = sech(y_s).^2;
        % sigma_dd = -2 * sech(y_s).^3 .* sinh(y_s);
        phi_ys = mean(-2*sech(y_s).*sinh(y_s));
        % x = mean(U_hat);
        g_s = -(Wt + phi_ys' * mean(U_hat));
    end

    g = params.alpha*g_s + (1-params.alpha)*g_t;
end
