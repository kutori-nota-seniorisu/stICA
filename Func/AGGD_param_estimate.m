function [a, bl, br] = AGGD_param_estimate(vec)
% % 峰值的位置
% [~,m] = max(vec);
% 左右比例参数的估计
sl = sqrt(mean((vec(vec<0)).^2));
sr = sqrt(mean((vec(vec>=0)).^2));
% 比值
gamma_hat = sl/sr;

% 一阶绝对矩
m1_hat = mean(abs(vec));
% 二阶原点矩
u2_hat = mean(vec.^2);
% 比值
r_hat = m1_hat^2/u2_hat;

% ρ(α)的估计值
R_hat = r_hat*(gamma_hat^3+1)*(gamma_hat+1)/((gamma_hat^2+1)^2);
% 给定候选的形状参数
alpha   = 0:0.001:10;
alpha(1) = [];
% 计算所有可能的ρ(α)
rho_alpha = ((gamma(2./alpha)).^2)./(gamma(1./alpha).*gamma(3./alpha));
% 寻找距离最近的α
[~, idx] = min((rho_alpha-R_hat).^2);
a = alpha(idx);
bl = sl * sqrt(gamma(3/a)/gamma(1/a));
br = sr * sqrt(gamma(3/a)/gamma(1/a));