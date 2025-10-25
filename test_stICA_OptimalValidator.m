% 适用于stICA迭代结果的优化结果测试
% 检测最优结果周围是否有更优结果
addpath('F:/EEEMG/stICA/Func');
% params.NMFt = 0;
% pts = -3 : 0.005 : 3;
% for i = 1:size(Xt,2)
%     [f, x] = ksdensity(Xt(:,i), pts);
%     pp = spline(x, f);
%     pp_d = fnder(pp, 1);
%     params.pp(i) = pp;
%     params.pp_d(i) = pp_d;
% end
fun = @(Ws) f_function(Ws,U_hat,V_hat,params);
fval_CG = fun(Ws);
for i = 1:1000
    % 这一步的步长越小，能产生的更优结果越多。
    % 但是同样的，步长越小，目标函数值的变化也越小。
    % 第一步：调整步长范围，看在哪个范围内有更优点。
    % 步长小于1e-6时，即使有更优结果，也基本可以认为优化到局部最优点。
    % 步长大于1e-6时，若有更优结果，再看函数变化量的大小。
    Ws_test = Ws+randn(params.k)*1e-3;
    fval_test = fun(Ws_test);
    if fval_test < fval_CG
        disp(['有更优的结果,i=' num2str(i)]);
        % 第二步：调整目标函数值变化范围。
        % 如果函数变化很小，则可以认为已经收敛。
        % 如果函数变化值大于1e-3，基本可以认为还未收敛。
        if abs(fval_test - fval_CG) > 1e-3
            disp(num2str(abs(fval_test - fval_CG)));
            break
        end
    end
end