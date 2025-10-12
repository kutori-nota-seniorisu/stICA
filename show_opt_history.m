%% 导入优化历史后，输出收敛指标（共轭梯度法）
disp(['结束时函数值变化量:' num2str(history.fval(end-1)-history.fval(end))]);
Ws_iter = reshape(history.x, 12, 12, []);
a = norm(Ws_iter(:,:,end-1)-Ws_iter(:,:,end), 'fro');
b = norm(Ws_iter(:,:,end), 'fro');
disp(['结束时自变量的相对变化量:' num2str(a/b*100) '%']);
disp(['结束时一阶最优性条件:' num2str(history.optimality(end))]);