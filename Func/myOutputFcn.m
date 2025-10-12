function stop = myOutputFcn(x, optimValues, state)
% 输出函数，用于在每次迭代时保存数据
% stop: 标志，如果为true则停止优化
% x: 当前迭代点
% optimValues: 包含当前迭代信息的结构体
% state: 当前迭代状态（'init', 'iter', 'interrupt', 'done'）

% 声明一个持久变量（persistent variable）来在函数调用间存储数据
persistent history

% 初始化停止标志
stop = false;

% 根据迭代状态进行操作
switch state
    case 'init'
        % 初始化：创建一个空结构来存储历史记录
        history = struct();
        history.iteration = [];
        history.funccount = [];
        history.fval = [];
        history.stepsize = [];
        history.optimality = [];
        history.gradient = [];
        history.x = []; % 如果需要也可以保存迭代点

    case 'iter'
        % 迭代中：将本次迭代的数据添加到history中
        history.iteration = [history.iteration; optimValues.iteration];
        history.funccount = [history.funccount; optimValues.funccount];
        history.fval = [history.fval; optimValues.fval];
        history.stepsize = [history.stepsize; optimValues.stepsize];
        history.optimality = [history.optimality; optimValues.firstorderopt]; % 注意字段名
        history.gradient = [history.gradient; optimValues.gradient];
        history.x = [history.x; x];

        % （可选）实时保存到mat文件，防止程序中断丢失数据
        % save('opt_history.mat', 'history');

    case 'done'
         % 优化完成：保存最终的历史记录
        history_file = dir(fullfile('./Log', '*.mat'));
        if isempty(history_file)
            save('./Log/cg_history_final1.mat', 'history');
        else
            save(['./Log/cg_history_final' num2str(size(history_file,1)+1) '.mat'], 'history');
        end
        % save('opt_history_final.mat', 'history');
        % 清除持久变量，释放内存
        clear history
        fprintf('最终历史记录已保存\n');
end

end