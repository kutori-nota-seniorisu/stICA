function stop = myOutputFcn2(x, optimValues, state, varargin)
% 输出函数，用于在每次迭代时保存数据
% stop: 标志，如果为true则停止优化
% x: 当前迭代点
% optimValues: 包含当前迭代信息的结构体
% state: 当前迭代状态（'init', 'iter', 'interrupt', 'done'）

% 声明一个持久变量（persistent variable）来在函数调用间存储数据
persistent history

% 初始化停止标志
stop = false;

if nargin > 3 && ischar(varargin{1})
    run_id = varargin{1};
else
    run_id = 'run_1';
end

if nargin > 4 && ischar(varargin{2})
    save_filename = varargin{2};
else
    save_filename = 'test_save.mat';
end

if isempty(history)
    history = struct();
end

if ~isfield(history, run_id)
    history.(run_id) = struct();
    history.(run_id).iteration = [];
    history.(run_id).funccount = [];
    history.(run_id).fval = [];
    history.(run_id).stepsize = [];
    history.(run_id).optimality = [];
    history.(run_id).gradient = [];
    history.(run_id).x = [];
end

current_history = history.(run_id);

% 根据迭代状态进行操作
switch state
    case 'init'
        % 初始化：创建一个空结构来存储历史记录
        fprintf('梯度下降运行 %s 开始...\n', run_id);
        % fprintf('初始点: %s\n', mat2str(x));
        fprintf('初始函数值: %.6f\n', optimValues.fval);
        % history = struct();
        % history.iteration = [];
        % history.funccount = [];
        % history.fval = [];
        % history.stepsize = [];
        % history.optimality = [];
        % history.gradient = [];
        % history.x = []; % 如果需要也可以保存迭代点

    case 'iter'
        % 迭代中：将本次迭代的数据添加到history中
        % current_history.iteration(end+1) = optimValues.iteration;
        % current_history.iteration(end+1) = optimValues.iteration;
        % current_history.funccount(end+1) = optimValues.funccount;
        % current_history.fval = [current_history.fval; optimValues.fval];
        % current_history.stepsize = [current_history.stepsize; optimValues.stepsize];
        % current_history.optimality = [current_history.optimality; optimValues.firstorderopt]; % 注意字段名
        % current_history.gradient = [current_history.gradient; optimValues.gradient];

        current_history.iteration = [current_history.iteration; optimValues.iteration];
        current_history.funccount = [current_history.funccount; optimValues.funccount];
        current_history.fval = [current_history.fval; optimValues.fval];
        current_history.stepsize = [current_history.stepsize; optimValues.stepsize];
        current_history.optimality = [current_history.optimality; optimValues.firstorderopt]; % 注意字段名
        current_history.gradient = [current_history.gradient; optimValues.gradient];
        current_history.x = [current_history.x; x];

        history.(run_id) = current_history;

    case 'done'
        save(save_filename, 'history');
end

end