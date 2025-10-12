his = struct();
his.iteration = [];
his.funccount = [];
his.fval = [];
his.stepsize = [];
his.optimality = [];
his.x = []; % 如果需要也可以保存迭代点

%% 迭代中：将本次迭代的数据添加到history中
his.iteration = [his.iteration; history.iteration];
his.funccount = [his.funccount; history.funccount];
his.fval = [his.fval; history.fval];
his.stepsize = [his.stepsize; history.stepsize];
his.optimality = [his.optimality; history.optimality]; % 注意字段名
his.x = [his.x; history.x]; % 将列向量转置为行向量存储

%% 保存
history = his;
save('history.mat','history');