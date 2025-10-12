%% 检测数据集是否满秩。太棒了，怎么都不满秩啊！
clear; clc;
for i = 1:10
    for k = 12%7:12
        path = ['F:/EEEMG/stICA/Data/simulation/datasets' num2str(i) '/UV_compo' num2str(k) '_NMFm.mat'];
        % path = ['G:/stICA_Data/Data/simulation/datasets' num2str(i) '/UV_compo' num2str(k) '_NMFm.mat'];
        load(path);
        r = rank(U_hat);
        if k==r
            disp('满秩');
        else
            disp(['不满秩,rank=' num2str(r) ',set=' num2str(i) ',k=' num2str(k)]);
        end
    end
end