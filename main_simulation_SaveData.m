%% 保存变量，在运行完main_simulation_multiChannel_iEMG后执行
savepath = ['F:/EEEMG/stICA/Data/simulation/datasets' datasets_num ];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end
% 保存生成的TVI数据
% save([savepath '/TVIdata_compo' num2str(params.k) '_NC.mat'], 'TVIdata');
% 保存分解结果
% save([savepath '/' load_NMF_data_name '.mat'], 'S','T','U_hat','V_hat','Xs','Xt','Ws0','Ws','Wt','params');
save([savepath '/UV_compo12_NC_NMFm0.1.mat'], 'S','T','U_hat','V_hat','Xs','Xt','Ws0','Ws','Wt','params');
disp('数据保存完成');