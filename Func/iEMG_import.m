clear;
addpath(genpath('D:/Code/Verasonics/运动单位解码/stICA_simple/Func/'))
fsampu = 2000; % 采样率
file_name = 'IEMG_SEMG_UU_10000_Trial1_1114_iPulses';
data = importdata(['D:/Code/Verasonics/运动单位解码/stICA_simple/Data/iEMG/' file_name '.eaf']); % 读取eaf文件
muNum = max(data.data(:,2)); % MU的个数
iPulses = {};
for mu = 1:muNum
    % iPulses就是这个eaf文件里分解得到的spike train，每个cell表示一个MU，里面的数字是该MU每次放电的时刻
    iPulses{mu} = round(data.data(find(data.data(:,2)==mu),1)'*fsampu); 
end


flag_feat = 0;
order = 0;%~ 是否按MU最开始激活的时间排序
plotDecomps(iPulses,[],fsampu,flag_feat,order,[]);