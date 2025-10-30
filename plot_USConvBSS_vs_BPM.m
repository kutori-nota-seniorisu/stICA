%% 将USConvBSS的结果与BPM进行对比
figure;
subplot(2,1,1);
hold on;
for i=1:10
    % load(['Results/result' num2str(i) '_Update.mat']);
    load(['Results/datasets' num2str(i) '_result.mat']);
    boxplot(matchresult_time.time, 'Positions', i-0.2, 'Colors', 'r');
    hold on;
    timeMean(i) = mean(matchresult_time.time);
    numMU(i) = size(matchresult_time, 1);
end
plot((1:10)-0.2, timeMean, 'r-*');
hold on;

% for i=1:10
%     % load(['Results/result' num2str(i) '_Update.mat']);
%     load(['Results/datasets' num2str(i) '_resultnew.mat']);
%     boxplot(matchresult_time.time, 'Positions', i, 'Colors', 'b');
%     hold on;
%     timeMean2(i) = mean(matchresult_time.time);
%     numMU2(i) = size(matchresult_time, 1);
% end
% plot((1:10), timeMean2, 'b-*');
% hold on;

load('Results/compo12_NC_NMFm0.9_BPM_10sets.mat');
for i = 1:10
    % space_mean_BPM(i) = mean(matchresult_final_all{i}.space);
    time_mean_BPM(i) = mean(matchresult_final_all{i}.time);
    No_BPM(i) = size(matchresult_final_all{i}, 1);
    % time_No_NMFm(i) = size(matchresult_time_all{i}, 1);
    % space_No_NMFm(i) = size(matchresult_space_all{i}, 1);
end
match_BPM = matchresult_final_all;
for i = 1:10
    boxplot(match_BPM{i}.time, 'Positions', i+0.2, 'Colors', 'g');
    hold on;
end
plot((1:10)+0.2,time_mean_BPM,'g-*');
ylabel('RoA (Time)');
xticks(1:10)
xticklabels(1:10);
legend('USCBSS', 'BPM');

subplot(2,1,2);
plot(1:10, numMU, 'r-*');
hold on;
% plot(1:10, numMU2, 'b-*');
% hold on;
plot(1:10, No_BPM, 'g-*');
ylabel('No. MU');
legend('USCBSS', 'BPM');

sgtitle('USCBSS vs BPM')
set(gcf,'unit','normalized','position',[0.05,0.1,0.9,0.6]);