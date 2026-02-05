%%% 针对 Peak2Peak 输出结构，绘制单探头、单方法（P2P）下每个 MU 的聚类区域与平均曲线
% 需要step1 和step2结果
% 输出每个MU的结果图到savefolder中
clearvars -except map_all curvesofmap_all meantwitch_all motionnum levelnum trialnum saveflag savefile threshold Sub property varFrames_all twitchFrames_all

methodName = 'P2P';  % 固定方法名
motionnum = 1;
[levelnum, trialnum] = size(map_all{motionnum});
Sub=1;
saveflag=1;
savefolder='Z:\Result\25-07-17前臂屈肌肌电超声\MUarea\median\ilim132jlim100';
% savefolder='D:\ICdata\muarea';

drawflag = true;
if drawflag
    for M = motionnum
        savefile=[savefolder  '\M' num2str(M)];%
        % loadfile=['Z:\Result\25-07-04\USSTA_5-35Hz_M' num2str(M) '_0s_15s_-49_to_100_Frame_2000Hz'];
        % load(loadfile)
        for L = 1:levelnum
            for T = 1:trialnum
                mapset = map_all{M}{L,T};
                varset=varFrames_all{M}{L,T};
                propertyuse=property{M}{L,T}.tfpp;
                % propertyuse=property{M}{L,T}.emma;
                if isempty(mapset), continue; end
                munum = size(mapset,2);
                for mu = 1:munum
                    % if ~isfield(mapset{1,mu}, methodName), continue; end
                    areaMap = mapset{1,mu}> 0;
                    varMap = mean(varset{mu},3);
                    curves = curvesofmap_all{M}{L,T}{1,mu};
                    meanCurves = meantwitch_all{M}{L,T}{1,mu};

                    if isempty(curves) || isempty(meanCurves), continue; end

                    % 构建聚类编号图
                    heatmap = zeros(size(areaMap));
                    for i = 1:numel(curves)
                        pos = ceil(curves(i).Position);
                        heatmap(pos(1), pos(2)) = curves(i).group;
                    end

                    % 绘图
                    figure('Visible','off','Name', sprintf('M%dL%dT%d MU%d P2P', M,L,T,mu));
                    set(gcf, 'Units','centimeters', ...
         'Position',[5, 5, 15, 15]); 
                    tiledlayout(2,2);
                

                    % 区域图
                    nexttile(1);
                    imagesc(heatmap);
                    % axis image off;
                    title(sprintf('M%dL%dT%d MU%d 区域 (P2P)', M,L,T,mu));

                    % 曲线图
                    nexttile(2);
                    hold on;
                    for k = 1:numel(curves)
                        plot(curves(k).Data, 'Color', [0.7,0.7,0.7], 'LineWidth', 0.3);
                    end
                    colors = lines(numel(meanCurves));
                    for g = 1:numel(meanCurves)
                        plot(meanCurves{g}, 'LineWidth', 2, 'Color', colors(g,:));
                    end
                    hold off;
                    xlabel('时间帧'); ylabel('响应速度');
                    title('Mean Twitch - P2P');
                    
                      %MAP
                    nexttile(3)
                    imagesc(propertyuse{mu}{1})
                    title('特征值')
                    
                    %Var
                    nexttile(4);
                    imagesc(varMap)
                    title('mean variation')
                    colorbar

                    sgtitle(sprintf('Sub%d M%dL%dT%d MU%d', Sub, M, L, T, mu));

                    if exist('saveflag','var') && saveflag
                        if ~exist(savefile,'dir'), mkdir(savefile); end
                        fn = sprintf('Sub%d_M%dL%dT%d_MU%d_%s.png', Sub, M, L, T, mu,methodName);
                        saveas(gcf, fullfile(savefile, fn));
                    end
                    close(gcf);
                end
            end
        end
    end
end