function stop = psoOutputFcn(optimValues, state, varargin)
% PSOOUTPUTFCN 粒子群算法输出函数，用于保存优化过程历史记录
% 
% 输入参数：
%   optimValues - 优化过程信息结构体
%   state       - 当前状态 ('init', 'iter', 'interrupt', 'done')
%   varargin    - 可选参数：历史记录结构体
%
% 输出参数：
%   stop - 是否停止优化 (true/false)

    persistent history % 持久变量保存历史记录
    
    stop = false; % 默认不停止优化
    
    % 处理可选输入参数
    if nargin > 2 && isstruct(varargin{1})
        history = varargin{1};
    end
    
    % 初始化历史记录结构
    if isempty(history)
        history = struct();
        history.iteration = [];
        history.globalBest = [];
        history.globalBestFval = [];
        history.meanFval = [];
        history.stdFval = [];
        % history.bestParticles = {};
        % history.allPositions = {};
        % history.allFvals = {};
        % history.timeStamps = [];
        history.gradientNorms = [];
        history.swarmDiversity = [];
        history.convergenceHistory = [];
    end
    
    switch state
        case 'init'
            % 初始化阶段
            fprintf('PSO 优化开始...\n');
            fprintf('种群大小: %d, 变量维度: %d\n', ...
                    size(optimValues.swarm, 1), size(optimValues.swarm, 2));

        case 'iter'
            % 迭代阶段 - 记录历史信息
            iter = optimValues.iteration;
            currentBestx = optimValues.bestx;
            currentBestFval = optimValues.bestfval;
            currentFunccount = optimValues.funccount;
            currentStall = optimValues.stalliterations;
            currentSwarm = optimValues.swarm;
            currentSwarmFvals = optimValues.swarmfvals;
            currentMeanFval = optimValues.meanfval;
            
            % 计算种群统计信息
            stdFval = std(currentSwarmFvals);
            
            % 计算种群多样性（平均距离）
            swarmDiversity = computeSwarmDiversity(currentSwarm);
            
            % 计算梯度范数（近似）
            if iter > 1
                gradNorm = norm(history.globalBest(:,end) - currentBestx);
            else
                gradNorm = NaN;
            end
            
            % 记录历史信息
            history.iteration(end+1) = iter;
            history.globalBest(:,end+1) = currentBestx;
            history.globalBestFval(end+1) = currentBestFval;
            history.meanFval(end+1) = currentMeanFval;
            history.stdFval(end+1) = stdFval;
            % history.bestParticles{end+1} = currentSwarm(1:min(5,end),:); % 保存前5个最佳粒子
            % history.allPositions{end+1} = currentSwarm;
            % history.timeStamps(end+1) = now;
            history.gradientNorms(end+1) = gradNorm;
            history.swarmDiversity(end+1) = swarmDiversity;
            
            % 记录收敛历史
            if iter > 1
                improvement = history.globalBestFval(end-1) - currentBestFval;
                history.convergenceHistory(end+1) = improvement;
            else
                history.convergenceHistory(end+1) = 0;
            end
            
            % 显示当前迭代信息
            % fprintf('%-8d %-12.6f %-12.6f %-12.6f %-12.6f\n', ...
            %         iter, currentBestFval, currentMeanFval, stdFval, swarmDiversity);
            
            % 每10次迭代保存一次完整数据
            % if mod(iter, 10) == 0
            %     save(sprintf('pso_history_iter_%d.mat', iter), 'history');
            %     fprintf('已保存第 %d 次迭代的历史记录\n', iter);
            % end
            
            % 检查收敛条件（可选）
            % if iter > 50 && std(history.globalBestFval(end-49:end)) < 1e-8
            %     fprintf('检测到收敛，停止优化\n');
            %     stop = true;
            % end
            
        case 'done'
            % 优化完成阶段
            fprintf('\nPSO 优化完成！\n');
            fprintf('最终最佳值: %.6f\n', optimValues.bestfval);
            fprintf('总迭代次数: %d\n', optimValues.iteration);
            fprintf('最佳解: %s\n', mat2str(optimValues.bestx, 4));
            
            % 保存最终历史记录
            history_file = dir(fullfile('./Log', '*.mat'));
            if isempty(history_file)
                save('./Log/pso_history_final1.mat', 'history');
            else
                save(['./Log/pso_history_final' num2str(size(history_file,1)+1) '.mat'], 'history');
            end

            % save('pso_final_history.mat', 'history');
            fprintf('最终历史记录已保存\n');
            
            % 绘制收敛曲线
            % plotConvergence(history);
    end
    
    % 将历史记录保存到基础工作区（可选）
    assignin('base', 'psoHistory', history);
end

function diversity = computeSwarmDiversity(swarm)
% 计算种群多样性 - 粒子间的平均欧几里得距离
    nParticles = size(swarm, 1);
    center = mean(swarm, 1);
    distances = sqrt(sum((swarm - center).^2, 2));
    diversity = mean(distances);
end

function plotConvergence(history)
% 绘制收敛曲线和其他分析图表
    if isempty(history.iteration)
        return;
    end
    
    ax = figure('Position', [100, 100, 1200, 800]);
    
    % 1. 最佳值收敛曲线
    subplot(2,3,1);
    semilogy(history.iteration, history.globalBestFval, 'b-', 'LineWidth', 2);
    hold on;
    plot(history.iteration, history.meanFval, 'r--', 'LineWidth', 1.5);
    title('适应度值收敛曲线');
    xlabel('迭代次数');
    ylabel('适应度值 (对数尺度)');
    legend('全局最佳', '种群平均', 'Location', 'best');
    grid on;
    
    % 2. 种群多样性
    subplot(2,3,2);
    plot(history.iteration, history.swarmDiversity, 'g-', 'LineWidth', 2);
    title('种群多样性');
    xlabel('迭代次数');
    ylabel('平均距离');
    grid on;
    
    % 3. 梯度范数变化
    subplot(2,3,3);
    semilogy(history.iteration(2:end), history.gradientNorms(2:end), 'm-', 'LineWidth', 2);
    title('梯度范数变化');
    xlabel('迭代次数');
    ylabel('梯度范数 (对数尺度)');
    grid on;
    
    % 4. 适应度值标准差
    subplot(2,3,4);
    plot(history.iteration, history.stdFval, 'c-', 'LineWidth', 2);
    title('适应度值标准差');
    xlabel('迭代次数');
    ylabel('标准差');
    grid on;
    
    % 5. 收敛历史
    subplot(2,3,5);
    plot(history.iteration, history.convergenceHistory, 'k-', 'LineWidth', 2);
    title('每次迭代的改进量');
    xlabel('迭代次数');
    ylabel('改进量');
    grid on;

    % % 6. 最终种群分布（2D示例）
    % if size(history.globalBest, 1) >= 2
    %     subplot(2,3,6);
    %     finalPositions = history.allPositions{end};
    %     scatter(finalPositions(:,1), finalPositions(:,2), 50, history.allFvals{end}, 'filled');
    %     hold on;
    %     plot(history.globalBest(1,end), history.globalBest(2,end), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    %     colorbar;
    %     title('最终种群分布');
    %     xlabel('维度 1');
    %     ylabel('维度 2');
    %     grid on;
    % end
    
    sgtitle('PSO 优化过程分析', 'FontSize', 14, 'FontWeight', 'bold');

    png_file = dir(fullfile('./Log', '*.png'));
    if isempty(history_file)
        saveas(ax, './Log/pso_history1', 'png');
    else
        aveas(ax, ['./Log/pso_history' num2str(size(png_file,1)+1)], 'png');
    end
end