function evaluate_fit_results(X_data, K_data, K_sym, rank, show_plots)
%   评估多项式拟合结果
%
%   输入参数:
%   X_data - 自变量数据，格式为 [total_samples × n_vars]
%   K_data - 原始K矩阵数据，格式为 [total_samples × N_u × N_x]
%   K_sym - 拟合出的符号表达式矩阵，[N_u × N_x]符号矩阵
%   rank - 拟合阶数
%   show_plots - 是否显示绘图，默认为true

    if nargin < 5
        show_plots = true;
    end
    
    fprintf("\n=== 正在评估拟合结果 ===\n");
    
    %% 基本参数和输入检查
    dims = size(K_data);
    total_samples = dims(1);
    N_u = dims(2);
    N_x = dims(3);
    poly_n = N_u * N_x;
    n_vars = size(X_data, 2);
    
    % 检查输入数据维度
    if total_samples ~= size(X_data, 1)
        error('X_data 和 K_data 的样本数量不一致');
    end
    
    fprintf('数据维度: %d 个样本, %d 个自变量\n', total_samples, n_vars);
    fprintf('K矩阵维度: %d × %d\n', N_u, N_x);
    fprintf('拟合阶数: %d\n', rank);
    
    %% 计算预测值和误差
    MSE = zeros(N_u, N_x);  % 均方误差
    MAE = zeros(N_u, N_x);  % 平均绝对误差
    MAX_AE = zeros(N_u, N_x);  % 最大绝对误差
    R_squared = zeros(N_u, N_x);  % 决定系数R²
    
    fprintf('计算拟合误差...\n');
    
    for i = 1:N_u
        for j = 1:N_x
            % 获取原始数据
            K_actual = squeeze(K_data(:, i, j));
            
            % 从符号表达式中提取变量名
            vars = symvar(K_sym(i, j));
            
            % 检查变量数量
            if length(vars) ~= n_vars
                warning('K(%d,%d): 符号表达式有 %d 个变量, 但 X_data 有 %d 列。尝试继续计算...', ...
                    i, j, length(vars), n_vars);
                % 尝试使用通用的变量名
                vars = sym(sprintf('x%d', 1:n_vars));
            end
            
            % 计算预测值（使用符号表达式）
            K_pred = zeros(total_samples, 1);
            for k = 1:total_samples
                % 动态替换变量
                subs_expr = K_sym(i, j);
                for v_idx = 1:min(length(vars), n_vars)
                    subs_expr = subs(subs_expr, vars(v_idx), X_data(k, v_idx));
                end
                K_pred(k) = double(subs_expr);
            end
            
            % 计算误差指标
            residuals = K_actual - K_pred;
            
            % 均方误差
            MSE(i, j) = mean(residuals.^2);
            
            % 平均绝对误差
            MAE(i, j) = mean(abs(residuals));
            
            % 最大绝对误差
            MAX_AE(i, j) = max(abs(residuals));
            
            % 决定系数R²
            SS_res = sum(residuals.^2);
            SS_tot = sum((K_actual - mean(K_actual)).^2);
            if SS_tot > 0
                R_squared(i, j) = 1 - SS_res / SS_tot;
            else
                R_squared(i, j) = 1;
            end
        end
    end
    
    %% 输出整体统计
    fprintf('拟合误差统计：\n');
    fprintf('==============\n');
    
    fprintf('整体均方误差 (MSE):\n');
    fprintf('  平均值: %.4e\n', mean(MSE(:)));
    [min_val, min_idx] = min(MSE(:));
    [i_min, j_min] = ind2sub(size(MSE), min_idx);
    fprintf('  最小值: %.4e (位置: [%d,%d])\n', min_val, i_min, j_min);
    [max_val, max_idx] = max(MSE(:));
    [i_max, j_max] = ind2sub(size(MSE), max_idx);
    fprintf('  最大值: %.4e (位置: [%d,%d])\n', max_val, i_max, j_max);
    fprintf('  中位数: %.4e\n', median(MSE(:)));
    
    fprintf('\n整体平均绝对误差 (MAE):\n');
    fprintf('  平均值: %.4e\n', mean(MAE(:)));
    [min_val, min_idx] = min(MAE(:));
    [i_min, j_min] = ind2sub(size(MAE), min_idx);
    fprintf('  最小值: %.4e (位置: [%d,%d])\n', min_val, i_min, j_min);
    [max_val, max_idx] = max(MAE(:));
    [i_max, j_max] = ind2sub(size(MAE), max_idx);
    fprintf('  最大值: %.4e (位置: [%d,%d])\n', max_val, i_max, j_max);
    
    fprintf('\n整体最大绝对误差:\n');
    fprintf('  平均值: %.4e\n', mean(MAX_AE(:)));
    fprintf('  最小值: %.4e\n', min(MAX_AE(:)));
    [max_val, max_idx] = max(MAX_AE(:));
    [i_max, j_max] = ind2sub(size(MAX_AE), max_idx);
    fprintf('  最大值: %.4e (位置: [%d,%d])\n', max_val, i_max, j_max);
    
    fprintf('\n决定系数 R²:\n');
    fprintf('  平均值: %.4f\n', mean(R_squared(:)));
    [min_val, min_idx] = min(R_squared(:));
    [i_min, j_min] = ind2sub(size(R_squared), min_idx);
    fprintf('  最小值: %.4f (位置: [%d,%d])\n', min_val, i_min, j_min);
    [max_val, max_idx] = max(R_squared(:));
    [i_max, j_max] = ind2sub(size(R_squared), max_idx);
    fprintf('  最大值: %.4f (位置: [%d,%d])\n', max_val, i_max, j_max);
    
    fprintf('  R² < 0.9 的元素数: %d/%d (%.1f%%)\n', ...
        sum(R_squared(:) < 0.9), poly_n, sum(R_squared(:) < 0.9)/poly_n*100);
    fprintf('  R² < 0.95 的元素数: %d/%d (%.1f%%)\n', ...
        sum(R_squared(:) < 0.95), poly_n, sum(R_squared(:) < 0.95)/poly_n*100);
    
    fprintf('==============\n\n');
    
    %% 绘制误差分布图
    if show_plots
        fprintf('生成误差分析图...\n');
        
        figure('Position', [100, 100, 933, 533]);
        
        % 子图1: R²分布直方图
        subplot(2, 3, 1);
        histogram(R_squared(:), 20);
        xlabel('决定系数 R²');
        ylabel('频数');
        title('R² 分布直方图');
        grid on;
        hold on;
        xline(0.9, 'r--', 'R²=0.9', 'LineWidth', 1.5);
        xline(0.95, 'g--', 'R²=0.95', 'LineWidth', 1.5);
        hold off;
        
        % 子图2: MAE分布直方图
        subplot(2, 3, 2);
        histogram(MAE(:), 20);
        xlabel('平均绝对误差 (MAE)');
        ylabel('频数');
        title('MAE 分布直方图');
        grid on;
        
        % 子图3: MSE分布直方图
        subplot(2, 3, 3);
        histogram(log10(MSE(:)), 20);
        xlabel('log_{10}(MSE)');
        ylabel('频数');
        title('MSE 对数分布直方图');
        grid on;
        
        % 子图4: R²热力图
        subplot(2, 3, 4);
        imagesc(R_squared);
        colorbar;
        xlabel('K矩阵列索引 (j)');
        ylabel('K矩阵行索引 (i)');
        title('R² 热力图');
        axis equal tight;
        
        % 子图5: 某个K元素的拟合效果示例
        subplot(2, 3, 5);
        % 选择R²最高的元素作为示例
        [max_r2, max_idx] = max(R_squared(:));
        [i_best, j_best] = ind2sub(size(R_squared), max_idx);
        
        K_actual_best = squeeze(K_data(:, i_best, j_best));
        
        % 获取变量名
        vars = symvar(K_sym(i_best, j_best));
        if length(vars) ~= n_vars
            vars = sym(sprintf('x%d', 1:n_vars));
        end
        
        % 计算预测值
        K_pred_best = zeros(total_samples, 1);
        for k = 1:total_samples
            subs_expr = K_sym(i_best, j_best);
            for v_idx = 1:min(length(vars), n_vars)
                subs_expr = subs(subs_expr, vars(v_idx), X_data(k, v_idx));
            end
            K_pred_best(k) = double(subs_expr);
        end
        
        % 绘制散点图
        scatter(K_actual_best, K_pred_best, 20, 'filled');
        hold on;
        % 绘制y=x参考线
        min_val = min([K_actual_best; K_pred_best]);
        max_val = max([K_actual_best; K_pred_best]);
        plot([min_val, max_val], [min_val, max_val], 'r--', 'LineWidth', 1.5);
        hold off;
        
        xlabel('实际值');
        ylabel('预测值');
        title(sprintf('拟合效果示例 (K(%d,%d), R²=%.4f)', i_best, j_best, max_r2));
        grid on;
        axis equal;
        
        % 子图6: 根据自变量数量选择不同的可视化方式
        subplot(2, 3, 6);
        
        if n_vars == 1
            % 单变量情况：绘制拟合曲线
            x_range = linspace(min(X_data(:,1)), max(X_data(:,1)), 100);
            K_curve = zeros(size(x_range));
            
            for k = 1:length(x_range)
                subs_expr = K_sym(i_best, j_best);
                subs_expr = subs(subs_expr, vars(1), x_range(k));
                K_curve(k) = double(subs_expr);
            end
            
            plot(x_range, K_curve, 'b-', 'LineWidth', 2);
            hold on;
            scatter(X_data(:,1), K_actual_best, 50, 'r', 'filled');
            hold off;
            
            xlabel('自变量 x1');
            ylabel(sprintf('K(%d,%d)', i_best, j_best));
            title('单变量拟合曲线');
            grid on;
            legend('拟合曲线', '原始数据', 'Location', 'best');
            
        elseif n_vars == 2
            % 双变量情况：绘制拟合曲面
            x1_range = linspace(min(X_data(:,1)), max(X_data(:,1)), 30);
            x2_range = linspace(min(X_data(:,2)), max(X_data(:,2)), 30);
            [X1_grid, X2_grid] = meshgrid(x1_range, x2_range);
            
            % 计算拟合曲面
            K_surf = zeros(size(X1_grid));
            for m = 1:size(X1_grid, 1)
                for n = 1:size(X1_grid, 2)
                    subs_expr = K_sym(i_best, j_best);
                    subs_expr = subs(subs_expr, vars(1), X1_grid(m, n));
                    subs_expr = subs(subs_expr, vars(2), X2_grid(m, n));
                    K_surf(m, n) = double(subs_expr);
                end
            end
            
            surf(X1_grid, X2_grid, K_surf);
            hold on;
            % 绘制原始采样点
            scatter3(X_data(:,1), X_data(:,2), K_actual_best, 50, 'r', 'filled');
            hold off;
            
            xlabel('自变量 x1');
            ylabel('自变量 x2');
            zlabel(sprintf('K(%d,%d)', i_best, j_best));
            title('双变量拟合曲面');
            grid on;
            
        else
            % 多变量情况：显示变量重要性或其他信息
            text(0.3, 0.5, sprintf('自变量数量: %d\n无法直接可视化', n_vars), ...
                'FontSize', 12, 'HorizontalAlignment', 'center');
            title('多变量拟合结果');
            axis off;
        end
        
        sgtitle(sprintf('多项式拟合评估 (阶数: %d, 元素总数: %d, 自变量: %d)', rank, poly_n, n_vars));
        
        fprintf('误差分析图已生成。\n');
    end
    
    %% 生成简要报告
    fprintf('\n拟合质量简要报告：\n');
    fprintf('==============\n');
    
    % 根据R²平均值评估拟合质量
    mean_r2 = mean(R_squared(:));
    if mean_r2 >= 0.99
        quality = '优秀';
    elseif mean_r2 >= 0.95
        quality = '良好';
    elseif mean_r2 >= 0.9
        quality = '一般';
    elseif mean_r2 >= 0.8
        quality = '较差';
    else
        quality = '很差';
    end
    
    fprintf('拟合质量评估: %s (平均R² = %.4f)\n', quality, mean_r2);
    fprintf('建议：\n');
    
    if mean_r2 < 0.9
        fprintf('  - 考虑增加拟合阶数（当前阶数: %d）\n', rank);
        fprintf('  - 检查数据范围是否合理\n');
        fprintf('  - 考虑增加采样点数量\n');
    elseif mean_r2 < 0.95
        fprintf('  - 拟合质量可以接受，但可以尝试增加阶数以获得更好效果\n');
    else
        fprintf('  - 拟合质量良好，可以用于实际应用\n');
    end
    
    fprintf('==============\n\n');
    
    %% 返回计算结果（可选）
    results.MSE = MSE;
    results.MAE = MAE;
    results.MAX_AE = MAX_AE;
    results.R_squared = R_squared;
    results.mean_r2 = mean_r2;
    results.quality = quality;
    results.n_vars = n_vars;
    
    % 如果有输出参数，则返回结果结构体
    if nargout > 0
        varargout{1} = results;
    end
end