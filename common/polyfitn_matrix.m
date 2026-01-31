function [K_sym, K_sim_param, results] = polyfitn_matrix(X_data, K_data, rank, is_evaluate)
%   根据给定的腿长和K矩阵序列，对K矩阵的每一个元素进行多项式拟合
%   支持任意维度的K_data数据域
%   增加了数值稳定性和错误处理机制

%   输入参数:
%   X_data/K_data - 采样点，应是一个 （X, K）的张量, 其中K是因变量张量，X是自变量张量
%       格式要求: 
%       sample_idx为采样点数总数
%       K_data(sample_idx, ...)  - 第一维是样本索引，后面的维度可以是任意维度
%       例如K是个n行m列的二维矩阵，即有n*m个因变量元素需要拟合
%       对于二维矩阵则此格式为
%       K_data = zeros(TOTAL_N_SAMPLES, N_u, N_x)
%       计算时K_data(sample_idx, :, :) = calc_K
%       三维、四维以此类推 
%       X_data(sample_idx, ...) 
%       例如需要进行左右腿的相关拟合，每次采样都有两个独立变量，可理解为1行2列的行向量
%       对于该行向量，格式为
%       X_data = zeros(TOTAL_N_SAMPLES, n_vars);
%       计算时X_data(sample_idx, :) = [var1, var2, ..., var_n];
%
%   rank - 拟合阶数，阶数越高，拟合越精准，算力资源消耗也越高 
%   is_evaluate - true or false，决定是否需要输出统计特征
%
%   输出参数:
%   K_sym - 拟合出的表达式，符号数组，维度与K_data(2:end)相同
%   K_sim_param - 多项式的所有系数，元胞数组
%   results - 可选的评估结果结构体

    fprintf("\n=== 正在运行" + mfilename + ".m ===\n");
    
    % 记录开始时间
    fit_start_time = tic;

    %% 数据预处理和检查
    fprintf('正在进行数据预处理和检查...\n');
    
    % 检查输入数据维度
    if size(X_data, 1) ~= size(K_data, 1)
        error('X_data 和 K_data 的样本数量不一致');
    end
    
    % 检查输入数据中的NaN和Inf
    nan_inf_check_X = any(any(isnan(X_data) | isinf(X_data)));
    nan_inf_check_K = any(any(isnan(K_data(:)) | isinf(K_data(:))));
    
    if nan_inf_check_X
        warning('X_data 中包含 NaN 或 Inf，这可能导致拟合失败');
        % 移除包含NaN或Inf的样本
        valid_samples = all(~isnan(X_data) & ~isinf(X_data), 2);
        X_data = X_data(valid_samples, :);
        K_data = K_data(valid_samples, :);
        fprintf('已移除 %d 个无效样本\n', sum(~valid_samples));
    end
    
    if nan_inf_check_K
        warning('K_data 中包含 NaN 或 Inf');
    end
    
    % 检查数据范围是否合理
    fprintf('数据统计信息:\n');
    for i = 1:size(X_data, 2)
        fprintf('  自变量 %d: 范围 [%.4g, %.4g], 均值 %.4g, 标准差 %.4g\n', ...
            i, min(X_data(:, i)), max(X_data(:, i)), ...
            mean(X_data(:, i)), std(X_data(:, i)));
    end
    
    %% 提取采样点和K数据维度
    dims = size(K_data); % 拟合维度
    total_samples = dims(1); % 采样点数，第一维必须是样本数
    
    % 处理K_dims（注意：当K_data只有2维时，dims(2:end)可能返回标量）
    if length(dims) == 2
        % 如果K_data只有2维，说明每个样本的K是一个标量
        K_dims = 1; % 标量情况
        n_elements = 1;
    else
        K_dims = dims(2:end); % K矩阵/张量的维度
        n_elements = prod(K_dims); % 计算总的需要拟合的元素数
    end
    
    n_vars = size(X_data, 2); % 自变量数量
    
    % 检查样本数量是否足够
    min_samples_needed = (rank+1)^n_vars;
    if total_samples < min_samples_needed
        warning('样本数量 %d 可能不足，建议至少 %d 个样本以获得稳定的 %d 阶 %d 元多项式拟合', ...
            total_samples, min_samples_needed, rank, n_vars);
    end
    
    % 将K_data重塑为二维矩阵，便于处理
    % 第一维：样本数，第二维：所有K元素的展平
    if length(dims) == 2
        % 标量情况，K_data已经是二维的 [samples, 1]
        K_flat = K_data;
    else
        K_flat = reshape(K_data, total_samples, n_elements);
    end
    
    %% 初始化输出
    % 创建符号表达式数组（将展平后再reshape回原始形状）
    K_sym_flat = cell(1, n_elements);
    K_sim_param = cell(1, n_elements); % 使用元胞数组存储每个元素的系数
    
    % 记录拟合失败的元素的索引
    failed_indices = [];
    failed_reasons = {};
    
    %% 用于显示进度
    progress_step = 0;
    last_percent = 0;
    if n_elements > 0
        progress_interval = max(1, round(n_elements / 100));
    else
        progress_interval = 1;
    end

    %% 开始拟合
    fprintf('K拟合进度:')
    
    % 预分配进度条空间
    fprintf(repmat(' ', 1, 50));
    fprintf(repmat('\b', 1, 50));
    
    for elem_idx = 1:n_elements
        try
            % 提取当前K元素的所有样本值
            K = K_flat(:, elem_idx);
            
            % 检查当前K元素是否包含无效值
            valid_samples = ~isnan(K) & ~isinf(K);
            valid_K = K(valid_samples);
            valid_X = X_data(valid_samples, :);
            
            if length(valid_K) < min_samples_needed
                warning('元素 %d: 有效样本数量不足 (%d < %d)，尝试降低阶数或增加数据', ...
                    elem_idx, length(valid_K), min_samples_needed);
                
                % 尝试使用较低的阶数
                if n_vars == 1
                    max_possible_rank = length(valid_K) - 1;
                else
                    max_possible_rank = floor(length(valid_K)^(1/n_vars)) - 1;
                end
                
                if max_possible_rank > 0
                    actual_rank = min(rank, max_possible_rank);
                    fprintf('元素 %d: 使用降阶拟合 (阶数从 %d 降为 %d)\n', ...
                        elem_idx, rank, actual_rank);
                    polyval = polyfitn(valid_X, valid_K, actual_rank);
                else
                    error('有效样本太少，无法拟合');
                end
            else
                % 正常拟合
                polyval = polyfitn(valid_X, valid_K, rank);
            end
            
            % 存储多项式系数到元胞数组
            K_sim_param{elem_idx} = polyval;
            
            % 生成符号表达式
            K_sym_flat{elem_idx} = simplify(polyn2sym(polyval));
            
        catch ME
            % 记录失败的拟合
            failed_indices(end+1) = elem_idx;
            failed_reasons{end+1} = ME.message;
            
            % 尝试使用简单的最小二乘拟合作为降级方案
            try
                warning('元素 %d: polyfitn失败，尝试使用最小二乘法降级拟合', elem_idx);
                
                % 构建多项式设计矩阵
                X_design = build_polynomial_matrix(valid_X, min(rank, max_possible_rank));
                
                % 使用最小二乘法
                coeffs = (X_design' * X_design) \ (X_design' * valid_K);
                
                % 创建polyfitn兼容的结构
                polyval = struct();
                polyval.Coefficients = coeffs;
                % 生成多项式项
                if n_vars == 1
                    polyval.ModelTerms = (0:min(rank, max_possible_rank))';
                else
                    % 对于多变量，简化处理
                    polyval.ModelTerms = nchoosek(1:n_vars+min(rank, max_possible_rank), n_vars) - 1;
                    polyval.ModelTerms = unique(polyval.ModelTerms, 'rows');
                    polyval.ModelTerms = polyval.ModelTerms(sum(polyval.ModelTerms, 2) <= min(rank, max_possible_rank), :);
                end
                
                % 存储结果
                K_sim_param{elem_idx} = polyval;
                K_sym_flat{elem_idx} = polyn2sym(polyval);
                
                fprintf('元素 %d: 降级拟合成功\n', elem_idx);
                
            catch ME2
                % 如果降级拟合也失败，则使用常数拟合
                warning('元素 %d: 降级拟合也失败，使用常数拟合: %s', elem_idx, ME2.message);
                
                polyval = struct();
                polyval.Coefficients = mean(valid_K);
                if n_vars == 1
                    polyval.ModelTerms = 0;
                else
                    polyval.ModelTerms = zeros(1, n_vars);
                end
                
                K_sim_param{elem_idx} = polyval;
                K_sym_flat{elem_idx} = sym(mean(valid_K));
            end
        end

        % 更新进度显示
        progress_step = progress_step + 1;
        current_percent = floor(progress_step / n_elements * 100);
        if current_percent > last_percent
            fprintf(repmat('#', 1, min(5, current_percent - last_percent)));
            last_percent = current_percent;
        end
    end

    if last_percent < 100
        fprintf(repmat('#', 1, 100 - last_percent));
    end
    
    % 报告拟合失败的情况
    if ~isempty(failed_indices)
        fprintf('\n警告: %d/%d 个元素在拟合过程中遇到问题:\n', ...
            length(failed_indices), n_elements);
        for i = 1:min(5, length(failed_indices)) % 只显示前5个错误
            fprintf('  元素 %d: %s\n', failed_indices(i), failed_reasons{i});
        end
        if length(failed_indices) > 5
            fprintf('  ... 还有 %d 个错误未显示\n', length(failed_indices)-5);
        end
    end
    
    %% 将K_sym从元胞数组转换回符号数组，并reshape为原始K的形状
    try
        % 创建适当大小的符号数组
        if length(dims) == 2
            % 标量情况
            K_sym = K_sym_flat{1};
        elseif length(K_dims) == 1
            % 向量情况
            K_sym = sym(zeros(K_dims, 1));
            for i = 1:K_dims
                K_sym(i) = K_sym_flat{i};
            end
        else
            % 矩阵或高维数组情况
            K_sym = sym(zeros(K_dims));
            
            % 将展平的索引映射回多维索引
            indices = cell(1, length(K_dims));
            [indices{:}] = ind2sub(K_dims, 1:n_elements);
            
            for elem_idx = 1:n_elements
                % 获取多维索引
                idx_cell = cell(1, length(K_dims));
                for dim_idx = 1:length(K_dims)
                    idx_cell{dim_idx} = indices{dim_idx}(elem_idx);
                end
                
                % 赋值
                K_sym(idx_cell{:}) = K_sym_flat{elem_idx};
            end
        end
    catch ME
        warning('创建符号数组失败: %s', ME.message);
        % 返回元胞数组
        K_sym = K_sym_flat;
    end

    %% 计算耗时
    fit_time = toc(fit_start_time);

    %% 输出拟合特征和统计信息
    fprintf('\n');
    fprintf('拟合完成！\n');
    fprintf('  - 拟合阶数: %d阶\n', rank);
    fprintf('  - 自变量维度: %d\n', n_vars);
    fprintf('  - 采样点数量: %d\n', total_samples);
    fprintf('  - K数据域维度: ');
    if length(dims) == 2
        fprintf('标量\n');
    else
        fprintf('%d ', K_dims);
        fprintf('\n');
    end
    fprintf('  - 待拟合元素总数: %d\n', n_elements);
    fprintf('  - 成功拟合元素: %d\n', n_elements - length(failed_indices));
    fprintf('  - 失败拟合元素: %d\n', length(failed_indices));
    
    % 计算多项式项数
    if n_vars == 1
        n_terms = rank + 1;
        fprintf('  - 多项式项数: %d项 (一元多项式)\n', n_terms);
    elseif n_vars == 2
        n_terms = (rank+1)*(rank+2)/2;
        fprintf('  - 多项式项数: %d项 (二元多项式)\n', n_terms);
    else
        fprintf('  - 变量数: %d\n', n_vars);
    end
    
    fprintf('  - 拟合总耗时: %.2f秒 (约%.2f分钟)\n', fit_time, fit_time/60);
    if n_elements > 0
        fprintf('  - 平均每个元素耗时: %.3f毫秒\n', fit_time/n_elements*1000);
    end
    fprintf('\n');
    
    %% 评估拟合结果（可选）
    results = [];
    if (is_evaluate == true)
        results = evaluate_fit_results(X_data, K_data, K_sym_flat, K_sim_param, rank, dims, true);
    end
end

%% 辅助函数：构建多项式设计矩阵
function X_design = build_polynomial_matrix(X, degree)
% 构建多项式设计矩阵
% X: [n_samples × n_vars] 自变量矩阵
% degree: 多项式阶数
% X_design: [n_samples × n_terms] 设计矩阵

    n_samples = size(X, 1);
    n_vars = size(X, 2);
    
    if n_vars == 1
        % 单变量情况
        X_design = zeros(n_samples, degree+1);
        for d = 0:degree
            X_design(:, d+1) = X.^d;
        end
    else
        % 多变量情况（简化处理）
        % 生成所有可能的指数组合
        exponents = {};
        for d = 0:degree
            combos = nchoosek(1:n_vars+d, n_vars) - 1;
            combos = unique(combos, 'rows');
            combos = combos(sum(combos, 2) == d, :);
            exponents = [exponents; num2cell(combos, 2)];
        end
        
        n_terms = length(exponents);
        X_design = ones(n_samples, n_terms);
        
        % 构建设计矩阵
        for term_idx = 1:n_terms
            for var_idx = 1:n_vars
                if exponents{term_idx}(var_idx) > 0
                    X_design(:, term_idx) = X_design(:, term_idx) .* ...
                        (X(:, var_idx) .^ exponents{term_idx}(var_idx));
                end
            end
        end
    end
end

%% 辅助函数：评估拟合结果（优化版）
function results = evaluate_fit_results(X_data, K_data, K_sym_flat, K_sim_param, rank, dims, show_plots)
%   评估多项式拟合结果，支持任意维度的K_data
%
%   输入参数:
%   X_data - 自变量数据，格式为 [total_samples × n_vars]
%   K_data - 原始K矩阵数据，格式为 [total_samples × K_dims...]
%   K_sym_flat - 拟合出的符号表达式元胞数组
%   K_sim_param - 多项式拟合参数元胞数组
%   rank - 拟合阶数
%   dims - K_data的原始维度
%   show_plots - 是否显示绘图，默认为true

    if nargin < 7
        show_plots = true;
    end
    
    fprintf("\n=== 正在评估拟合结果 ===\n");
    
    %% 基本参数和输入检查
    total_samples = dims(1);
    
    % 处理K_dims
    if length(dims) == 2
        % 标量情况
        K_dims = 1;
        n_elements = 1;
    else
        K_dims = dims(2:end);
        n_elements = prod(K_dims);
    end
    
    n_vars = size(X_data, 2);
    
    % 检查输入数据维度
    if total_samples ~= size(X_data, 1)
        error('X_data 和 K_data 的样本数量不一致');
    end
    
    % 检查K_sym_flat的维度
    if length(K_sym_flat) ~= n_elements
        error('K_sym_flat 的长度与K数据域元素数量不一致');
    end
    
    fprintf('数据维度: %d 个样本, %d 个自变量\n', total_samples, n_vars);
    fprintf('K数据域维度: ');
    if length(dims) == 2
        fprintf('标量\n');
    else
        fprintf('%d ', K_dims);
        fprintf('\n');
    end
    fprintf('拟合阶数: %d\n', rank);
    
    %% 计算预测值和误差
    % 将K_data展平为二维矩阵以便处理
    if length(dims) == 2
        K_flat = K_data;
    else
        K_flat = reshape(K_data, total_samples, n_elements);
    end
    
    % 预分配误差矩阵
    MSE = zeros(1, n_elements);  % 均方误差
    MAE = zeros(1, n_elements);  % 平均绝对误差
    MAX_AE = zeros(1, n_elements);  % 最大绝对误差
    R_squared = zeros(1, n_elements);  % 决定系数R²
    
    fprintf('计算拟合误差...\n');
    
    % 进度显示
    fprintf('进度: [');
    if n_elements > 0
        progress_interval = max(1, round(n_elements / 50));
    else
        progress_interval = 1;
    end
    
    for elem_idx = 1:n_elements
        try
            % 获取原始数据
            K_actual = K_flat(:, elem_idx);
            
            % 检查数据有效性
            valid_samples = ~isnan(K_actual) & ~isinf(K_actual);
            if sum(valid_samples) == 0
                MSE(elem_idx) = NaN;
                MAE(elem_idx) = NaN;
                MAX_AE(elem_idx) = NaN;
                R_squared(elem_idx) = NaN;
                continue;
            end
            
            K_actual_valid = K_actual(valid_samples);
            X_data_valid = X_data(valid_samples, :);
            
            % 使用多项式参数计算预测值
            if ~isempty(K_sim_param{elem_idx}) && isfield(K_sim_param{elem_idx}, 'Coefficients')
                % 使用polyvaln函数（如果可用）或自定义计算
                K_pred = compute_polyval(K_sim_param{elem_idx}, X_data_valid);
            else
                % 如果多项式参数不可用，尝试使用符号表达式
                K_pred = zeros(size(K_actual_valid));
                for k = 1:length(K_actual_valid)
                    try
                        % 创建符号变量
                        x_sym = sym('x', [1, n_vars]);
                        subs_list = num2cell(X_data_valid(k, :));
                        K_pred(k) = double(subs(K_sym_flat{elem_idx}, x_sym, subs_list));
                    catch
                        K_pred(k) = NaN;
                    end
                end
            end
            
            % 检查预测值的有效性
            valid_pred = ~isnan(K_pred) & ~isinf(K_pred);
            if sum(valid_pred) < 2  % 至少需要2个有效点才能计算R²
                MSE(elem_idx) = NaN;
                MAE(elem_idx) = NaN;
                MAX_AE(elem_idx) = NaN;
                R_squared(elem_idx) = NaN;
                continue;
            end
            
            K_actual_valid = K_actual_valid(valid_pred);
            K_pred_valid = K_pred(valid_pred);
            
            % 计算误差指标
            residuals = K_actual_valid - K_pred_valid;
            
            % 均方误差
            MSE(elem_idx) = mean(residuals.^2);
            
            % 平均绝对误差
            MAE(elem_idx) = mean(abs(residuals));
            
            % 最大绝对误差
            MAX_AE(elem_idx) = max(abs(residuals));
            
            % 决定系数R²
            SS_res = sum(residuals.^2);
            SS_tot = sum((K_actual_valid - mean(K_actual_valid)).^2);
            if SS_tot > 0
                R_squared(elem_idx) = 1 - SS_res / SS_tot;
            else
                R_squared(elem_idx) = 1;
            end
            
        catch ME
            % 如果计算过程中出现错误，记录NaN
            fprintf('\n元素 %d 评估失败: %s\n', elem_idx, ME.message);
            MSE(elem_idx) = NaN;
            MAE(elem_idx) = NaN;
            MAX_AE(elem_idx) = NaN;
            R_squared(elem_idx) = NaN;
        end
        
        % 显示进度
        if mod(elem_idx, progress_interval) == 0 || elem_idx == n_elements
            fprintf('=');
        end
    end
    
    fprintf(']\n');
    
    % 移除NaN值进行统计
    MSE_valid = MSE(~isnan(MSE));
    MAE_valid = MAE(~isnan(MAE));
    MAX_AE_valid = MAX_AE(~isnan(MAX_AE));
    R_squared_valid = R_squared(~isnan(R_squared));
    
    %% 将误差指标reshape回原始K的形状
    % 注意：当K_dims是标量时，不需要reshape
    if length(dims) == 2
        % 标量情况
        MSE_reshaped = MSE;
        MAE_reshaped = MAE;
        MAX_AE_reshaped = MAX_AE;
        R_squared_reshaped = R_squared;
    elseif length(K_dims) == 1
        % 向量情况
        MSE_reshaped = reshape(MSE, [K_dims, 1]);
        MAE_reshaped = reshape(MAE, [K_dims, 1]);
        MAX_AE_reshaped = reshape(MAX_AE, [K_dims, 1]);
        R_squared_reshaped = reshape(R_squared, [K_dims, 1]);
    else
        % 矩阵或高维数组情况
        MSE_reshaped = reshape(MSE, K_dims);
        MAE_reshaped = reshape(MAE, K_dims);
        MAX_AE_reshaped = reshape(MAX_AE, K_dims);
        R_squared_reshaped = reshape(R_squared, K_dims);
    end

    %% 输出整体统计
    fprintf('拟合误差统计：\n');
    fprintf('==============\n');
    
    if ~isempty(MSE_valid)
        fprintf('整体均方误差 (MSE):\n');
        fprintf('  平均值: %.4e\n', mean(MSE_valid));
        fprintf('  最小值: %.4e\n', min(MSE_valid));
        fprintf('  最大值: %.4e\n', max(MSE_valid));
        fprintf('  中位数: %.4e\n', median(MSE_valid));
        
        % 显示MSE最差和最好的5个元素（如果元素数量足够）
        if length(MSE_valid) >= 5
            [sorted_MSE, sorted_MSE_idx] = sort(MSE_valid);
            fprintf('  最小的5个MSE: ');
            fprintf('%.4e ', sorted_MSE(1:min(5, length(sorted_MSE))));
            fprintf('\n');
            fprintf('  最大的5个MSE: ');
            fprintf('%.4e ', sorted_MSE(max(1, end-4):end));
            fprintf('\n');
        end
    else
        fprintf('整体均方误差 (MSE):\n');
        fprintf('  所有MSE值均为NaN\n');
    end
    
    if ~isempty(MAE_valid)
        fprintf('\n整体平均绝对误差 (MAE):\n');
        fprintf('  平均值: %.4e\n', mean(MAE_valid));
        fprintf('  最小值: %.4e\n', min(MAE_valid));
        fprintf('  最大值: %.4e\n', max(MAE_valid));
        fprintf('  中位数: %.4e\n', median(MAE_valid));
        
        % 显示MAE分布的分位数
        if length(MAE_valid) >= 10
            quantiles = quantile(MAE_valid, [0.25, 0.5, 0.75, 0.9, 0.95]);
            fprintf('  分布分位数:\n');
            fprintf('    25%%分位数: %.4e\n', quantiles(1));
            fprintf('    中位数: %.4e\n', quantiles(2));
            fprintf('    75%%分位数: %.4e\n', quantiles(3));
            fprintf('    90%%分位数: %.4e\n', quantiles(4));
            fprintf('    95%%分位数: %.4e\n', quantiles(5));
        end
    else
        fprintf('\n整体平均绝对误差 (MAE):\n');
        fprintf('  所有MAE值均为NaN\n');
    end
    
    if ~isempty(MAX_AE_valid)
        fprintf('\n整体最大绝对误差:\n');
        fprintf('  平均值: %.4e\n', mean(MAX_AE_valid));
        fprintf('  最小值: %.4e\n', min(MAX_AE_valid));
        fprintf('  最大值: %.4e\n', max(MAX_AE_valid));
        
        % 统计超过阈值的元素数量
        if ~isempty(MAX_AE_valid)
            thresholds = [1e-3, 1e-2, 1e-1, 1, 10];
            fprintf('  超过阈值的元素数量:\n');
            for thresh = thresholds
                count = sum(MAX_AE_valid > thresh);
                if count > 0
                    fprintf('    > %.1e: %d/%d (%.1f%%)\n', ...
                        thresh, count, length(MAX_AE_valid), count/length(MAX_AE_valid)*100);
                end
            end
        end
    else
        fprintf('\n整体最大绝对误差:\n');
        fprintf('  所有MAX_AE值均为NaN\n');
    end
    
    if ~isempty(R_squared_valid)
        fprintf('\n决定系数 R²:\n');
        fprintf('  平均值: %.4f\n', mean(R_squared_valid));
        fprintf('  最小值: %.4f\n', min(R_squared_valid));
        fprintf('  最大值: %.4f\n', max(R_squared_valid));
        fprintf('  中位数: %.4f\n', median(R_squared_valid));
        
        % R²分布统计
        r2_ranges = [0, 0.5, 0.8, 0.9, 0.95, 0.99, 1.0];
        fprintf('  R²分布:\n');
        for i = 1:length(r2_ranges)-1
            lower_bound = r2_ranges(i);
            upper_bound = r2_ranges(i+1);
            if i == length(r2_ranges)-1
                % 最后一个区间是闭区间
                count = sum(R_squared_valid >= lower_bound & R_squared_valid <= upper_bound);
            else
                count = sum(R_squared_valid >= lower_bound & R_squared_valid < upper_bound);
            end
            
            if i == length(r2_ranges)-1
                range_str = sprintf('[%.2f, %.2f]', lower_bound, upper_bound);
            else
                range_str = sprintf('[%.2f, %.2f)', lower_bound, upper_bound);
            end
            
            percentage = count/length(R_squared_valid)*100;
            fprintf('    %s: %d/%d (%.1f%%)\n', ...
                range_str, count, length(R_squared_valid), percentage);
        end
        
        % 常见阈值统计
        fprintf('  R² < 0.5 的元素数: %d/%d (%.1f%%)\n', ...
            sum(R_squared_valid < 0.5), length(R_squared_valid), ...
            sum(R_squared_valid < 0.5)/length(R_squared_valid)*100);
        fprintf('  R² < 0.8 的元素数: %d/%d (%.1f%%)\n', ...
            sum(R_squared_valid < 0.8), length(R_squared_valid), ...
            sum(R_squared_valid < 0.8)/length(R_squared_valid)*100);
        fprintf('  R² < 0.9 的元素数: %d/%d (%.1f%%)\n', ...
            sum(R_squared_valid < 0.9), length(R_squared_valid), ...
            sum(R_squared_valid < 0.9)/length(R_squared_valid)*100);
        fprintf('  R² < 0.95 的元素数: %d/%d (%.1f%%)\n', ...
            sum(R_squared_valid < 0.95), length(R_squared_valid), ...
            sum(R_squared_valid < 0.95)/length(R_squared_valid)*100);
        fprintf('  R² >= 0.99 的元素数: %d/%d (%.1f%%)\n', ...
            sum(R_squared_valid >= 0.99), length(R_squared_valid), ...
            sum(R_squared_valid >= 0.99)/length(R_squared_valid)*100);
        
        % 显示R²最好和最差的元素
        if length(R_squared_valid) >= 5
            [sorted_R2, sorted_R2_idx] = sort(R_squared_valid);
            fprintf('  最差的5个R²: ');
            fprintf('%.4f ', sorted_R2(1:min(5, length(sorted_R2))));
            fprintf('\n');
            fprintf('  最好的5个R²: ');
            fprintf('%.4f ', sorted_R2(max(1, end-4):end));
            fprintf('\n');
        end
        
        % 计算拟合优度的其他指标
        mean_abs_error = mean(MAE_valid);
        std_abs_error = std(MAE_valid);
        fprintf('  平均绝对误差均值: %.4e\n', mean_abs_error);
        fprintf('  平均绝对误差标准差: %.4e\n', std_abs_error);
        fprintf('  变异系数(CV): %.2f%%\n', std_abs_error/mean_abs_error*100);
        
    else
        fprintf('\n决定系数 R²:\n');
        fprintf('  所有R²值均为NaN\n');
    end
    
    % 数据有效性统计
    valid_count = length(R_squared_valid);
    total_count = n_elements;
    invalid_count = total_count - valid_count;
    
    fprintf('\n数据有效性统计:\n');
    fprintf('  有效元素数量: %d/%d (%.1f%%)\n', ...
        valid_count, total_count, valid_count/total_count*100);
    fprintf('  无效元素数量(NaN): %d/%d (%.1f%%)\n', ...
        invalid_count, total_count, invalid_count/total_count*100);
    
    % 如果存在无效元素，显示原因
    if invalid_count > 0
        nan_indices = find(isnan(R_squared));
        fprintf('  无效元素索引(前10个): ');
        if length(nan_indices) > 10
            fprintf('%d ', nan_indices(1:10));
            fprintf('... (还有 %d 个)\n', length(nan_indices)-10);
        else
            fprintf('%d ', nan_indices);
            fprintf('\n');
        end
    end
    
    % 计算拟合优度的综合评分
    if ~isempty(R_squared_valid)
        % 基于R²的评分系统
        mean_r2 = mean(R_squared_valid);
        median_r2 = median(R_squared_valid);
        min_r2 = min(R_squared_valid);
        
        % 基于多个指标的综合评分（0-100分）
        r2_score = mean_r2 * 100;  % R²贡献（最高100分）
        consistency_score = min(100, 100 * (median_r2 / mean_r2));  % 一致性贡献
        robustness_score = min(100, 100 * (1 - sum(isnan(R_squared))/total_count));  % 鲁棒性贡献
        
        % 加权综合评分
        overall_score = 0.6 * r2_score + 0.2 * consistency_score + 0.2 * robustness_score;
        
        fprintf('\n拟合优度综合评分:\n');
        fprintf('  平均R²得分: %.1f/100\n', r2_score);
        fprintf('  一致性得分: %.1f/100 (中位数/平均值: %.3f)\n', consistency_score, median_r2/mean_r2);
        fprintf('  鲁棒性得分: %.1f/100 (有效拟合比例: %.1f%%)\n', robustness_score, valid_count/total_count*100);
        fprintf('  综合得分: %.1f/100\n', overall_score);
        
        % 根据综合评分给出评级
        if overall_score >= 90
            rating = '★★★★★ (优秀)';
        elseif overall_score >= 80
            rating = '★★★★ (良好)';
        elseif overall_score >= 70
            rating = '★★★ (中等)';
        elseif overall_score >= 60
            rating = '★★ (合格)';
        else
            rating = '★ (需要改进)';
        end
        fprintf('  评级: %s\n', rating);
    else
        fprintf('\n拟合优度综合评分:\n');
        fprintf('  无法计算综合评分（无有效数据）\n');
    end
    
    fprintf('==============\n\n');
    
    %% 绘制误差分布图
    if show_plots && ~isempty(R_squared_valid) && length(R_squared_valid) > 1
        try
            fprintf('生成误差分析图...\n');
            
            % 创建主要图形窗口
            main_fig = figure('Name', '多项式拟合误差分析', ...
                'Position', [100, 100, 1200, 600], ...
                'Color', 'w');
            
            % 根据K的维度决定如何显示
            if length(dims) == 2 || length(K_dims) == 1
                % 标量或向量情况，使用简化绘图
                plot_scalar_vector_results(X_data, K_data, K_sym_flat, K_sim_param, ...
                    MSE_valid, MAE_valid, R_squared_valid, rank, n_elements, n_vars, dims);
            elseif length(K_dims) == 2
                % 二维矩阵情况，使用详细绘图
                plot_matrix_results(X_data, K_data, K_sym_flat, K_sim_param, ...
                    MSE, MAE, R_squared, MSE_reshaped, MAE_reshaped, R_squared_reshaped, ...
                    rank, n_elements, n_vars, K_dims);
            else
                % 高维张量情况，使用简化的三维可视化
                plot_high_dimensional_results(X_data, K_flat, K_sym_flat, K_sim_param, ...
                    MSE, MAE, R_squared, rank, n_elements, n_vars, K_dims);
            end
            
            fprintf('误差分析图已生成。\n');
            
            % 保存图形（可选）
            try
                save_path = fullfile(pwd, 'fitting_results');
                if ~exist(save_path, 'dir')
                    mkdir(save_path);
                end
                
                timestamp = datestr(now, 'yyyymmdd_HHMMSS');
                filename = sprintf('fitting_analysis_%s.png', timestamp);
                saveas(main_fig, fullfile(save_path, filename));
                fprintf('图形已保存至: %s\n', fullfile(save_path, filename));
            catch
                % 如果保存失败，忽略错误
            end
            
        catch ME
            fprintf('绘图过程中出现错误: %s\n', ME.message);
            fprintf('跳过绘图步骤...\n');
        end
    elseif show_plots
        if isempty(R_squared_valid)
            fprintf('无有效数据，无法生成误差分析图。\n');
        elseif length(R_squared_valid) <= 1
            fprintf('有效数据点不足，无法生成完整的误差分析图。\n');
        end
    end
    
    %% 返回计算结果
    results.MSE = MSE_reshaped;
    results.MAE = MAE_reshaped;
    results.MAX_AE = MAX_AE_reshaped;
    results.R_squared = R_squared_reshaped;
    
    if ~isempty(R_squared_valid)
        results.mean_r2 = mean(R_squared_valid);
        % 根据R²平均值评估拟合质量
        mean_r2 = results.mean_r2;
        if mean_r2 >= 0.99
            results.quality = '优秀';
        elseif mean_r2 >= 0.95
            results.quality = '良好';
        elseif mean_r2 >= 0.9
            results.quality = '一般';
        elseif mean_r2 >= 0.8
            results.quality = '较差';
        else
            results.quality = '很差';
        end
    else
        results.mean_r2 = NaN;
        results.quality = '未知';
    end
    
    results.n_vars = n_vars;
    if length(dims) == 2
        results.K_dims = 1;
    else
        results.K_dims = K_dims;
    end
    results.n_elements = n_elements;
end

%% 辅助函数：计算多项式值
function y_pred = compute_polyval(poly_struct, X)
% 使用多项式结构计算预测值
    
    if isempty(poly_struct) || ~isfield(poly_struct, 'Coefficients')
        y_pred = zeros(size(X, 1), 1);
        return;
    end
    
    try
        % 尝试使用polyvaln（如果可用）
        if exist('polyvaln', 'file') == 2
            y_pred = polyvaln(poly_struct, X);
        else
            % 自定义计算
            n_samples = size(X, 1);
            n_terms = length(poly_struct.Coefficients);
            y_pred = zeros(n_samples, 1);
            
            for i = 1:n_samples
                term_value = 0;
                for j = 1:n_terms
                    term = 1;
                    for k = 1:size(X, 2)
                        if poly_struct.ModelTerms(j, k) > 0
                            term = term * (X(i, k) ^ poly_struct.ModelTerms(j, k));
                        end
                    end
                    term_value = term_value + poly_struct.Coefficients(j) * term;
                end
                y_pred(i) = term_value;
            end
        end
    catch
        % 如果计算失败，返回零向量
        y_pred = zeros(size(X, 1), 1);
    end
end

%% 子函数：绘制标量/向量结果的图形
function plot_scalar_vector_results(X_data, K_data, K_sym_flat, K_sim_param, ...
    MSE_valid, MAE_valid, R_squared_valid, rank, n_elements, n_vars, dims)
    
    % 创建4x3的子图布局
    subplot(4, 3, 1);
    histogram(R_squared_valid, 20, 'FaceColor', [0.2, 0.4, 0.8], 'EdgeColor', 'k');
    xlabel('决定系数 R²');
    ylabel('频数');
    title('R² 分布直方图', 'FontWeight', 'bold');
    grid on;
    hold on;
    xline(0.9, 'r--', 'R²=0.9', 'LineWidth', 1.5);
    xline(0.95, 'g--', 'R²=0.95', 'LineWidth', 1.5);
    hold off;
    
    subplot(4, 3, 2);
    histogram(MAE_valid, 20, 'FaceColor', [0.8, 0.4, 0.2], 'EdgeColor', 'k');
    xlabel('平均绝对误差 (MAE)');
    ylabel('频数');
    title('MAE 分布直方图', 'FontWeight', 'bold');
    grid on;
    
    subplot(4, 3, 3);
    histogram(log10(MSE_valid), 20, 'FaceColor', [0.4, 0.8, 0.2], 'EdgeColor', 'k');
    xlabel('log_{10}(MSE)');
    ylabel('频数');
    title('MSE 对数分布直方图', 'FontWeight', 'bold');
    grid on;
    
    % 选择R²最高的元素作为示例
    [max_r2, max_idx] = max(R_squared_valid);
    
    % 获取原始K_data中对应的元素
    if length(dims) == 2
        K_actual_best = K_data;
    else
        K_flat = reshape(K_data, size(K_data, 1), []);
        K_actual_best = K_flat(:, max_idx);
    end
    
    % 计算预测值
    X_data_valid = X_data(~isnan(K_actual_best) & ~isinf(K_actual_best), :);
    K_actual_valid = K_actual_best(~isnan(K_actual_best) & ~isinf(K_actual_best));
    
    if ~isempty(K_sim_param{max_idx})
        K_pred_best = compute_polyval(K_sim_param{max_idx}, X_data_valid);
    else
        K_pred_best = zeros(size(K_actual_valid));
    end
    
    % 子图4: 拟合效果散点图
    subplot(4, 3, [4, 5, 6]);
    scatter(K_actual_valid, K_pred_best, 40, 'filled', ...
        'MarkerFaceColor', [0.2, 0.6, 0.8], 'MarkerEdgeColor', 'k');
    hold on;
    % 绘制y=x参考线
    min_val = min([K_actual_valid; K_pred_best]);
    max_val = max([K_actual_valid; K_pred_best]);
    plot([min_val, max_val], [min_val, max_val], 'r--', 'LineWidth', 2);
    
    % 添加回归线
    if length(K_actual_valid) > 1
        p = polyfit(K_actual_valid, K_pred_best, 1);
        x_fit = linspace(min_val, max_val, 100);
        y_fit = polyval(p, x_fit);
        plot(x_fit, y_fit, 'g-', 'LineWidth', 2);
        
        % 计算回归线的R²
        y_pred_reg = polyval(p, K_actual_valid);
        SS_res_reg = sum((K_pred_best - y_pred_reg).^2);
        SS_tot_reg = sum((K_pred_best - mean(K_pred_best)).^2);
        r2_reg = 1 - SS_res_reg/SS_tot_reg;
        
        legend({'数据点', 'y=x参考线', sprintf('回归线 (R²=%.4f)', r2_reg)}, ...
            'Location', 'best', 'FontSize', 9);
    else
        legend({'数据点', 'y=x参考线'}, 'Location', 'best', 'FontSize', 9);
    end
    hold off;
    
    xlabel('实际值', 'FontSize', 11);
    ylabel('预测值', 'FontSize', 11);
    title(sprintf('拟合效果散点图 (最佳元素, R²=%.4f)', max_r2), 'FontWeight', 'bold');
    grid on;
    axis equal tight;
    
    % 子图7-9: 残差分析
    subplot(4, 3, 7);
    if length(K_actual_valid) > 0
        residuals = K_actual_valid - K_pred_best;
        histogram(residuals, 20, 'FaceColor', [0.8, 0.2, 0.2], 'EdgeColor', 'k');
        xlabel('残差');
        ylabel('频数');
        title('残差分布直方图', 'FontWeight', 'bold');
        grid on;
    else
        text(0.5, 0.5, '无有效数据', 'HorizontalAlignment', 'center');
        title('残差分布直方图', 'FontWeight', 'bold');
        axis off;
    end
    
    subplot(4, 3, 8);
    if length(K_actual_valid) > 0
        plot(K_pred_best, residuals, 'o', 'MarkerFaceColor', [0.8, 0.2, 0.2], ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 6);
        hold on;
        xlims = xlim;
        plot(xlims, [0, 0], 'k--', 'LineWidth', 1.5);
        hold off;
        xlabel('预测值');
        ylabel('残差');
        title('残差 vs 预测值', 'FontWeight', 'bold');
        grid on;
    else
        text(0.5, 0.5, '无有效数据', 'HorizontalAlignment', 'center');
        title('残差 vs 预测值', 'FontWeight', 'bold');
        axis off;
    end
    
    subplot(4, 3, 9);
    if length(K_actual_valid) > 0
        normplot(residuals);
        title('残差正态概率图', 'FontWeight', 'bold');
        grid on;
    else
        text(0.5, 0.5, '无有效数据', 'HorizontalAlignment', 'center');
        title('残差正态概率图', 'FontWeight', 'bold');
        axis off;
    end
    
    % 子图10-12: 根据自变量数量选择可视化方式
    subplot_positions = {[0.05, 0.05, 0.3, 0.2], [0.36, 0.05, 0.3, 0.2], [0.67, 0.05, 0.3, 0.2]};
    
    for i = 1:3
        subplot('Position', subplot_positions{i});
        
        if i == 1
            % 箱线图：误差统计
            if ~isempty(MAE_valid)
                box_data = [MAE_valid(:), MSE_valid(:)];
                boxplot(box_data, 'Labels', {'MAE', 'MSE'});
                ylabel('误差值');
                title('误差统计箱线图', 'FontWeight', 'bold');
                grid on;
            else
                text(0.5, 0.5, '无有效数据', 'HorizontalAlignment', 'center');
                title('误差统计箱线图', 'FontWeight', 'bold');
                axis off;
            end
            
        elseif i == 2 && n_vars == 1
            % 单变量拟合曲线
            try
                x_range = linspace(min(X_data(:,1)), max(X_data(:,1)), 100);
                K_curve = zeros(size(x_range));
                
                for k = 1:length(x_range)
                    K_curve(k) = eval_polynomial(K_sim_param{max_idx}, [x_range(k)]);
                end
                
                plot(x_range, K_curve, 'b-', 'LineWidth', 2);
                hold on;
                scatter(X_data_valid(:,1), K_actual_valid, 50, 'r', 'filled');
                hold off;
                
                xlabel('自变量 x1');
                ylabel('K值');
                title('单变量拟合曲线', 'FontWeight', 'bold');
                grid on;
                legend('拟合曲线', '原始数据', 'Location', 'best', 'FontSize', 8);
            catch
                text(0.5, 0.5, '无法绘制拟合曲线', 'HorizontalAlignment', 'center');
                title('单变量拟合曲线', 'FontWeight', 'bold');
                axis off;
            end
            
        elseif i == 2 && n_vars == 2
            % 双变量拟合曲面
            try
                x1_range = linspace(min(X_data(:,1)), max(X_data(:,1)), 20);
                x2_range = linspace(min(X_data(:,2)), max(X_data(:,2)), 20);
                [X1_grid, X2_grid] = meshgrid(x1_range, x2_range);
                
                % 计算拟合曲面
                K_surf = zeros(size(X1_grid));
                for m = 1:size(X1_grid, 1)
                    for n = 1:size(X1_grid, 2)
                        K_surf(m, n) = eval_polynomial(K_sim_param{max_idx}, ...
                            [X1_grid(m, n), X2_grid(m, n)]);
                    end
                end
                
                surf(X1_grid, X2_grid, K_surf, 'FaceAlpha', 0.8);
                hold on;
                scatter3(X_data_valid(:,1), X_data_valid(:,2), K_actual_valid, ...
                    50, 'r', 'filled');
                hold off;
                
                xlabel('自变量 x1');
                ylabel('自变量 x2');
                zlabel('K值');
                title('双变量拟合曲面', 'FontWeight', 'bold');
                grid on;
            catch
                text(0.5, 0.5, '无法绘制拟合曲面', 'HorizontalAlignment', 'center');
                title('双变量拟合曲面', 'FontWeight', 'bold');
                axis off;
            end
            
        elseif i == 2
            % 多变量情况
            text(0.5, 0.5, sprintf('自变量数: %d\n无法直接可视化', n_vars), ...
                'FontSize', 10, 'HorizontalAlignment', 'center');
            title('多变量拟合结果', 'FontWeight', 'bold');
            axis off;
            
        elseif i == 3
            % 拟合质量总结
            valid_count = length(R_squared_valid);
            total_count = n_elements;
            
            text_str = {
                sprintf('拟合质量总结');
                sprintf('==============');
                sprintf('拟合阶数: %d', rank);
                sprintf('自变量数: %d', n_vars);
                sprintf('元素总数: %d', total_count);
                sprintf('有效拟合: %d (%.1f%%)', valid_count, valid_count/total_count*100);
                sprintf('平均 R²: %.4f', mean(R_squared_valid));
                sprintf('最佳 R²: %.4f', max_r2);
                sprintf('最差 R²: %.4f', min(R_squared_valid));
                sprintf('==============');
            };
            
            text(0.05, 0.95, text_str, 'VerticalAlignment', 'top', ...
                'FontSize', 9, 'FontName', 'Monospaced');
            axis off;
        end
    end
    
    % 添加总标题
    sgtitle_str = sprintf('多项式拟合误差分析 (阶数: %d, 元素数: %d, 自变量: %d)', ...
        rank, n_elements, n_vars);
    sgtitle(sgtitle_str, 'FontSize', 14, 'FontWeight', 'bold');
end

%% 子函数：绘制矩阵结果的图形
function plot_matrix_results(X_data, K_data, K_sym_flat, K_sim_param, ...
    MSE, MAE, R_squared, MSE_reshaped, MAE_reshaped, R_squared_reshaped, ...
    rank, n_elements, n_vars, K_dims)
    
    % 创建3x3的子图布局
    subplot(3, 3, 1);
    histogram(R_squared, 20, 'FaceColor', [0.2, 0.4, 0.8], 'EdgeColor', 'k');
    xlabel('决定系数 R²');
    ylabel('频数');
    title('R² 分布直方图', 'FontWeight', 'bold');
    grid on;
    hold on;
    xline(0.9, 'r--', 'R²=0.9', 'LineWidth', 1.5);
    xline(0.95, 'g--', 'R²=0.95', 'LineWidth', 1.5);
    hold off;
    
    subplot(3, 3, 2);
    histogram(MAE, 20, 'FaceColor', [0.8, 0.4, 0.2], 'EdgeColor', 'k');
    xlabel('平均绝对误差 (MAE)');
    ylabel('频数');
    title('MAE 分布直方图', 'FontWeight', 'bold');
    grid on;
    
    subplot(3, 3, 3);
    histogram(log10(MSE), 20, 'FaceColor', [0.4, 0.8, 0.2], 'EdgeColor', 'k');
    xlabel('log_{10}(MSE)');
    ylabel('频数');
    title('MSE 对数分布直方图', 'FontWeight', 'bold');
    grid on;
    
    % 子图4: R²热力图
    subplot(3, 3, 4);
    imagesc(R_squared_reshaped);
    colorbar;
    colormap(jet);
    xlabel('列索引');
    ylabel('行索引');
    title('R² 热力图', 'FontWeight', 'bold');
    axis equal tight;
    
    % 添加数值标签（如果矩阵不太大）
    if K_dims(1) <= 20 && K_dims(2) <= 20
        [X, Y] = meshgrid(1:K_dims(2), 1:K_dims(1));
        for i = 1:numel(R_squared_reshaped)
            if ~isnan(R_squared_reshaped(i))
                text(X(i), Y(i), sprintf('%.2f', R_squared_reshaped(i)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 7, ...
                    'Color', 'w', 'FontWeight', 'bold');
            end
        end
    end
    
    % 子图5: MAE热力图
    subplot(3, 3, 5);
    imagesc(MAE_reshaped);
    colorbar;
    colormap(jet);
    xlabel('列索引');
    ylabel('行索引');
    title('MAE 热力图', 'FontWeight', 'bold');
    axis equal tight;
    
    % 子图6: MSE热力图
    subplot(3, 3, 6);
    mse_log = log10(MSE_reshaped);
    mse_log(isinf(mse_log)) = min(mse_log(~isinf(mse_log)));
    imagesc(mse_log);
    colorbar;
    colormap(jet);
    xlabel('列索引');
    ylabel('行索引');
    title('log_{10}(MSE) 热力图', 'FontWeight', 'bold');
    axis equal tight;
    
    % 选择R²最高的元素作为示例
    [max_r2, max_idx] = max(R_squared);
    [i_best, j_best] = ind2sub(K_dims, max_idx);
    
    % 获取原始K_data中对应的元素
    K_flat = reshape(K_data, size(K_data, 1), []);
    K_actual_best = K_flat(:, max_idx);
    
    % 计算预测值
    X_data_valid = X_data(~isnan(K_actual_best) & ~isinf(K_actual_best), :);
    K_actual_valid = K_actual_best(~isnan(K_actual_best) & ~isinf(K_actual_best));
    
    if ~isempty(K_sim_param{max_idx})
        K_pred_best = compute_polyval(K_sim_param{max_idx}, X_data_valid);
    else
        K_pred_best = zeros(size(K_actual_valid));
    end
    
    % 子图7: 最佳元素拟合效果
    subplot(3, 3, 7);
    scatter(K_actual_valid, K_pred_best, 30, 'filled', ...
        'MarkerFaceColor', [0.2, 0.6, 0.8], 'MarkerEdgeColor', 'k');
    hold on;
    min_val = min([K_actual_valid; K_pred_best]);
    max_val = max([K_actual_valid; K_pred_best]);
    plot([min_val, max_val], [min_val, max_val], 'r--', 'LineWidth', 2);
    hold off;
    
    xlabel('实际值');
    ylabel('预测值');
    title(sprintf('最佳元素 [%d,%d] (R²=%.4f)', i_best, j_best, max_r2), ...
        'FontWeight', 'bold');
    grid on;
    axis equal tight;
    
    % 子图8: 根据自变量数量选择可视化方式
    subplot(3, 3, 8);
    if n_vars == 1
        % 单变量拟合曲线
        try
            x_range = linspace(min(X_data(:,1)), max(X_data(:,1)), 100);
            K_curve = zeros(size(x_range));
            
            for k = 1:length(x_range)
                K_curve(k) = eval_polynomial(K_sim_param{max_idx}, [x_range(k)]);
            end
            
            plot(x_range, K_curve, 'b-', 'LineWidth', 2);
            hold on;
            scatter(X_data_valid(:,1), K_actual_valid, 40, 'r', 'filled');
            hold off;
            
            xlabel('自变量 x1');
            ylabel('K值');
            title('拟合曲线示例', 'FontWeight', 'bold');
            grid on;
            legend('拟合曲线', '原始数据', 'Location', 'best', 'FontSize', 8);
        catch
            text(0.5, 0.5, '无法绘制拟合曲线', 'HorizontalAlignment', 'center');
            title('拟合曲线示例', 'FontWeight', 'bold');
            axis off;
        end
    elseif n_vars == 2
        % 双变量拟合曲面轮廓
        try
            x1_range = linspace(min(X_data(:,1)), max(X_data(:,1)), 30);
            x2_range = linspace(min(X_data(:,2)), max(X_data(:,2)), 30);
            [X1_grid, X2_grid] = meshgrid(x1_range, x2_range);
            
            % 计算拟合曲面
            K_surf = zeros(size(X1_grid));
            for m = 1:size(X1_grid, 1)
                for n = 1:size(X1_grid, 2)
                    K_surf(m, n) = eval_polynomial(K_sim_param{max_idx}, ...
                        [X1_grid(m, n), X2_grid(m, n)]);
                end
            end
            
            % 创建3D曲面图
            surf(X1_grid, X2_grid, K_surf, 'FaceAlpha', 0.8, ...
                'EdgeColor', 'none', 'FaceColor', 'interp');
            colormap(jet);
            hold on;

            % 绘制原始采样点
            scatter3(X_data_valid(:,1), X_data_valid(:,2), K_actual_valid, ...
                60, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
            
            % 添加网格和标签
            grid on;
            xlabel('自变量 x1', 'FontSize', 11, 'FontWeight', 'bold');
            ylabel('自变量 x2', 'FontSize', 11, 'FontWeight', 'bold');
            zlabel(sprintf('K(%d,%d)', i_best, j_best), 'FontSize', 11, 'FontWeight', 'bold');
            title('3D拟合曲面图', 'FontSize', 12, 'FontWeight', 'bold');
            
            % 设置视角
            view(45, 30); % 方位角45度，俯仰角30度
            
            % 添加颜色条
            colorbar;
            
            % 添加光照效果
            lighting gouraud;
            light('Position', [1, 1, 1], 'Style', 'infinite');

            % 添加图例
            legend('拟合曲面', '原始数据点', 'Location', 'best');
            
            hold off;
        catch
            text(0.5, 0.5, '无法绘制拟合曲面', 'HorizontalAlignment', 'center');
            title('拟合曲面轮廓图', 'FontWeight', 'bold');
            axis off;
        end
    else
        text(0.5, 0.5, sprintf('自变量数: %d\n无法直接可视化', n_vars), ...
            'FontSize', 10, 'HorizontalAlignment', 'center');
        title('多变量拟合结果', 'FontWeight', 'bold');
        axis off;
    end
    
    % 子图9: 误差统计总结
    subplot(3, 3, 9);
    
    valid_R2 = R_squared(~isnan(R_squared));
    if ~isempty(valid_R2)
        mean_r2 = mean(valid_R2);
        median_r2 = median(valid_R2);
        min_r2 = min(valid_R2);
        max_r2 = max(valid_R2);
        
        text_str = {
            sprintf('矩阵拟合质量总结');
            sprintf('====================');
            sprintf('矩阵维度: %d×%d', K_dims(1), K_dims(2));
            sprintf('有效拟合: %d/%d', length(valid_R2), n_elements);
            sprintf('平均 R²: %.4f', mean_r2);
            sprintf('中位数 R²: %.4f', median_r2);
            sprintf('最佳 R²: %.4f', max_r2);
            sprintf('最差 R²: %.4f', min_r2);
            sprintf('R² < 0.9: %.1f%%', sum(valid_R2<0.9)/length(valid_R2)*100);
            sprintf('R² < 0.95: %.1f%%', sum(valid_R2<0.95)/length(valid_R2)*100);
            sprintf('====================');
        };
    else
        text_str = {
            sprintf('矩阵拟合质量总结');
            sprintf('====================');
            sprintf('矩阵维度: %d×%d', K_dims(1), K_dims(2));
            sprintf('有效拟合: 0/%d', n_elements);
            sprintf('无有效拟合结果');
            sprintf('====================');
        };
    end
    
    text(0.05, 0.95, text_str, 'VerticalAlignment', 'top', ...
        'FontSize', 9, 'FontName', 'Monospaced');
    axis off;
    
    % 添加总标题
    sgtitle_str = sprintf('矩阵多项式拟合误差分析 (阶数: %d, 矩阵: %d×%d, 自变量: %d)', ...
        rank, K_dims(1), K_dims(2), n_vars);
    sgtitle(sgtitle_str, 'FontSize', 14, 'FontWeight', 'bold');
end

%% 子函数：绘制高维张量结果的图形
function plot_high_dimensional_results(X_data, K_flat, K_sym_flat, K_sim_param, ...
    MSE, MAE, R_squared, rank, n_elements, n_vars, K_dims)
    
    % 创建2x3的子图布局
    subplot(2, 3, 1);
    valid_R2 = R_squared(~isnan(R_squared));
    if ~isempty(valid_R2)
        histogram(valid_R2, 20, 'FaceColor', [0.2, 0.4, 0.8], 'EdgeColor', 'k');
        xlabel('决定系数 R²');
        ylabel('频数');
        title('R² 分布直方图', 'FontWeight', 'bold');
        grid on;
        hold on;
        xline(0.9, 'r--', 'R²=0.9', 'LineWidth', 1.5);
        xline(0.95, 'g--', 'R²=0.95', 'LineWidth', 1.5);
        hold off;
    else
        text(0.5, 0.5, '无有效R²数据', 'HorizontalAlignment', 'center');
        title('R² 分布直方图', 'FontWeight', 'bold');
        axis off;
    end
    
    subplot(2, 3, 2);
    valid_MAE = MAE(~isnan(MAE));
    if ~isempty(valid_MAE)
        histogram(valid_MAE, 20, 'FaceColor', [0.8, 0.4, 0.2], 'EdgeColor', 'k');
        xlabel('平均绝对误差 (MAE)');
        ylabel('频数');
        title('MAE 分布直方图', 'FontWeight', 'bold');
        grid on;
    else
        text(0.5, 0.5, '无有效MAE数据', 'HorizontalAlignment', 'center');
        title('MAE 分布直方图', 'FontWeight', 'bold');
        axis off;
    end
    
    subplot(2, 3, 3);
    valid_MSE = MSE(~isnan(MSE));
    if ~isempty(valid_MSE)
        histogram(log10(valid_MSE), 20, 'FaceColor', [0.4, 0.8, 0.2], 'EdgeColor', 'k');
        xlabel('log_{10}(MSE)');
        ylabel('频数');
        title('MSE 对数分布直方图', 'FontWeight', 'bold');
        grid on;
    else
        text(0.5, 0.5, '无有效MSE数据', 'HorizontalAlignment', 'center');
        title('MSE 对数分布直方图', 'FontWeight', 'bold');
        axis off;
    end
    
    % 选择R²最高的元素作为示例
    if ~isempty(valid_R2)
        [max_r2, max_idx] = max(valid_R2);
        
        % 获取原始K_data中对应的元素
        K_actual_best = K_flat(:, max_idx);
        
        % 计算预测值
        X_data_valid = X_data(~isnan(K_actual_best) & ~isinf(K_actual_best), :);
        K_actual_valid = K_actual_best(~isnan(K_actual_best) & ~isinf(K_actual_best));
        
        if ~isempty(K_sim_param{max_idx})
            K_pred_best = compute_polyval(K_sim_param{max_idx}, X_data_valid);
        else
            K_pred_best = zeros(size(K_actual_valid));
        end
        
        % 子图4: 拟合效果散点图
        subplot(2, 3, 4);
        scatter(K_actual_valid, K_pred_best, 30, 'filled', ...
            'MarkerFaceColor', [0.2, 0.6, 0.8], 'MarkerEdgeColor', 'k');
        hold on;
        min_val = min([K_actual_valid; K_pred_best]);
        max_val = max([K_actual_valid; K_pred_best]);
        plot([min_val, max_val], [min_val, max_val], 'r--', 'LineWidth', 2);
        hold off;
        
        xlabel('实际值');
        ylabel('预测值');
        title(sprintf('最佳元素 (R²=%.4f)', max_r2), 'FontWeight', 'bold');
        grid on;
        axis equal tight;
        
        % 子图5: 残差分析
        subplot(2, 3, 5);
        if length(K_actual_valid) > 0
            residuals = K_actual_valid - K_pred_best;
            
            % 绘制残差直方图
            yyaxis left;
            histogram(residuals, 20, 'FaceColor', [0.8, 0.2, 0.2], 'EdgeColor', 'k');
            xlabel('残差');
            ylabel('频数');
            
            % 添加正态分布曲线
            yyaxis right;
            mu = mean(residuals);
            sigma = std(residuals);
            x_vals = linspace(min(residuals), max(residuals), 100);
            y_vals = normpdf(x_vals, mu, sigma) * length(residuals) * (max(residuals)-min(residuals))/20;
            plot(x_vals, y_vals, 'b-', 'LineWidth', 2);
            ylabel('正态分布');
            
            title('残差分布与正态拟合', 'FontWeight', 'bold');
            grid on;
        else
            text(0.5, 0.5, '无有效数据', 'HorizontalAlignment', 'center');
            title('残差分布与正态拟合', 'FontWeight', 'bold');
            axis off;
        end
        
    else
        % 如果没有有效数据，显示空白图
        subplot(2, 3, 4);
        text(0.5, 0.5, '无有效拟合数据', 'HorizontalAlignment', 'center');
        title('拟合效果散点图', 'FontWeight', 'bold');
        axis off;
        
        subplot(2, 3, 5);
        text(0.5, 0.5, '无有效拟合数据', 'HorizontalAlignment', 'center');
        title('残差分析', 'FontWeight', 'bold');
        axis off;
    end
    
    % 子图6: 数据总结和信息显示
    subplot(2, 3, 6);
    
    dim_str = sprintf('%d', K_dims(1));
    for i = 2:length(K_dims)
        dim_str = sprintf('%s×%d', dim_str, K_dims(i));
    end
    
    if ~isempty(valid_R2)
        mean_r2 = mean(valid_R2);
        median_r2 = median(valid_R2);
        min_r2 = min(valid_R2);
        max_r2 = max(valid_R2);
        
        text_str = {
            sprintf('高维张量拟合总结');
            sprintf('====================');
            sprintf('张量维度: %s', dim_str);
            sprintf('元素总数: %d', n_elements);
            sprintf('有效拟合: %d (%.1f%%)', length(valid_R2), length(valid_R2)/n_elements*100);
            sprintf('平均 R²: %.4f', mean_r2);
            sprintf('中位数 R²: %.4f', median_r2);
            sprintf('最佳 R²: %.4f', max_r2);
            sprintf('最差 R²: %.4f', min_r2);
            sprintf('====================');
            sprintf('拟合阶数: %d', rank);
            sprintf('自变量数: %d', n_vars);
            sprintf('====================');
        };
    else
        text_str = {
            sprintf('高维张量拟合总结');
            sprintf('====================');
            sprintf('张量维度: %s', dim_str);
            sprintf('元素总数: %d', n_elements);
            sprintf('有效拟合: 0 (0.0%%)');
            sprintf('====================');
            sprintf('拟合阶数: %d', rank);
            sprintf('自变量数: %d', n_vars);
            sprintf('====================');
        };
    end
    
    text(0.05, 0.95, text_str, 'VerticalAlignment', 'top', ...
        'FontSize', 9, 'FontName', 'Monospaced');
    axis off;
    
    % 添加总标题
    sgtitle_str = sprintf('高维张量多项式拟合分析 (阶数: %d, 维度: %s, 自变量: %d)', ...
        rank, dim_str, n_vars);
    sgtitle(sgtitle_str, 'FontSize', 14, 'FontWeight', 'bold');
end

function y = eval_polynomial(poly_struct, x)
% 直接计算多项式值（单个样本或多个样本）
% 输入:
%   poly_struct - 多项式结构体，包含Coefficients和ModelTerms字段
%   x - 自变量值，单个样本为行向量 [1×n_vars]，多个样本为矩阵 [n_samples×n_vars]
% 输出:
%   y - 多项式值，单个样本为标量，多个样本为列向量

    if isempty(poly_struct) || ~isfield(poly_struct, 'Coefficients')
        if size(x, 1) == 1
            y = 0;
        else
            y = zeros(size(x, 1), 1);
        end
        return;
    end
    
    n_samples = size(x, 1);
    n_terms = length(poly_struct.Coefficients);
    
    if n_samples == 1
        % 单个样本情况
        y = 0;
        for j = 1:n_terms
            term = 1;
            for k = 1:size(x, 2)
                if poly_struct.ModelTerms(j, k) > 0
                    term = term * (x(k) ^ poly_struct.ModelTerms(j, k));
                end
            end
            y = y + poly_struct.Coefficients(j) * term;
        end
    else
        % 多个样本情况
        y = zeros(n_samples, 1);
        for i = 1:n_samples
            term_value = 0;
            for j = 1:n_terms
                term = 1;
                for k = 1:size(x, 2)
                    if poly_struct.ModelTerms(j, k) > 0
                        term = term * (x(i, k) ^ poly_struct.ModelTerms(j, k));
                    end
                end
                term_value = term_value + poly_struct.Coefficients(j) * term;
            end
            y(i) = term_value;
        end
    end
end