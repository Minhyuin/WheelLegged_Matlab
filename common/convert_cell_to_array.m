function param_array = convert_cell_to_array(param_cell, N_u, N_x, rank, n_vars)
%   将元胞数组格式的K_sim_param转换为数组格式
%   适用于Simulink中的Constant模块
%
%   输入参数:
%   param_cell - 元胞数组格式的多项式系数，每个元素是一个结构体
%   N_u - K矩阵行数
%   N_x - K矩阵列数
%   rank - 多项式阶数
%   n_vars - 自变量数量（1或2）
%
%   输出参数:
%   param_array - 数组格式：
%       当 n_vars = 1 时: 3维数组 [N_u, N_x, rank+1]
%       当 n_vars = 2 时: 4维数组 [N_u, N_x, rank+1, rank+1]

    fprintf('\n=== 正在转换系数格式 ===\n');
    
    % 检查输入
    if nargin < 5
        error('需要5个输入参数: param_cell, N_u, N_x, rank, n_vars');
    end
    
    if ~ismember(n_vars, [1, 2])
        error('自变量数量 n_vars 必须为 1 或 2，当前为 %d', n_vars);
    end
    
    if length(param_cell) ~= N_u * N_x
        error('元胞数组长度 %d 与矩阵大小 %d×%d=%d 不匹配', ...
            length(param_cell), N_u, N_x, N_u*N_x);
    end
    
    % 创建数组
    if n_vars == 1
        % 一元多项式：3维数组 [N_u, N_x, rank+1]
        param_array = zeros(N_u, N_x, rank+1);
        fprintf('创建3维数组: [%d, %d, %d] (一元多项式)\n', N_u, N_x, rank+1);
    else
        % 二元多项式：4维数组 [N_u, N_x, rank+1, rank+1]
        param_array = zeros(N_u, N_x, rank+1, rank+1);
        fprintf('创建4维数组: [%d, %d, %d, %d] (二元多项式)\n', N_u, N_x, rank+1, rank+1);
    end
    
    % 用于显示进度
    total_elements = N_u * N_x;
    fprintf('转换进度: ');
    progress_interval = max(1, round(total_elements / 50));
    
    for elem_idx = 1:total_elements
        % 获取二维索引
        [i, j] = ind2sub([N_u, N_x], elem_idx);
        
        % 获取当前元素的多项式结构
        poly_struct = param_cell{elem_idx};
        
        if isempty(poly_struct) || ~isfield(poly_struct, 'Coefficients')
            % 如果结构体为空或缺少Coefficients字段，跳过
            continue;
        end
        
        % 提取多项式的项和系数
        if isfield(poly_struct, 'ModelTerms')
            expon = poly_struct.ModelTerms;
        else
            % 如果没有ModelTerms字段，使用默认值（常数项）
            expon = zeros(1, n_vars);
        end
        
        coef = poly_struct.Coefficients;
        nterms = length(coef);
        
        % 将系数存储到数组
        for k = 1:nterms
            if n_vars == 1
                % 一元多项式情况
                if size(expon, 2) >= 1
                    exp1 = expon(k, 1) + 1;  % MATLAB索引从1开始
                    
                    % 检查指数是否在范围内
                    if exp1 <= rank + 1
                        param_array(i, j, exp1) = coef(k);
                    else
                        warning('元素 [%d,%d] 的第 %d 项指数 %d 超过了阶数限制 %d', ...
                            i, j, k, exp1-1, rank);
                        
                        % 如果指数超过了范围，将其放到最后一个位置
                        exp1_adj = min(exp1, rank+1);
                        param_array(i, j, exp1_adj) = coef(k);
                    end
                else
                    % 如果没有指数信息，假设是常数项
                    param_array(i, j, 1) = coef(k);
                end
            else
                % 二元多项式情况
                if size(expon, 2) >= 2
                    exp1 = expon(k, 1) + 1;  % MATLAB索引从1开始
                    exp2 = expon(k, 2) + 1;
                    
                    % 检查指数是否在范围内
                    if exp1 <= rank + 1 && exp2 <= rank + 1
                        param_array(i, j, exp1, exp2) = coef(k);
                    else
                        warning('元素 [%d,%d] 的第 %d 项指数 (%d,%d) 超过了阶数限制 %d', ...
                            i, j, k, exp1-1, exp2-1, rank);
                        
                        % 如果指数超过了范围，将其放到最后一个位置
                        exp1_adj = min(exp1, rank+1);
                        exp2_adj = min(exp2, rank+1);
                        param_array(i, j, exp1_adj, exp2_adj) = coef(k);
                    end
                else
                    % 如果指数维度不足，假设是常数项
                    param_array(i, j, 1, 1) = coef(k);
                end
            end
        end
        
        % 显示进度
        if mod(elem_idx, progress_interval) == 0 || elem_idx == total_elements
            fprintf('#');
        end
    end
    
    fprintf('\n转换完成！\n');
    
    % 验证转换结果
    nonzero_count = nnz(param_array);
    total_size = numel(param_array);
    fprintf('非零系数比例: %.2f%% (%d/%d)\n', ...
        nonzero_count/total_size*100, nonzero_count, total_size);
    
    % 显示数组维度
    array_dims = size(param_array);
    fprintf('输出数组维度: [');
    for d = 1:length(array_dims)
        fprintf('%d', array_dims(d));
        if d < length(array_dims)
            fprintf(', ');
        end
    end
    fprintf(']\n');
end