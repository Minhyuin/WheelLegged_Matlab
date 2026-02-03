classdef WBR < handle
    %WBR WBR模型控制器设计类
    %   基于LQR/MPC计算WBR模型的轮腿控制增益K矩阵，并转换为C++风格代码
    
    properties
        % 控制器配置
        CALC_MPC logical = false   % 是否需要计算MPC
        FIT_K logical = true       % 是否需要拟合K
        EVALUATE_FITTING_STATS logical = false; % 是否需要评估拟合结果
        
        % 机器人参数
        robot_stats struct         % 机器人参数结构体
        
        % 状态矩阵数据
        status_matrix struct       % 从文件加载的状态矩阵

        % 计算结果
        K                          % LQR增益矩阵
        K_data                     % 拟合数据
        X_data                     % 腿长采样数据
        K_sym                      % 符号化拟合结果
        K_sim_param          % 拟合参数
    end
    
    properties (SetAccess = private)
        % 模型维度
        N_x double = 10            % 状态变量维度
        N_u double = 4             % 控制输入维度
        
        % 权重矩阵
        Q                          % 状态权重矩阵
        R                          % 控制权重矩阵

        % 控制频率
        control_freq double = 1000

        % 拟合参数
        N_SAMPLES double = 5       % 每个维度的采样点数
        poly_rank double = 3       % 多项式拟合阶数

        % 默认腿长
        default_length = 0.2;
    end
    
    methods
        function obj = WBR(robot_stats)
            %WBR 构造函数
            %   obj = WBR(robot_stats) 初始化控制器设计
            
            fprintf("\n=== 初始化WBR控制器设计 ===\n");
            
            % 保存机器人参数
            obj.robot_stats = robot_stats;
            
            % 初始化权重矩阵
            obj.initializeWeightMatrices();
            
            % 加载或生成状态矩阵
            obj.loadOrGenerateStatusMatrix();
        end
        
        function initializeWeightMatrices(obj)
            %INITIALIZEWEIGHTMATRICES 初始化权重矩阵
            obj.Q = diag([ ...
                1000, 2, ...       % 轮位置
                250, 200, ...       % 底盘偏航角
                5000, 600, ...       % 左腿倾角
                5000, 600, ...       % 右腿倾角
                100000, 1]);           % 底盘俯仰角

            obj.R = diag([ ...
                9e1, ...    % 左轮能量惩罚
                9e1, ...    % 右轮能量惩罚
                1e1, ...    % 左腿能量惩罚
                1e1]);      % 右腿能量惩罚
        end
        
        function loadOrGenerateStatusMatrix(obj)
            %LOADORGENERATESTATUSMATRIX 加载或生成状态矩阵文件
            
            status_matrix_file = 'A10_B4.mat';
            
            if exist(status_matrix_file, 'file')
                fprintf('检测到已存在 %s...\n', status_matrix_file);
                obj.status_matrix = load(status_matrix_file);
            else
                fprintf('%s 不存在，尝试生成...\n', status_matrix_file);
                WBR_generate_state_space();
                obj.status_matrix = load(status_matrix_file);
                fprintf('%s 生成完毕\n', status_matrix_file);
            end
        end
        
        function K = computeLQRForSingleLength(obj, length_left, length_right)
            %COMPUTELQRFORSINGLELENGTH 计算单腿长下的LQR增益
            
            fprintf("计算LQR增益 - 左腿长: %.3f, 右腿长: %.3f\n", ...
                length_left, length_right);
            
            K = obj.calc_LQR_K(obj.Q, obj.R, ...
                length_left, length_right, obj.control_freq, ...
                obj.robot_stats, obj.status_matrix);
        end
        
        function computeAndFitK(obj)
            %COMPUTEANDFITK 计算并拟合K矩阵
            
            if ~obj.FIT_K
                fprintf("\n=== 直接计算固定腿长K矩阵 ===\n");
                length_left = obj.default_length;
                length_right = obj.default_length;
                fprintf("左腿长: %.3f, 右腿长: %.3f\n", length_left, length_right);
                
                obj.K = obj.computeLQRForSingleLength(length_left, length_right);
                
                % 转换为C++代码
                obj.matrixToCpp(obj.K, 'LQR_K_DATA');
                return;
            end
            
            fprintf("\n=== 开始拟合K矩阵 ===\n");
            
            % 生成腿长采样点
            L_left_list = linspace(obj.robot_stats.leg_length_min, obj.robot_stats.leg_length_max, obj.N_SAMPLES);
            L_right_list = linspace(obj.robot_stats.leg_length_min, obj.robot_stats.leg_length_max, obj.N_SAMPLES);
            
            TOTAL_N_SAMPLES = numel(L_left_list) * numel(L_right_list);
            
            % 预分配存储空间
            obj.K_data = zeros(TOTAL_N_SAMPLES, obj.N_u, obj.N_x);
            obj.X_data = zeros(TOTAL_N_SAMPLES, 2);
            
            % 计算进度显示
            sample_idx = 0;
            calc_start_time = tic;
            
            fprintf('K计算进度: 0%%');
            
            % 遍历所有采样点
            for l_idx = 1:numel(L_left_list)
                length_left = L_left_list(l_idx);
                
                for r_idx = 1:numel(L_right_list)
                    length_right = L_right_list(r_idx);
                    sample_idx = sample_idx + 1;
                    
                    % 计算当前腿长下的K矩阵
                    current_K = obj.calc_LQR_K(obj.Q, obj.R, ...
                        length_left, length_right, obj.control_freq, ...
                        obj.robot_stats, obj.status_matrix);
                    
                    % 保存数据
                    obj.X_data(sample_idx, :) = [length_left, length_right];
                    obj.K_data(sample_idx, :, :) = current_K;
                    
                    % 更新进度显示
                    progress_percent = floor(sample_idx / TOTAL_N_SAMPLES * 100);
                    if mod(sample_idx, max(1, round(TOTAL_N_SAMPLES/100))) == 0
                        fprintf('\b\b\b\b%3d%%', progress_percent);
                    end
                end
            end
            
            calc_time = toc(calc_start_time);
            fprintf('\n');
            
            % 显示计算统计
            obj.displayComputationStats(calc_time, TOTAL_N_SAMPLES);
            
            % 多项式拟合
            obj.fitKMatrix();
            
            % 生成C++代码
            obj.polyfitnToCppOptimal(obj.K_sym, 'LQR_K_DATA');
        end
        
        function displayComputationStats(~, calc_time, total_samples)
            %DISPLAYCOMPUTATIONSTATS 显示计算统计信息
            
            fprintf('计算完成！\n');
            fprintf('  - 拟合总耗时: %.2f秒 (约%.2f分钟)\n', ...
                calc_time, calc_time/60);
            fprintf('  - 总采样点数: %d\n', total_samples);
            fprintf('  - 平均每个点耗时: %.3f毫秒\n', ...
                calc_time/total_samples*1000);
        end
        
        function fitKMatrix(obj)
            
            %FITKMATRIX 使用多项式拟合K矩阵
            
            fprintf('正在进行多项式拟合 (阶数: %d)...\n', obj.poly_rank);
            
            % 调用拟合函数
            [obj.K_sym, K_sim_cell] = ...
                obj.polyfitn_matrix(obj.X_data, obj.K_data, obj.poly_rank, obj.EVALUATE_FITTING_STATS);
            
            obj.K_sim_param = convert_cell_to_array(K_sim_cell, obj.N_u, obj.N_x, obj.poly_rank, 2);

            fprintf('多项式拟合完成！\n');
        end

        function saveSimParam(obj, filename)
            %SAVESIMPARAM 保存仿真所需的拟合系数计算结果
            WBR_SIM = struct();
            WBR_SIM.timestamp = datetime('now');
            WBR_SIM.robot_stats = obj.robot_stats;
            WBR_SIM.K_sim_param = obj.K_sim_param;
            
            save(filename, 'WBR_SIM');
            fprintf('拟合系数参数已保存到: %s\n', filename);
        end
    end
    
    methods (Static)
        function matrixToCpp(matrix, var_name)
            %MATRIXTOCPP 将矩阵转换为C++代码
            %   保持与原始函数相同的接口
            
            fprintf('\n转换为C++代码:\n');
            fprintf('变量名: %s\n', var_name);
            fprintf('矩阵维度: %d x %d\n', size(matrix, 1), size(matrix, 2));
            
            matrix_to_cpp(matrix, var_name);
        end
        
        function polyfitnToCppOptimal(K_sym, var_name)
            %POLYFITNTOCPPOPTIMAL 将多项式拟合结果转换为C++代码
            %   保持与原始函数相同的接口
            
            fprintf('\n转换多项式拟合结果为C++代码:\n');
            fprintf('变量名: %s\n', var_name);
            
            polyfitn_to_cpp_optimal(K_sym, var_name);
        end
    end
    
    methods (Access = private)
        % 这些方法假设已经在路径中定义
        function K = calc_LQR_K(~, Q, R, length_left, length_right, freq, robot_stats, status_matrix)
            K = WBR_calc_LQR_K(Q, R, length_left, length_right, ...
                freq, robot_stats, status_matrix);
        end
        
        function [K_sym, K_sim_param] = polyfitn_matrix(~, X_data, K_data, rank, verbose)
            %POLYFITN_K 多项式拟合K矩阵（包装外部函数）
            [K_sym, K_sim_param] = polyfitn_matrix(X_data, K_data, rank, verbose);
        end
    end
end

