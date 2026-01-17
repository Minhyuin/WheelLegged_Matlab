function startup_WBR(robot_stats)
    % 基于LQR/MPC计算WBR模型的轮腿控制增益K矩阵，并转换为cpp/c风格代码
    % Author: MYUIN
    
    fprintf("\n=== 正在运行" + mfilename + ".m ===\n");
    
    %% 配置项 %%
    CALC_MPC = false; % 是否需要计算MPC
    
    FIT_K = true; % 是否需要拟合K
    
    %% 初始化变量 %%
    N_x = 10; % 本模型状态变量
    N_u = 4; % 本模型输出变量

    status_matrix_file = 'A10_B4.mat';
    if exist(status_matrix_file, 'file')
        fprintf('检测到已存在 %s...\n', status_matrix_file);
    else
        fprintf('%s 不存在，尝试生成...\n', status_matrix_file);
        generate_A10_B4();
        fprintf('%s 生成完毕\n', status_matrix_file);
    end
    status_matrix = load(status_matrix_file);
    
    % 通用调节矩阵
     Q=diag([ ...
                1e3,    1e0, ...   % 轮位置
                250,    200,   ...     % 底盘偏航角
                1e2,    1e2, ...       % 左腿倾角
                1e2,    1e2, ...       % 右腿倾角
                1e5,    1]);         % 底盘俯仰角
    R=diag([9e1, 9e1, 1e1, 1e1]);

    %% 拟合 %%
    if (FIT_K == true)
        fprintf("\n=== 正在拟合K矩阵" + mfilename + ".m ===\n");
        L_MIN = 0.1;  % 最小腿长
        L_MAX = 0.3;  % 最大腿长
        N_SAMPLES = 5;  % 每个维度的采样点数
        L_left_list = linspace(L_MIN, L_MAX, N_SAMPLES);
        L_right_list =  linspace(L_MIN, L_MAX, N_SAMPLES);
        TOTAL_N_SAMPLES = numel(L_left_list) * numel(L_right_list);
        K_data = zeros(TOTAL_N_SAMPLES, N_u, N_x); %  存储所有采样点对应的K值
        X_data = zeros(TOTAL_N_SAMPLES, 2);

        % 计算采样点上的K矩阵
        sample_idx = max(1, round(TOTAL_N_SAMPLES / 100));
        last_percent = 0;
        calc_start_time = tic;
        fprintf('K计算进度:')
        for l_idx = 1 : numel(L_left_list)
            length_left = L_left_list(l_idx);
            for r_idx = 1 : numel(L_right_list)
                length_right = L_right_list(r_idx);
                sample_idx = sample_idx + 1;
                current_K = calc_LQR_K40(Q, R, length_left, length_right, 1000, robot_stats, status_matrix);
                X_data(sample_idx, :) = [length_left, length_right];
                K_data(sample_idx, :, :) = current_K;

                % 更新进度显示
                current_percent = floor(sample_idx / TOTAL_N_SAMPLES * 100);
                if current_percent > last_percent
                    fprintf(repmat('#', 1, current_percent - last_percent));
                    last_percent = current_percent;
                end
            end
        end
        calc_time = toc(calc_start_time);
        fprintf('\n');
        fprintf('计算完成！\n');
        fprintf('  - 拟合总耗时: %.2f秒 (约%.2f分钟)\n', calc_time, calc_time/60);
        fprintf('  - 平均每个元素耗时: %.3f毫秒\n', calc_time/TOTAL_N_SAMPLES*1000);

        % 使用polyfitn进行拟合
        rank = 3;
        [K_sym, K_sim_param] = polyfitn_K(X_data, K_data, rank, false);
        
        % 输出代码
        polyfitn_to_cpp_optimal(K_sym, 'LQR_K_DATA');
    else
        %% 不拟合 %%
        length_left_constant = 0.2;
        length_right_constant = 0.2;
        fprintf("\n=== 正在计算定腿长K矩阵" + mfilename + ".m ===\n");
        fprintf("左腿长为" + length_left_constant + "，右腿长为" + length_right_constant + "\n");
        % 计算K
        K = calc_LQR_K40(Q, R, length_left_constant, length_right_constant, 1000, robot_stats, status_matrix);
        matrix_to_cpp(K, 'LQR_K_DATA');
    end
end