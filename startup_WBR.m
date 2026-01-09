function startup_WBR(robot_stats)
    % 基于LQR/MPC计算WBR模型的轮腿控制增益K矩阵，并转换为cpp/c风格代码
    % Author: MYUIN
    
    fprintf("\n=== 正在运行" + mfilename + ".m ===\n");
    
    currentDir = fileparts(mfilename('fullpath'));  % 获取当前脚本所在目录
    commonPath = fullfile(currentDir, 'common');   % 构建 common 文件夹路径
    addpath(commonPath);
    
    %% 配置项 %%
    CALC_MPC = false; % 是否需要计算MPC
    
    FIT_K = false; % 是否需要拟合K
    
    %% 加载状态矩阵 %%
    status_matrix_file = 'A10_B4.mat';
    if exist(status_matrix_file, 'file')
        fprintf('检测到已存在 %s...\n', status_matrix_file);
    else
        fprintf('%s 不存在，尝试生成...\n', status_matrix_file);
        generate_A10_B4();
        fprintf('%s 生成完毕\n', status_matrix_file);
    end
    status_matrix = load(status_matrix_file);
    
    %% 测试项 %%
    Q=diag([ ...
                1e3,    1e0, ...   % 轮位置
                250,    200,   ...     % 底盘偏航角
                1e2,    1e2, ...       % 左腿倾角
                1e2,    1e2, ...       % 右腿倾角
                1e5,    1]);         % 底盘俯仰角
    R=diag([9e1, 9e1, 1e1, 1e1]);
    
    K = calc_LQR_K40(Q, R, 0.2, 0.2, 1000, robot_stats, status_matrix);
    disp(K);
end