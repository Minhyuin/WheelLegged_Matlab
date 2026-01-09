% 基于LQR/MPC计算WBR模型的轮腿控制增益K矩阵，并转换为cpp/c风格代码
% Author: MYUIN

%% START %%
tic
clear
clc
fprintf("\n=== 正在运行" + mfilename + ".m ===\n");
disp("清空命令行窗口");
disp("开始计时");

%% 配置项 %%
MODEL = 0; % 模型选择，0: WBR(全身建模) | 1: WL(单腿建模) | 2: ST(板凳建模)

%% 加载必要物理参数 %%
robot_stats_file  = 'robot_stats.mat';
if exist(robot_stats_file, 'file')
    fprintf('检测到已存在 %s...\n', robot_stats_file);
    robot_stats = load(robot_stats_file);
else
    fprintf('%s 不存在，尝试生成...\n', robot_stats_file);
    stats();
    fprintf('%s 生成完毕\n', robot_stats_file);
end
robot_stats = load(robot_stats_file);

if MODEL == 0
    startup_WBR(robot_stats);
end

%% END %%
fprintf("\n=== 解算耗时:" + toc + "===\n");
