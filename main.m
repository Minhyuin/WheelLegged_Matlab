% 主程序，按配置调用不同模型以生成能够在嵌入式端使用的增益参数
% Author: MYUIN

%% START %%
tic
clear
clc
fprintf("\n=== 正在运行" + mfilename + ".m ===\n");

currentDir = fileparts(mfilename('fullpath'));  % 获取当前脚本所在目录
commonPath = fullfile(currentDir, 'common');   % 构建 common 文件夹路径
addpath(commonPath);
simPath = fullfile(currentDir, 'sim');   % 构建 common 文件夹路径
addpath(simPath);

disp("清空命令行窗口");
disp("开始计时");

%% 配置项 %%
MODEL = 0; % 模型选择，0: WBR(全身建模) | 1: WL(单腿建模) | 2: ST(板凳建模) | 3：WBR2(优化平衡点的全身建模)
SIM_ENABLE = true; % 是否需要参与仿真

%% 加载必要物理参数 %%
robot_stats_file  = 'robot_stats.mat';
if exist(robot_stats_file, 'file')
    fprintf('检测到已存在 %s...\n', robot_stats_file);
else
    fprintf('%s 不存在，尝试生成...\n', robot_stats_file);
    stats();
    fprintf('%s 生成完毕\n', robot_stats_file);
end
robot_stats = load(robot_stats_file);

%% 启动模型 %%
if MODEL == 0
    startup_WBR_class(robot_stats);
elseif MODEL == 1
    startup_WL_class(robot_stats);
elseif MODEL == 2
    startup_ST_class(robot_stats);
elseif MODEL == 3
    startup_WBR2_class(robot_stats);
else
    error('\n！未知模型！\n');
end

%% 加载仿真参数 %%
if SIM_ENABLE == true
    SIMPARAM = load_sim(MODEL, false);
end

%% END %%
fprintf("\n=== 总耗时:" + toc + "===\n");

%% WBR启动函数
function startup_WBR_class(robot_stats)
    controller = WBR(robot_stats);
    controller.computeAndFitK();
    controller.saveSimParam('WBR_sim_param.mat');
    fprintf('\n=== WBR控制器设计完成 ===\n');
end

%% WL启动函数
function startup_WL_class(robot_stats)
    controller = WL(robot_stats);
    controller.computeAndFitK();
    controller.saveSimParam('WL_sim_param.mat');
    fprintf('\n=== WL控制器设计完成 ===\n');
end

%% ST启动函数
function startup_ST_class(robot_stats)
    controller = ST(robot_stats);
    controller.computeAndFitK();
    controller.saveSimParam('ST_sim_param.mat');
    fprintf('\n=== ST控制器设计完成 ===\n');
end
 
%% WBR2启动函数
function startup_WBR2_class(robot_stats)
    controller = WBR2(robot_stats);
    controller.computeAndFitK();
    controller.saveSimParam('WBR2_sim_param.mat');
    fprintf('\n=== WBR2控制器设计完成 ===\n');
end