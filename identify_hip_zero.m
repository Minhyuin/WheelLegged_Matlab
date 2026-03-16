% 脚本：计算应填入底盘代码中关节电机的零值
tic
clear
clc
fprintf("\n=== 生成HIP1与HIP4的电机零值 ===\n");

%% 导入工具
currentDir = fileparts(mfilename('fullpath'));  % 获取当前脚本所在目录
commonPath = fullfile(currentDir, 'common');   % 构建 common 文件夹路径
addpath(commonPath);

%% 导入真实数据（图纸数据）
FLOOR_TO_HIP1 = 29.54 / 180 * pi; % HIP1与地面的夹角，由地面指向HIP1
HIP4_TO_FLOOR = 65.95 / 180 * pi; % HIP4与地面的夹角，由HIP4指向地面
PHI0_TO_FLOOR = 18.21 / 180 * pi; % 等效摆杆与地面夹角，由等效摆杆指向地面

%% 导入现实数据
% 并使两腿伸到最长（限位位置）, 保证底盘整体与地面平行，底盘最底部导轮与地面相切，轮面与地面相切
REAL_RF_ANGLE = 5.12457561; % 右HIP4
REAL_LF_ANGLE = 4.23826456; % 左HIP4
REAL_LB_ANGLE = 2.71684599; % 左HIP1
REAL_RB_ANGLE = 0.699921608; % 右HIP1

%% 计算零值
% 已知视线平视面向关节电机输出轴时，输出轴逆时针旋转，使得角度编码器负方向变化
% 限定坐标系，视线平视面向底盘左侧，以逆时针为图纸数据的正方向

% 对于左腿而言，HIP1电机需要顺时针旋转FLOOR_TO_HIP1, HIP4电机需要顺时针旋转PI - HIP4_TO_FLOOR
% 即HIP1与HIP4电机角度需减小
OUTPUT_LF_ANGLE = REAL_LF_ANGLE - (pi - HIP4_TO_FLOOR);
OUTPUT_LB_ANGLE = REAL_LB_ANGLE - FLOOR_TO_HIP1;
% 对于右腿而言，HIP1电机需要逆时针旋转FLOOR_TO_HIP1, HIP4电机需要逆时针旋转PI - HIP4_TO_FLOOR
% 即HIP1与HIP4电机角度需增大，但由于代码中开启is_reverse为1，所以角度需减小
OUTPUT_RF_ANGLE = REAL_RF_ANGLE - (pi - HIP4_TO_FLOOR);
OUTPUT_RB_ANGLE = REAL_RB_ANGLE - FLOOR_TO_HIP1;

% 最后需要将输出角度重新归一到[0, 2pi]
RF_0 = rangeMap(OUTPUT_RF_ANGLE, 0, 2*pi);
LF_0 = rangeMap(OUTPUT_LF_ANGLE, 0, 2*pi);
LB_180 = rangeMap(OUTPUT_LB_ANGLE, 0, 2*pi);
RB_180 = rangeMap(OUTPUT_RB_ANGLE, 0, 2*pi);

fprintf("(LF_HIP)左HIP4在0度时编码值: %f\n", LF_0);
fprintf("(LB_HIP)左HIP1在180度时编码值: %f\n", LB_180);
fprintf("(RF_HIP)右HIP4在0度时编码值: %f\n", RF_0);
fprintf("(RB_HIP)右HIP1在180度时编码值: %f\n", RB_180);
fprintf("验证PHI0值(零点位于竖直方向, 逆时针增长): %f\n", pi / 2 - PHI0_TO_FLOOR);
