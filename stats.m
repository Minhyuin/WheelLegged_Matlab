% 本文件将会生成所有轮腿机器人的物理参量及数据矩阵
% 同时，本文件会拟合出腿部物理量随腿长变化的轨迹
% 工作区内容为空时，应先运行此文件
% Author: MYUIN
% Note : 全部为国际单位

%% START %%
tic
clear
clc
fprintf("\n=== 正在运行" + mfilename + ".m ===\n");
disp("清空命令行窗口");
disp("开始计时");

%% 常数物理量 %%
g_ = 9.788; % 在广州，重力加速度为9.788m/s^2，加入_后缀以区分matlab定义的g

%% 机器人整体 %%
% 轮腿机器人可以髋部转轴为关节点区分开【腿】和【机体】两部分
% 【腿】包含了从髋关节开始以下的所有部分，包含大腿、膝关节、小腿、驱动轮（踝关节）
% 【机体】即机器人整机除去腿外的部分，包含底盘、云台、发射机构

mass_total = 15; % 整机称重可得

inertia_z = 1; % 机器人绕z轴（中心）转动惯量 （实际在现实中并非常量，这里简化为常量）

distance_wheel = 1; % 左右驱动轮轮距
distance_half_w = distance_wheel / 2; % 左右驱动轮轮距 / 2

%% 腿 %%
%%% 非变化常量 %%%
mass_thigh_aver = 1; % 单根大腿质量（几何上等效）
mass_calf_aver = 1; % 单根小腿质量 (几何上等效)
mass_leg_other_aver = 1; % 除主结构大腿和小腿外的质量（如轴承、额外连杆等）
mass_leg_aver = mass_thigh_aver + mass_calf_aver + mass_leg_other_aver; % 假设两腿等价，单腿的质量
mass_leg_left = mass_leg_aver; % 左腿整腿质量
mass_leg_right = mass_leg_aver; % 右腿整腿质量

mass_wheel_motor_aver = 1; % 轮毂电机质量
mass_wheel_other_aver = 1; % 轮毂其他加工件质量
mass_wheel_aver = mass_wheel_motor_aver + mass_wheel_other_aver; % 假设两轮等价，单轮的质量 
mass_wheel_left = mass_wheel_aver; % 左轮毂质量
mass_wheel_right = mass_wheel_aver; % 右轮毂质量

mass_leg_total = mass_leg_left + mass_leg_right + mass_wheel_left + mass_wheel_right;

diameter_wheel = 1; % 驱动轮直径
radius_wheel = diameter_wheel / 2; % 驱动轮半径

% Inertial_wheel_aver =  0.5 * mass_wheel_aver * radius_wheel^2;
inertia_wheel_aver = 1; % 驱动轮绕轮新转动惯量

%%% 随腿运动学变化而变化的量 %%%

% 以下变量会因腿长变化而变，故采样拟合
% 驱动轮质心到腿摆杆质心距离 dis_w_l
% 机体转轴到腿摆杆质心距离 dis_b_l
% 腿摆杆绕髋关节转动惯量 I_l

leg_fit_n  = 3; % 拟合三阶

% 选定一只腿（假设左腿），根据不同腿长长度，在solidworks上量测出以上值
% 第一列为腿长，步进为0.01m
% 第二列为distance_wheel_leg
% 第三列为distance_body_leg
% 第四列为inertia_leg
leg_data_left = [0.11,  0.0480,  0.0620,  0.01819599;
                      0.12,  0.0470,  0.0730,  0.01862845;
                      0.13,  0.0476,  0.0824,  0.01898641;
                      0.14,  0.0480,  0.0920,  0.01931342;
                      0.15,  0.0490,  0.1010,  0.01962521;
                      0.16,  0.0500,  0.1100,  0.01993092;
                      0.17,  0.0510,  0.1190,  0.02023626;
                      0.18,  0.0525,  0.1275,  0.02054500;
                      0.19,  0.0539,  0.1361,  0.02085969;
                      0.20,  0.0554,  0.1446,  0.02118212;
                      0.21,  0.0570,  0.1530,  0.02151357;
                      0.22,  0.0586,  0.1614,  0.02185496;
                      0.23,  0.0600,  0.1700,  0.02220695;
                      0.24,  0.0621,  0.1779,  0.02256999;
                      0.25,  0.0639,  0.1861,  0.02294442;
                      0.26,  0.0657,  0.1943,  0.02333041;
                      0.27,  0.0676,  0.2024,  0.02372806;
                      0.28,  0.0700,  0.2100,  0.02413735;
                      0.29,  0.0713,  0.2187,  0.02455817;
                      0.30,  0.0733,  0.2267,  0.02499030];

% 假设两腿完全等价，拟合数据不分左右
[coeff_distance_wheel_leg, coeff_distance_body_leg, coeff_inertia_leg] = polyfit_leg_data(leg_data_left, ...
    'DegreeDisWl', leg_fit_n, ...
    'DegreeDisBl', leg_fit_n, ... 
    'DegreeIl', leg_fit_n, ...
    'ShowPlot', true, ...
    'SaveCoeff', false, ...
    'FileName', 'leg_coeff_left.mat'); 

%% 机体 %%
mass_body = mass_total - mass_leg_total; % 机体质量
inertia_body = 1; % 机体绕髋关节转动惯量

% 底盘三维
chassis_length = 0.3;
chassis_width = 0.3;
chassis_height = 0.3;

% 云台三维
gimbal_length = 0.3;
gimbal_width = 0.3;
gimbal_height = 0.3;


%% 存参 %%
fprintf("\n=== 解算耗时:" + toc + "===\n");

save robot_stats.mat 