function stats()
    % 本文件将会生成所有轮腿机器人的物理参量及数据矩阵
    % 同时，本文件会拟合出腿部物理量随腿长变化的轨迹
    % 工作区内容为空时，应先运行此文件
    % Author: MYUIN
    % Note : 全部为国际单位

    fprintf("\n=== 正在运行" + mfilename + ".m ===\n");
    
    %% 常数物理量 %%
    g_ = 9.788; % 在广州，重力加速度为9.788m/s^2，加入_后缀以区分matlab定义的g
    
    %% 机器人整体 %%
    % 轮腿机器人可以髋部转轴为关节点区分开【腿】和【机体】两部分
    % 【腿】包含了从髋关节开始以下的所有部分，包含大腿、膝关节、小腿、驱动轮（踝关节）
    % 【机体】即机器人整机除去腿外的部分，包含底盘、云台、发射机构
    
    mass_total = 23.735; % 整机称重可得
    
    inertia_z = 0.664166057; % 机器人绕z轴（中心）转动惯量 （实际在现实中并非常量，这里简化为常量
    
    distance_body_hip = -0.00345; % 机体质心到髋关节的距离
    
    distance_wheel = 0.4611; % 左右驱动轮轮距
    distance_half_wheel = distance_wheel / 2; % 左右驱动轮轮距 / 2
    
    %% 腿 %%
    %%% 非变化常量 %%%
    length_thigh = 0.210; % 大腿长;
    length_calf = 0.250; % 小腿长;

    mass_thigh_aver = 0.37838157; % 单根大腿质量（几何上等效）
    mass_calf_aver = 0.36700000; % 单根小腿质量 (几何上等效)
    mass_leg_aver = 1.43744054; % 假设两腿等价，单腿的质量
    mass_leg_other_aver = mass_leg_aver - mass_thigh_aver - mass_calf_aver; % 除主结构大腿和小腿外的质量（如轴承、额外连杆等）
    mass_leg_left = mass_leg_aver; % 左腿整腿质量
    mass_leg_right = mass_leg_aver; % 右腿整腿质量
    
    mass_wheel_motor_aver = 0.334; % 轮毂电机质量
    mass_wheel_other_aver = 0.179; % 轮毂其他加工件质量
    mass_wheel_aver = mass_wheel_motor_aver + mass_wheel_other_aver; % 假设两轮等价，单轮的质量 
    mass_wheel_left = mass_wheel_aver; % 左轮毂质量
    mass_wheel_right = mass_wheel_aver; % 右轮毂质量
    
    mass_leg_total = mass_leg_left + mass_leg_right + mass_wheel_left + mass_wheel_right;
    
    diameter_wheel = 0.11704; % 驱动轮直径
    radius_wheel = diameter_wheel / 2; % 驱动轮半径
    
    % inertia_wheel_aver =  0.5 * mass_wheel_aver * radius_wheel^2;
    inertia_wheel_aver = 0.00063738; % 驱动轮绕轮心转动惯量
    
    %%% 随腿运动学变化而变化的量 %%%
    
    leg_fit_n  = 3; % 拟合三阶
    
    % 选定一只腿（假设左腿），根据不同腿长长度，在solidworks上量测出以上值
    % 第一列为腿长，步进为0.01m
    % 第二列为驱动轮质心到腿摆杆质心距离 distance_wheel_leg
    % 第三列为机体转轴到腿摆杆质心距离 distance_hip_leg
    % 第四列为腿摆杆绕髋关节转动惯量 inertia_leg
    leg_data_left = [   0.16051326, 0.12104874, 0.09606230, 0.04553921;
                        0.161,  0.12116894,  0.09626696,  0.04563842;
                        0.171,  0.12366419,  0.10052335,  0.04774519;
                        0.181,  0.12621225,  0.10486884,  0.04998269;
                        0.191,  0.12881462,  0.10929144,  0.05235107;
                        0.201,  0.13146901,  0.11378131,  0.05485053;
                        0.211,  0.13417079,  0.11833030,  0.05748124;
                        0.221,  0.13691375,  0.12293159,  0.06024338;
                        0.231,  0.13969060,  0.12757942,  0.06313715;
                        0.241,  0.14249327,  0.13226890,  0.06616272;
                        0.251,  0.14531306,  0.13699585,  0.06932031;
                        0.261,  0.14814065,  0.14175666,  0.07261011;
                        0.271,  0.15096613,  0.14654822,  0.07603232;
                        0.281,  0.15377891,  0.15136782,  0.07958718;
                        0.291,  0.15656758,  0.15621312,  0.08327491;
                        0.301,  0.15931976,  0.16108205,  0.08709576;
                        0.311,  0.16202185,  0.16597283,  0.09105000;
                        0.321,  0.16465872,  0.17088389,  0.09513793;
                        0.331,  0.16721331,  0.17581389,  0.09935988;
                        0.341,  0.16966607,  0.18076166,  0.10371623;
                        0.351,  0.17199420,  0.18572622,  0.10820739;
                        0.361,  0.17417062,  0.19070678,  0.11283389;
                        0.36998937, 0.17597058, 0.19519711, 0.11710883];
    
    % 假设两腿完全等价，拟合数据不分左右
    [coeff_distance_wheel_leg, coeff_distance_hip_leg, coeff_inertia_leg] = polyfit_leg_data(leg_data_left, ...
        'DegreeDisWl', leg_fit_n, ...
        'DegreeDisHl', leg_fit_n, ... 
        'DegreeIl', leg_fit_n, ...
        'ShowPlot', true, ...
        'SaveCoeff', false, ...
        'FileName', 'leg_coeff_left.mat'); 
    
    %% 机体 %%
    mass_body = mass_total - mass_leg_total; % 机体质量
    inertia_body = 0.397912180; % 机体绕髋关节转动惯量
    body_centroid_offset_ang = 1.71356828; % 机体质心偏移角度（从机体方向指向质心，逆时针正方向）
    
    % 底盘三维
    chassis_length = 0.51879;
    chassis_width = 0.56580;
    chassis_height = 0.19300;
    
    % 云台三维
    gimbal_length = 0.3;
    gimbal_width = 0.3;
    gimbal_height = 0.3;
    
    %% 物理极限 %%
    max_torque_wheel = 3.94; % 轮毂最大连续输出力矩
    max_rpm_wheel = 280; % 轮毂最大瞬时转速

    max_torque_joint = 40; % 腿关节最大连续输出力矩
    max_rpm_joint = 300; % 腿关节最大瞬时转速

    leg_length_max = 0.36998937;
    leg_length_min = 0.16051326;
    
    %% 存参 %%
    save robot_stats.mat 
end