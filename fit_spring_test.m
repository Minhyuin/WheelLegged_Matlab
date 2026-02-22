% 脚本：气弹簧夹角随腿长变化的3阶多项式拟合
tic
clear
clc
fprintf("\n=== 气弹簧夹角随腿长变化分析 ===\n");

%% 导入参数
currentDir = fileparts(mfilename('fullpath'));  % 获取当前脚本所在目录
commonPath = fullfile(currentDir, 'common');   % 构建 common 文件夹路径
addpath(commonPath);

robot_stats_file  = 'robot_stats.mat';
if exist(robot_stats_file, 'file')
    fprintf('检测到已存在 %s...\n', robot_stats_file);
else
    fprintf('%s 不存在，尝试生成...\n', robot_stats_file);
    stats();
    fprintf('%s 生成完毕\n', robot_stats_file);
end
robot_stats = load(robot_stats_file);

%% 数据项
L1 = robot_stats.length_thigh; % F2_s_tran 的力臂
L2 = robot_stats.length_calf;
L3 = robot_stats.length_calf;
L4 = robot_stats.length_thigh;
L5 = 0;

L1_ = 0.05035871; % F1_h 的力臂
L2_ = 0.07224957; % F2_h 的力臂
delta_alpha = 0.08314123; % 气弹簧在l2杆固定点与B点连线和l2杆的夹角

% 第一列为腿长L (m)
% 第二列为气弹簧与tau1切向方向的夹角 alpha1 (rad)
% 第三列为气弹簧与tau2切向方向的夹角 alpha2 (rad)
spring_data = [...
    0.16051326, 0.05961389, 1.02158091;
    0.161, 0.06060368, 1.01825911;
    0.171, 0.08039549, 0.95042936;
    0.181, 0.09911020, 0.88338844;
    0.191, 0.11669365, 0.81712573;
    0.201, 0.13309997, 0.75162446;
    0.211, 0.14829151, 0.68685997;
    0.221, 0.16223837, 0.62279866;
    0.231, 0.17491761, 0.55939746;
    0.241, 0.18631230, 0.49660349;
    0.251, 0.19641037, 0.43435383;
    0.261, 0.20520347, 0.37257553;
    0.271, 0.21268562, 0.31118522;
    0.281, 0.21885193, 0.25008550;
    0.291, 0.22369719, 0.18917203;
    0.301, 0.22721440, 0.12832304;
    0.311, 0.22939311, 0.06740213;
    0.3135, 0.22972696, 0.05214295;
    0.3160, 0.22997590, 0.03686707;
    0.3185, 0.23013959, 0.02157178;
    0.321, 0.23021765, 0.00625424;
    0.3235, 0.23020969, 0.00908847;
    0.3260, 0.23011524, 0.02445937;
    0.3285, 0.22993383, 0.03986161;
    0.331, 0.22966490, 0.05529848;
    0.341, 0.22770161, 0.11746259;
    0.351, 0.22428094, 0.18048176;
    0.361, 0.21933763, 0.24464896;
    0.36998937, 0.21352046, 0.30360653;
];

%% 数据处理（拟合）
L = spring_data(:, 1);          % 腿长数据 (m)
alpha1_rad = spring_data(:, 2); % tau1切向夹角数据 (rad)
alpha2_rad = spring_data(:, 3); % tau2切向夹角数据 (rad)
degree = 5;
fprintf("\n=== %d阶多项式拟合 ===\n", degree);
p1 = polyfit(L, alpha1_rad, degree);
p2 = polyfit(L, alpha2_rad, degree);
% 计算拟合值
L_fit = linspace(min(L), max(L), 200);
alpha1_fit_rad = polyval(p1, L_fit);
alpha2_fit_rad = polyval(p2, L_fit);
% 计算拟合质量指标
alpha1_pred_rad = polyval(p1, L);
alpha2_pred_rad = polyval(p2, L);
residuals1 = alpha1_rad - alpha1_pred_rad;
residuals2 = alpha2_rad - alpha2_pred_rad;
SS_res1 = sum(residuals1.^2);
SS_res2 = sum(residuals2.^2);
SS_tot1 = sum((alpha1_rad - mean(alpha1_rad)).^2);
SS_tot2 = sum((alpha2_rad - mean(alpha2_rad)).^2);
r2_1 = 1 - SS_res1 / SS_tot1;
r2_2 = 1 - SS_res2 / SS_tot2;
rmse_1 = sqrt(mean(residuals1.^2));
rmse_2 = sqrt(mean(residuals2.^2));
fprintf("\n拟合质量指标:\n");
fprintf("  R² (决定系数): %.6f\n", r2_1);
fprintf("  RMSE (均方根误差): %.6f rad\n", rmse_1);
fprintf("  R² (决定系数): %.6f\n", r2_2);
fprintf("  RMSE (均方根误差): %.6f rad\n", rmse_2);
% 绘图
figure;
plot(L, alpha1_rad, 'o', L_fit, alpha1_fit_rad, '-', 'LineWidth', 2);
hold on;
plot(L, alpha2_rad, 'o', L_fit, alpha2_fit_rad, '-', 'LineWidth', 2);
title('气弹簧夹角随腿长L变化关系', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('L (m)', 'FontSize', 12);
ylabel('\alpha (rad)', 'FontSize', 12);
legend('alpha1_dat', 'alpha1_fit', 'alpha2_dat', 'alpha2_fit', 'FontSize', 11);
fprintf("\n=== 拟合完成 ===\n");

%% 前馈分析
N = 450; % 气弹簧推力为150N
F1 = N;
F2 = N;
delta_tau2 = F2 * L2_ * cos(alpha2_fit_rad);
delta_tau1 = zeros(size(delta_tau2));
delta_torque1 = zeros(size(delta_tau2)); % 气弹簧施加到hip1的附加力矩
delta_torque4 = zeros(size(delta_tau2)); % 气弹簧施加到hip4的附加力矩
delta_tp = zeros(size(delta_tau2)); % 气弹簧引起的附加tp
delta_fd = zeros(size(delta_tau2)); % 气弹簧引起的附加fd
for i = 1:numel(L_fit)
    L0 = L_fit(i);
    alpha = acos((L1^2 + L2^2 - L0^2) / (2 * L1 * L2));
    alpha3 = alpha - delta_alpha;
    delta_tau1(i) = F2 * L1 * sin(alpha2_fit_rad(i)) * sin(alpha3) - F1 * L1_ * sin(alpha1_fit_rad(i));
    
    % 正向运动学：求解关节角度（取phi1大的无机械干涉，即第二组解）
    all_solutions = VMC_position_forward_kinematics(L0, 0, L1, L2, L3, L4, L5);
    phi1 = all_solutions(2, 1);
    phi2 = all_solutions(2, 2);
    phi3 = all_solutions(2, 3);
    phi4 = all_solutions(2, 4); 
    a21 = - (L1 * sin(phi3 - phi1)) / (L2 * sin(phi3 - phi2));
    a24 = - (L4 * sin(phi3 - phi4)) / (L2 * sin(phi3 - phi2));
    delta_torque1(i) = delta_tau1(i) + delta_tau2(i) * a21;
    delta_torque4(i) = delta_tau2(i) * a24;
    [delta_tp(i), delta_fd(i)] = VMC_torque_inverse_kinematics(delta_torque1(i), delta_torque4(i), phi1, phi4, L1, L2, L3, L4, L5);
end
% 绘制前馈力矩曲线
figure;
plot(L_fit, delta_torque1, 'r', 'LineWidth', 2);
hold on;
plot(L_fit, delta_torque4, 'b', 'LineWidth', 2);
title('气弹簧前馈力矩', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('L (m)', 'FontSize', 12);
ylabel('Torque (Nm)', 'FontSize', 12);
legend('△torque1', '△torque4', 'FontSize', 11);
% 绘制前馈tp和fd曲线
figure;
subplot(2,1,1);
plot(L_fit, delta_tp, 'm', 'LineWidth', 2);
title('气弹簧引起的附加tp', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('L (m)', 'FontSize', 12);
ylabel('tp (Nm)', 'FontSize', 12);
subplot(2,1,2);
plot(L_fit, delta_fd, 'c', 'LineWidth', 2);
title('气弹簧引起的附加fd', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('L (m)', 'FontSize', 12);
ylabel('fd (Nm)', 'FontSize', 12);

fprintf("\n=== 前馈分析完成 ===\n");

%% 输出表达式
prec = 6;   % 精度
expr1 = poly2cstring(p1, 'l', prec);
expr2 = poly2cstring(p2, 'l', prec);
fprintf('float alpha1 = %s\n', expr1);
fprintf('float alpha2 = %s\n', expr2);

%% 显示执行时间
toc
fprintf("\n=== 分析完成 ===\n");

%% 辅助打印输出函数
function expr = poly2cstring(p, var, prec)
    % p       : 多项式系数，从高次到低次（polyfit 输出格式）
    % var     : 自变量符号，例如 'l'
    % prec    : 小数位数
    % expr    : 生成的表达式字符串，例如 '1.23f*lll - 0.45f*ll + 6.78f*l - 9.01f'
    n = length(p) - 1;          % 多项式的最高次数
    terms = {};                  % 存储各项字符串
    for k = 0:n
        coef = p(n - k + 1);     % 取出 k 次项系数（k=0 为常数项）
        if abs(coef) < 1e-12     % 忽略极小的系数（可调整阈值）
            continue;
        end
        % 格式化系数
        if abs(round(coef) - coef) < 1e-12    % 判断是否为整数
            coef_str = sprintf('%d', round(coef));
        else
            coef_str = sprintf(['%.', num2str(prec), 'f'], coef);
            coef_str = [coef_str, 'f'];       %#ok<AGROW> % 末尾加 f
        end
        % 构造项字符串
        if k == 0
            term = coef_str;                   % 常数项
        else
            pow_str = repmat(var, 1, k);       % 幂次表示：l, ll, lll, ...
            if abs(coef - 1) < 1e-12
                term = pow_str;                 % 系数为 1 时省略
            elseif abs(coef + 1) < 1e-12
                term = ['-', pow_str];          % 系数为 -1 时保留负号
            else
                term = [coef_str, '*', pow_str];
            end
        end
        terms{end+1} = term; %#ok<AGROW>
    end
    % 将各项用加号连接，处理负号
    expr = terms{1};
    for i = 2:length(terms)
        if terms{i}(1) == '-'
            expr = [expr, terms{i}];            %#ok<AGROW> % 负项直接拼接
        else
            expr = [expr, '+', terms{i}]; %#ok<AGROW>
        end
    end
end