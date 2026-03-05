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
% 第一列为腿长L (m)
% 第二列为气弹簧与腿长方向夹角 alpha (rad)
spring_data = [...
    0.16051326, 1.52115247;
    0.161,  1.51816648;
    0.171,  1.45830378;
    0.181,  1.40106219;
    0.191,  1.34618555;
    0.201,  1.29346758;
    0.211,  1.24273685;
    0.221,  1.19384638;
    0.231,  1.14666630;
    0.241,  1.10107873;
    0.251,  1.05697401;
    0.261,  1.01424804;
    0.271,  0.97280025;
    0.281,  0.93253200;
    0.291,  0.89334522;
    0.301,  0.85514114;
    0.311,  0.81781900;
    0.321,  0.78127446;
    0.331,  0.74539776;
    0.341,  0.71007126;
    0.351,  0.67516606;
    0.361,  0.64053724;
    0.36998937, 0.60950614];

%% 数据处理
L = spring_data(:, 1);          % 腿长数据 (m)
alpha_rad = spring_data(:, 2);  % 夹角数据 (rad)

%% 3阶多项式拟合
fprintf("\n=== 3阶多项式拟合 ===\n");
degree = 3;
p = polyfit(L, alpha_rad, degree);

% 计算拟合值
L_fit = linspace(min(L), max(L), 200);
alpha_fit_rad = polyval(p, L_fit);

% 计算拟合质量指标
alpha_pred_rad = polyval(p, L);
residuals = alpha_rad - alpha_pred_rad;
SS_res = sum(residuals.^2);
SS_tot = sum((alpha_rad - mean(alpha_rad)).^2);
r2 = 1 - SS_res / SS_tot;
rmse = sqrt(mean(residuals.^2));

fprintf("\n拟合质量指标:\n");
fprintf("  R² (决定系数): %.6f\n", r2);
fprintf("  RMSE (均方根误差): %.6f rad\n", rmse);

%% 绘制拟合曲线和散点图
figure('Position', [100, 100, 800, 600]);
hold on;
grid on;
box on;

% 绘制拟合曲线
plot(L_fit, alpha_fit_rad, 'r-', 'LineWidth', 2.5, 'DisplayName', '3阶拟合曲线');

% 绘制原始数据散点
scatter(L, alpha_rad, 80, 'b', 'filled', 'DisplayName', '原始数据');

title('气弹簧夹角α随腿长L变化关系', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('腿长 L (m)', 'FontSize', 12);
ylabel('夹角 α (rad)', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 11);

% 显示R²值
text(0.05, 0.88, sprintf('R² = %.6f\nRMSE = %.6f rad', r2, rmse), ...
    'Units', 'normalized', 'FontSize', 10, ...
    'BackgroundColor', [1, 1, 1, 0.9], 'EdgeColor', 'black');

%% 绘制气弹簧分力 Fx 和 Fy 曲线
N = 250;            % 气弹簧推力 (N)

% 根据拟合夹角计算分力
alpha_fit_forTpFd = polyval(p, L_fit);   % 拟合得到的连续夹角
Fx_fit = N * sin(alpha_fit_forTpFd);     % 切向分力
Fy_fit = N * cos(alpha_fit_forTpFd);     % 法向分力

% 同时计算原始数据点对应的分力（用于散点验证）
alpha_data = alpha_rad;                  % 原始夹角
Fx_data = N * sin(alpha_data);
Fy_data = N * cos(alpha_data);

% 新建分力图窗
figure('Position', [150, 100, 900, 500]);
subplot(1,2,1);
hold on; grid on; box on;
plot(L_fit, Fx_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fx (拟合)');
scatter(L, Fx_data, 60, 'b', 'filled', 'DisplayName', '原始点');
xlabel('腿长 L (m)'); ylabel('切向分力 Fx (N)');
title('Fx = N0·sin(α)'); legend('Location', 'best');

subplot(1,2,2);
hold on; grid on; box on;
plot(L_fit, Fy_fit, 'b-', 'LineWidth', 2, 'DisplayName', 'Fy (拟合)');
scatter(L, Fy_data, 60, 'r', 'filled', 'DisplayName', '原始点');
xlabel('腿长 L (m)'); ylabel('法向分力 Fy (N)');
title('Fy = N0·cos(α)'); legend('Location', 'best');
sgtitle(sprintf('气弹簧分力随腿长变化 (推力 N0 = %d N)', N), 'FontSize', 14, 'FontWeight', 'bold');

fprintf("\n分力范围: Fx [%.2f, %.2f] N, Fy [%.2f, %.2f] N\n", ...
    min(Fx_fit), max(Fx_fit), min(Fy_fit), max(Fy_fit));

%% ========== 等效关节力矩 tau1 / tau4 计算 ==========
% 预分配数组
tau1_data = zeros(size(L));
tau4_data = zeros(size(L));

fprintf("\n正在计算等效关节力矩...\n");
L1 = robot_stats.length_thigh;
L2 = robot_stats.length_calf;
L3 = robot_stats.length_calf;
L4 = robot_stats.length_thigh;
L5 = 0;
for i = 1:numel(L_fit)
    L0 = L_fit(i);  % 当前腿长

    % 正向运动学：求解关节角度（取phi1大的无机械干涉，即第二组解）
    all_solutions = VMC_position_forward_kinematics(L0, 0, L1, L2, L3, L4, L5);
    phi1 = all_solutions(2, 1);   % 髋关节角度
    phi4 = all_solutions(2, 4);   % 踝关节角度

    % 力矩前向动力学：计算等效关节力矩（Fx为不需关注的内力）
    [tau1_data(i), tau4_data(i)] = VMC_torque_forward_kinematics(...
        0, Fy_fit(i), phi1, phi4, L1, L2, L3, L4, L5);
end
fprintf("力矩计算完成。\n");

% 绘制 tau1 / tau4 曲线
figure('Position', [200, 100, 800, 500]);
hold on; grid on; box on;

% tau1 (髋关节) —— 蓝色方块线
plot(L_fit, tau1_data, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'b', 'DisplayName', '\tau_1 (髋关节1)');

% tau4 (踝关节) —— 红色圆点线
plot(L_fit, tau4_data, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'r', 'DisplayName', '\tau_4 (髋关节4)');

xlabel('腿长 L (m)', 'FontSize', 12);
ylabel('等效关节力矩 (N·m)', 'FontSize', 12);
title('气弹簧推力引起的关节等效力矩', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;

% 添加零线参考（若力矩有正负变化）
yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

% 输出力矩范围
fprintf("\n关节力矩范围:\n");
fprintf("  tau1: [%.3f, %.3f] N·m\n", min(tau1_data), max(tau1_data));
fprintf("  tau4: [%.3f, %.3f] N·m\n", min(tau4_data), max(tau4_data));

%% 将拟合结果转换为浮点数表达式
prec = 6;   % 精度
expr = poly2cstring(p, 'l', prec);
fprintf('float alpha = %s\n', expr);

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