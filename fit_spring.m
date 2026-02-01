% 脚本：气弹簧夹角随腿长变化的3阶多项式拟合
tic
clear
clc
fprintf("\n=== 气弹簧夹角随腿长变化分析 ===\n");

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

% 显示拟合方程和R²值
eq_str = '拟合方程: α = ';
for i = 1:length(p)
    power = length(p) - i;
    if power > 1
        eq_str = [eq_str, sprintf('%.4fL^{%d} + ', p(i), power)];
    elseif power == 1
        eq_str = [eq_str, sprintf('%.4fL + ', p(i))];
    else
        eq_str = [eq_str, sprintf('%.4f', p(i))];
    end
end
text(0.05, 0.95, eq_str, 'Units', 'normalized', 'FontSize', 10, ...
    'BackgroundColor', [1, 1, 1, 0.9], 'EdgeColor', 'black', 'FontWeight', 'bold');

text(0.05, 0.88, sprintf('R² = %.6f\nRMSE = %.6f rad', r2, rmse), ...
    'Units', 'normalized', 'FontSize', 10, ...
    'BackgroundColor', [1, 1, 1, 0.9], 'EdgeColor', 'black');

%% 生成拟合函数
fprintf("\n=== 拟合函数 ===\n");
fprintf("使用以下函数计算任意腿长L对应的夹角α(rad):\n\n");

% 显示函数代码
function_code = ['function alpha = spring_angle(L)\n' ...
    '    % 计算气弹簧与腿长方向的夹角(rad)\n' ...
    '    % 输入: L - 腿长(m)\n' ...
    '    % 输出: alpha - 夹角(rad)\n' ...
    '    \n' ...
    '    % 3阶多项式拟合系数\n'];
for i = 1:length(p)
    function_code = [function_code sprintf('    p(%d) = %.8f; %% L^%d系数\n', i, p(i), length(p)-i)];
end
function_code = [function_code '    \n' ...
    '    alpha = polyval(p, L);\n' ...
    'end\n'];

fprintf(function_code);

%% 保存结果
save('spring_fit_result.mat', 'L', 'alpha_rad', 'p', 'r2', 'rmse');
fprintf('\n结果已保存到 spring_fit_result.mat\n');

%% 显示执行时间
toc
fprintf("\n=== 分析完成 ===\n");