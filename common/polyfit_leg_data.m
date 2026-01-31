function [coeff_dis_w_l, coeff_dis_h_l, coeff_I_l] = polyfit_leg_data(leg_data, varargin)
% POLYFIT_LEG_DATA 拟合随腿长变化的机器人参数
% 
% 输入参数:
%   leg_data - 腿长数据矩阵，格式为 [腿长, dis_w_l, dis_h_l, I_l]
%   varargin   - 可选参数对:
%       'DegreeDisWl' - dis_w_l的拟合阶数，默认3
%       'DegreeDisHl' - dis_b_l的拟合阶数，默认3 
%       'DegreeIl'    - I_l的拟合阶数，默认3
%       'ShowPlot'    - 是否显示图形，默认true
%       'SaveCoeff'   - 是否保存系数到文件，默认false
%       'FileName'    - 保存的文件名，默认'leg_fit_coeff.mat'
%
% 输出参数:
%   coeff_dis_w_l - dis_w_l的拟合系数(从高次到低次)
%   coeff_dis_h_l - dis_b_l的拟合系数(从高次到低次)
%   coeff_I_l     - I_l的拟合系数(从高次到低次)
%
% 示例:
%   [coeff1, coeff2, coeff3] = polyfit_leg_data(leg_data);
%   [coeff1, coeff2, coeff3] = polyfit_leg_data(leg_data, 'ShowPlot', false);

% 解析可选参数
p = inputParser;
addRequired(p, 'leg_data', @(x) validateattributes(x, {'numeric'}, {'2d', 'ncols', 4}));
addParameter(p, 'DegreeDisWl', 3, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(p, 'DegreeDisHl', 3, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(p, 'DegreeIl', 3, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(p, 'ShowPlot', true, @islogical);
addParameter(p, 'SaveCoeff', false, @islogical);
addParameter(p, 'FileName', 'leg_fit_coeff.mat', @ischar);
parse(p, leg_data, varargin{:});

% 提取数据
L = leg_data(:, 1);        % 腿长
dis_w_l = leg_data(:, 2);  % 驱动轮质心到腿摆杆质心距离
dis_h_l = leg_data(:, 3);  % 机体转轴到腿摆杆质心距离
I_l = leg_data(:, 4);      % 腿摆杆绕髋关节转动惯量

% 检查数据有效性
if any(L <= 0)
    error('腿长数据必须为正数');
end
if size(leg_data, 1) < max([p.Results.DegreeDisWl, p.Results.DegreeDisHl, p.Results.DegreeIl]) + 1
    error('数据点数量不足以进行指定阶数的多项式拟合');
end

coeff_dis_w_l = polyfit(L, dis_w_l, p.Results.DegreeDisWl);
fit_dis_w_l = @(x) polyval(coeff_dis_w_l, x);

coeff_dis_h_l = polyfit(L, dis_h_l, p.Results.DegreeDisHl);
fit_dis_h_l = @(x) polyval(coeff_dis_h_l, x);

coeff_I_l = polyfit(L, I_l, p.Results.DegreeIl);
fit_I_l = @(x) polyval(coeff_I_l, x);

if p.Results.SaveCoeff
    save(p.Results.FileName, 'coeff_dis_w_l', 'coeff_dis_h_l', 'coeff_I_l', ...
         'DegreeDisWl', 'DegreeDisHl', 'DegreeIl');
    fprintf('拟合系数已保存到 %s\n', p.Results.FileName);
end

if p.Results.ShowPlot
    display_fit_results(L, dis_w_l, dis_h_l, I_l, ...
                       coeff_dis_w_l, coeff_dis_h_l, coeff_I_l, ...
                       p.Results.DegreeDisWl, p.Results.DegreeDisHl, p.Results.DegreeIl);
end

end

%% 辅助函数：显示拟合结果
function display_fit_results(L, dis_w_l, dis_h_l, I_l, ...
                            coeff_dis_w_l, coeff_dis_h_l, coeff_I_l, ...
                            degree1, degree2, degree3)
% 生成更密集的点用于绘制拟合曲线
L_fine = linspace(min(L), max(L), 200);

% 计算拟合值
fit_dis_w_l = @(x) polyval(coeff_dis_w_l, x);
fit_dis_h_l = @(x) polyval(coeff_dis_h_l, x);
fit_I_l = @(x) polyval(coeff_I_l, x);

% 计算拟合误差
dis_w_l_pred = fit_dis_w_l(L);
dis_w_l_rmse = sqrt(mean((dis_w_l_pred - dis_w_l).^2));

dis_h_l_pred = fit_dis_h_l(L);
dis_h_l_rmse = sqrt(mean((dis_h_l_pred - dis_h_l).^2));

I_l_pred = fit_I_l(L);
I_l_rmse = sqrt(mean((I_l_pred - I_l).^2));

% 创建图形窗口
figure('Position', [100, 100, 900, 600], 'Name', '腿长相关参数拟合效果');

% 子图1: dis_w_l 拟合
subplot(2, 2, 1);
plot(L, dis_w_l, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', '原始数据');
hold on;
plot(L_fine, fit_dis_w_l(L_fine), 'r-', 'LineWidth', 2, 'DisplayName', '拟合曲线');
xlabel('腿长 L (m)', 'FontSize', 12);
ylabel('dis\_w\_l (m)', 'FontSize', 12);
title(sprintf('dis\\_w\\_l拟合 (阶数: %d)', degree1), 'FontSize', 14);
legend('Location', 'best');
grid on;
hold off;

% 子图2: dis_h_l 拟合
subplot(2, 2, 2);
plot(L, dis_h_l, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', '原始数据');
hold on;
plot(L_fine, fit_dis_h_l(L_fine), 'r-', 'LineWidth', 2, 'DisplayName', '拟合曲线');
xlabel('腿长 L (m)', 'FontSize', 12);
ylabel('dis\_h\_l (m)', 'FontSize', 12);
title(sprintf('dis\\_h\\_l拟合 (阶数: %d)', degree2), 'FontSize', 14);
legend('Location', 'best');
grid on;
hold off;

% 子图3: I_l 拟合
subplot(2, 2, 3);
plot(L, I_l, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', '原始数据');
hold on;
plot(L_fine, fit_I_l(L_fine), 'r-', 'LineWidth', 2, 'DisplayName', '拟合曲线');
xlabel('腿长 L (m)', 'FontSize', 12);
ylabel('I\_l (kg·m²)', 'FontSize', 12);
title(sprintf('I\\_l拟合 (阶数: %d)', degree3), 'FontSize', 14);
legend('Location', 'best');
grid on;
hold off;

% 子图4: 拟合误差对比
subplot(2, 2, 4);
bar_data = [dis_w_l_rmse, dis_h_l_rmse, I_l_rmse*1000]; % I_l误差放大1000倍显示
bar(bar_data);
set(gca, 'XTickLabel', {'dis\_w\_l', 'dis\_b\_l', 'I\_l(×10^{-3})'});
ylabel('均方根误差 (RMSE)', 'FontSize', 12);
title('各参数拟合误差对比', 'FontSize', 14);

% 在柱状图上添加数值标签
for i = 1:length(bar_data)
    if i == 3
        text(i, bar_data(i), sprintf('%.6f', I_l_rmse), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    else
        text(i, bar_data(i), sprintf('%.6f', bar_data(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end
grid on;

sgtitle('腿长相关参数多项式拟合效果', 'FontSize', 16, 'FontWeight', 'bold');

fprintf('\n=== 拟合误差统计 ===\n');
fprintf('dis_w_l - 最大误差: %.6f m, RMSE: %.6f m\n', ...
        max(abs(dis_w_l_pred - dis_w_l)), dis_w_l_rmse);
fprintf('dis_h_l - 最大误差: %.6f m, RMSE: %.6f m\n', ...
        max(abs(dis_h_l_pred - dis_h_l)), dis_h_l_rmse);
fprintf('I_l     - 最大误差: %.6f kg·m², RMSE: %.6f kg·m²\n', ...
        max(abs(I_l_pred - I_l)), I_l_rmse);

fprintf('\n=== 拟合函数表达式 ===\n');
display_poly_expression('dis_w_l(L)', coeff_dis_w_l);
display_poly_expression('dis_h_l(L)', coeff_dis_h_l);
display_poly_expression('I_l(L)', coeff_I_l);

end

%% 辅助函数：显示多项式表达式
function display_poly_expression(name, coeff)
fprintf('\n%s = ', name);
terms = length(coeff);
for i = 1:terms
    power = terms - i;
    if power > 1
        term_str = sprintf('%.6f·L^%d', coeff(i), power);
    elseif power == 1
        term_str = sprintf('%.6f·L', coeff(i));
    else % power == 0
        term_str = sprintf('%.6f', coeff(i));
    end
    
    if i == 1
        fprintf('%s', term_str);
    else
        if coeff(i) >= 0
            fprintf(' + %s', term_str);
        else
            fprintf(' - %s', term_str(2:end)); % 去掉负号
        end
    end
end
fprintf('\n');
end