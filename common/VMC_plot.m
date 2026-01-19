function [L0, phi0, phi2, phi3] = VMC_plot(phi1, phi4, L1, L2, L3, L4, L5, show_plot)
    % 绘制并联机构的2D关节连杆关系图
    % 坐标系：原点在A点，x轴向右为正，y轴向上为正
    % 所有角度从x轴正方向开始，逆时针增长
    
    % 计算正运动学
    [L0, phi0, phi2, phi3] = VMC_position_inverse_kinematics(phi1, phi4, L1, L2, L3, L4, L5);
    
    if show_plot == false
        return
    end
    
    % 计算各关键点坐标
    sin1 = sin(phi1);
    cos1 = cos(phi1);
    sin2 = sin(phi2);
    cos2 = cos(phi2);
    sin3 = sin(phi3);
    cos3 = cos(phi3);
    sin4 = sin(phi4);
    cos4 = cos(phi4);
    
    % 固定点
    A = [0; 0];  % 左边固定点（原点）
    E = [L5; 0]; % 右边固定点
    
    % 各关节点
    B = A + L1 * [cos1; sin1];  % 左摆杆末端
    D = E + L4 * [cos4; sin4];  % 右摆杆末端
    C = B + L2 * [cos2; sin2];  % 中心点
    F = C + L0 * [cos(phi0); sin(phi0)];  % 输出杆末端
    
    % 定义各角度颜色
    phi1_color = [0.8, 0.2, 0.2];    % 红色系
    phi2_color = [0.2, 0.6, 0.2];    % 绿色系
    phi3_color = [0.2, 0.2, 0.8];    % 蓝色系
    phi4_color = [0.6, 0.2, 0.6];    % 紫色系
    phi0_color = [0.9, 0.5, 0.1];    % 橙色系（输出角度）
    
    % 创建图形
    figure('Position', [100, 100, 1400, 600]);
    
    % 子图1：机构示意图
    subplot(1, 3, [1, 2]);
    hold on;
    grid on;
    box on;
    axis equal;
    title('并联机构示意图', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('X (mm 或 m)', 'FontSize', 12);
    ylabel('Y (mm 或 m)', 'FontSize', 12);
    
    % 绘制连杆（灰色，减少视觉干扰）
    plot([A(1), B(1)], [A(2), B(2)], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2.5);  % 连杆 L1
    plot([B(1), C(1)], [B(2), C(2)], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2.5);  % 连杆 L2
    plot([C(1), D(1)], [C(2), D(2)], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2.5);  % 连杆 L3
    plot([D(1), E(1)], [D(2), E(2)], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2.5);  % 连杆 L4
    
    % 绘制输出杆（黑色，粗一些）
    plot([C(1), F(1)], [C(2), F(2)], 'k-', 'LineWidth', 3.5);
    
    % 绘制固定点
    plot(A(1), A(2), 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    plot(E(1), E(2), 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    
    % 绘制活动关节点
    plot(B(1), B(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', phi1_color, 'MarkerEdgeColor', 'k');
    plot(D(1), D(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', phi4_color, 'MarkerEdgeColor', 'k');
    plot(C(1), C(2), 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0.3, 0.3, 0.3], 'MarkerEdgeColor', 'k');
    plot(F(1), F(2), 'o', 'MarkerSize', 14, 'MarkerFaceColor', phi0_color, 'MarkerEdgeColor', 'k');
    
    % 添加标签
    % text(A(1), A(2)-0.15, 'A', 'FontSize', 11, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    % text(E(1), E(2)-0.15, 'E', 'FontSize', 11, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    % text(B(1), B(2)+0.15, 'B', 'FontSize', 11, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    % text(D(1), D(2)+0.15, 'D', 'FontSize', 11, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    % text(C(1), C(2)+0.15, 'C', 'FontSize', 11, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    % text(F(1), F(2)+0.15, 'F', 'FontSize', 11, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    
    % 添加连杆长度标注
    label_offset = 0.08;  % 标注偏移距离
    font_size = 11;
    
    % 标注L1
    mid_AB = (A + B) / 2;
    text(mid_AB(1), mid_AB(2), sprintf('L1=%.3f', L1), ...
         'FontSize', font_size, 'FontWeight', 'bold', 'Color', 'b', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'BackgroundColor', [1, 1, 1, 0.7], 'EdgeColor', 'b');
    
    % 标注L2
    mid_BC = (B + C) / 2;
    text(mid_BC(1), mid_BC(2), sprintf('L2=%.3f', L2), ...
         'FontSize', font_size, 'FontWeight', 'bold', 'Color', 'b', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'BackgroundColor', [1, 1, 1, 0.7], 'EdgeColor', 'b');
    
    % 标注L3
    mid_CD = (C + D) / 2;
    text(mid_CD(1), mid_CD(2), sprintf('L3=%.3f', L3), ...
         'FontSize', font_size, 'FontWeight', 'bold', 'Color', 'b', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'BackgroundColor', [1, 1, 1, 0.7], 'EdgeColor', 'b');
    
    % 标注L4
    mid_DE = (D + E) / 2;
    text(mid_DE(1), mid_DE(2), sprintf('L4=%.3f', L4), ...
         'FontSize', font_size, 'FontWeight', 'bold', 'Color', 'b', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'BackgroundColor', [1, 1, 1, 0.7], 'EdgeColor', 'b');
    
    % 标注L0
    mid_CF = (C + F) / 2;
    text(mid_CF(1), mid_CF(2), sprintf('L0=%.3f', L0), ...
         'FontSize', font_size, 'FontWeight', 'bold', 'Color', 'r', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'BackgroundColor', [1, 1, 1, 0.7], 'EdgeColor', 'r');
    
    % 添加角度标记（从x轴正方向开始），使用不同的颜色
    add_angle_marker(A, [A(1)+0.5, A(2)], B, 0.2, phi1, '\phi_1', phi1_color);
    add_angle_marker(B, [B(1)+0.5, B(2)], C, 0.2, phi2, '\phi_2', phi2_color);
    add_angle_marker(D, [D(1)+0.5, D(2)], C, 0.2, phi3, '\phi_3', phi3_color);
    add_angle_marker(E, [E(1)+0.5, E(2)], D, 0.2, phi4, '\phi_4', phi4_color);
    add_angle_marker(C, [C(1)+0.5, C(2)], F, 0.2, phi0, '\phi_0', phi0_color);
    
    % 绘制坐标轴
    x_min = -max([L1, L2, L3, L4, L0]) * 0.3;
    x_max = max([L1+L2, L5+L4+L3, F(1)]) * 1.1;
    y_min = min([0, B(2), C(2), D(2), F(2)]) * 1.3;
    y_max = max([0, B(2), C(2), D(2), F(2)]) * 1.3;
    
    % 绘制x轴和y轴
    plot([x_min, x_max], [0, 0], 'k--', 'LineWidth', 0.8); % x轴
    plot([0, 0], [y_min, y_max], 'k--', 'LineWidth', 0.8); % y轴
    
    % 添加坐标系指示
    arrow_length = (x_max - x_min) * 0.06;
    quiver(0, 0, arrow_length, 0, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 2);
    quiver(0, 0, 0, arrow_length, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 2);
    text(arrow_length*1.2, -0.02, 'X', 'FontSize', 12, 'FontWeight', 'bold');
    text(-0.02, arrow_length*1.2, 'Y', 'FontSize', 12, 'FontWeight', 'bold');
    
    % 设置坐标范围
    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    
    hold off;
    
    % 子图3：机构参数显示（密排）
    subplot(1, 3, 3);
    axis off;
    title('机构参数与计算结果', 'FontSize', 14, 'FontWeight', 'bold');
    
    % 文本位置参数
    text_x = 0.05;
    text_y = 0.98;
    line_height = 0.035;  % 行高
    section_gap = 0.02;   % 节间间隙
    
    % 1. 机构参数
    text(text_x, text_y, '机构参数:', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0, 0, 0.7]);
    text_y = text_y - line_height;
    
    params = {
        sprintf('L1 = %.4f (大腿)', L1);
        sprintf('L2 = %.4f (小腿)', L2);
        sprintf('L3 = %.4f (小腿)', L3);
        sprintf('L4 = %.4f (大腿)', L4);
        sprintf('L5 = %.4f (固定间距)', L5);
    };
    
    for i = 1:length(params)
        text(text_x, text_y - (i-1)*line_height, params{i}, 'FontSize', 10);
    end
    text_y = text_y - (length(params))*line_height - section_gap;
    
    % 2. 输入角度
    text(text_x, text_y, '输入角度:', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0, 0, 0.7]);
    text_y = text_y - line_height;
    
    inputs = {
        sprintf('φ₁ = %.4f rad  (%.2f°)', phi1, rad2deg(phi1));
        sprintf('φ₄ = %.4f rad  (%.2f°)', phi4, rad2deg(phi4));
    };
    
    % 使用对应颜色绘制输入角度
    text(text_x, text_y - (0)*line_height, inputs{1}, 'FontSize', 10, 'Color', phi1_color, 'FontWeight', 'bold');
    text(text_x, text_y - (1)*line_height, inputs{2}, 'FontSize', 10, 'Color', phi4_color, 'FontWeight', 'bold');
    text_y = text_y - (length(inputs))*line_height - section_gap;
    
    % 3. 计算结果
    text(text_x, text_y, '计算结果:', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0, 0, 0.7]);
    text_y = text_y - line_height;
    
    calculations = {
        sprintf('φ₂ = %.4f rad  (%.2f°)', phi2, rad2deg(phi2));
        sprintf('φ₃ = %.4f rad  (%.2f°)', phi3, rad2deg(phi3));
        sprintf('φ₀ = %.4f rad  (%.2f°)', phi0, rad2deg(phi0));
        sprintf('L₀ = %.4f', L0);
    };
    
    % 使用对应颜色绘制计算出的角度
    text(text_x, text_y - (0)*line_height, calculations{1}, 'FontSize', 10, 'Color', phi2_color, 'FontWeight', 'bold');
    text(text_x, text_y - (1)*line_height, calculations{2}, 'FontSize', 10, 'Color', phi3_color, 'FontWeight', 'bold');
    text(text_x, text_y - (2)*line_height, calculations{3}, 'FontSize', 10, 'Color', phi0_color, 'FontWeight', 'bold');
    text(text_x, text_y - (3)*line_height, calculations{4}, 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold'); % L0用黑色
    text_y = text_y - (length(calculations))*line_height - section_gap;
    
    % 4. 坐标信息
    text(text_x, text_y, '关键点坐标:', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0, 0, 0.7]);
    text_y = text_y - line_height;
    
    coords = {
        sprintf('A: (%.4f, %.4f)', A(1), A(2));
        sprintf('B: (%.4f, %.4f)', B(1), B(2));
        sprintf('C: (%.4f, %.4f)', C(1), C(2));
        sprintf('D: (%.4f, %.4f)', D(1), D(2));
        sprintf('E: (%.4f, %.4f)', E(1), E(2));
        sprintf('F: (%.4f, %.4f)', F(1), F(2));
    };
    
    for i = 1:length(coords)
        text(text_x, text_y - (i-1)*line_height, coords{i}, 'FontSize', 9, 'Color', [0.3, 0.3, 0.3]);
    end
    text_y = text_y - (length(coords))*line_height - section_gap;
    
    % 5. 几何验证
    text(text_x, text_y, '几何验证:', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0, 0, 0.7]);
    text_y = text_y - line_height;
    
    % 验证C点从B和D计算的两种方式
    C_from_D = D + L3 * [cos3; sin3];
    C_error = norm(C - C_from_D);
    BD = norm(D - B);
    
    % 余弦定理验证
    cos_beta = (L2^2 + BD^2 - L3^2) / (2 * L2 * BD);
    if abs(cos_beta) <= 1
        beta = acos(cos_beta);
        alpha = atan2(D(2)-B(2), D(1)-B(1));
        phi2_calc = beta + alpha;
        phi3_calc = pi - beta + alpha;
        phi2_error = abs(phi2 - phi2_calc);
        phi3_error = abs(phi3 - phi3_calc);
        
        verifications = {
            sprintf('C点误差: %.6f', C_error);
            sprintf('φ₂误差: %.6f rad', phi2_error);
            sprintf('φ₃误差: %.6f rad', phi3_error);
            sprintf('BD距离: %.4f', BD);
        };
    else
        verifications = {
            sprintf('C点误差: %.6f', C_error);
            sprintf('余弦超限: %.4f', cos_beta);
            sprintf('BD距离: %.4f', BD);
        };
    end
    
    for i = 1:length(verifications)
        text(text_x, text_y - (i-1)*line_height, verifications{i}, 'FontSize', 9);
    end
    
end
    
function add_angle_marker(vertex, start_point, end_point, radius, angle, label, color)
    % 添加角度标记
    % vertex: 角度顶点
    % start_point: 起始参考点（定义x轴方向）
    % end_point: 结束点（定义角度的结束边）
    % radius: 圆弧半径
    % angle: 角度值（弧度）
    % label: 角度标签
    % color: 角度颜色
    
    % 计算起始向量和结束向量
    start_vec = start_point(:) - vertex(:);
    end_vec = end_point(:) - vertex(:);
    
    % 计算起始角度和结束角度（从x轴正方向开始）
    start_angle = atan2(start_vec(2), start_vec(1));
    end_angle = atan2(end_vec(2), end_vec(1));
    
    % 确保角度在[0, 2π)范围内
    start_angle = mod(start_angle, 2*pi);
    end_angle = mod(end_angle, 2*pi);
    
    % 绘制圆弧
    theta = linspace(start_angle, end_angle, 30);
    x = vertex(1) + radius * cos(theta);
    y = vertex(2) + radius * sin(theta);
    plot(x, y, '-', 'LineWidth', 2, 'Color', color);
    
    % 绘制角度边（从顶点到圆弧起点和终点）
    plot([vertex(1), vertex(1) + radius*cos(start_angle)], ...
         [vertex(2), vertex(2) + radius*sin(start_angle)], ...
         '--', 'LineWidth', 1.2, 'Color', color);
    
    plot([vertex(1), vertex(1) + radius*cos(end_angle)], ...
         [vertex(2), vertex(2) + radius*sin(end_angle)], ...
         '--', 'LineWidth', 1.2, 'Color', color);
    
    % 添加角度标签
    mid_angle = (start_angle + end_angle) / 2;
    label_radius = radius * 1.35;
    label_x = vertex(1) + label_radius * cos(mid_angle);
    label_y = vertex(2) + label_radius * sin(mid_angle);
    text(label_x, label_y, label, 'FontSize', 11, 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'Color', color, 'BackgroundColor', 'white');

end

function deg = rad2deg(rad)
    % 弧度转角度
    deg = rad * 180 / pi;
end