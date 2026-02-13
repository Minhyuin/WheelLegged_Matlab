function [all_solutions] = VMC_position_forward_kinematics(L0, phi0, L1, L2, L3, L4, L5)
    % 寻找所有可能的正解
    %   注意，摆杆长度为正数，摆杆角度零点位置为摆杆竖直，详细见
    %   https://jjyczursn9.feishu.cn/wiki/IW2BwNm7GiU5gZkxX96c0rd5nxe?from=from_copylink
    
    % 首先计算C点坐标
    phiC = phi0 + pi/2;
    xC = L0 * cos(phiC);
    yC = L0 * sin(phiC);
    
    % 方法1: 几何法求解可能的phi1和phi2
    all_solutions = [];
    
    % 对于phi1的可能解 (通过求解B点到C点的距离约束)
    % B点坐标: xB = -L5/2 + L1*cos(phi1), yB = L1*sin(phi1)
    % 距离BC = L2
    
    % 重新整理方程: (xC - (-L5/2 + L1*cos(phi1)))^2 + (yC - L1*sin(phi1))^2 = L2^2
    % 令: dx = xC + L5/2, dy = yC
    % 方程: (dx - L1*cos(phi1))^2 + (dy - L1*sin(phi1))^2 = L2^2
    % 展开: dx^2 + dy^2 + L1^2 - 2*L1*(dx*cos(phi1) + dy*sin(phi1)) = L2^2
    % 令: R = sqrt(dx^2 + dy^2), theta = atan2(dy, dx)
    % 则: R*cos(phi1 - theta) = (dx^2 + dy^2 + L1^2 - L2^2)/(2*L1)
    
    dx = xC + L5/2;
    dy = yC;
    R = sqrt(dx^2 + dy^2);
    theta = atan2(dy, dx);
    A = (R^2 + L1^2 - L2^2) / (2 * L1 * R);
    
    if abs(A) <= 1
        delta = acos(A);
        % 两个可能的phi1解
        phi1_1 = theta + delta;
        phi1_2 = theta - delta;
        
        % 对于每个phi1，计算phi2
        for phi1 = [phi1_1, phi1_2]
            % 从phi1计算B点
            xB = -L5/2 + L1 * cos(phi1);
            yB = L1 * sin(phi1);
            
            % 从B和C计算phi2
            phi2 = atan2(yC - yB, xC - xB);
            
            % 现在从phi4和phi3计算
            % 同样方法求解phi4
            % D点坐标: xD = L5/2 + L4*cos(phi4), yD = L4*sin(phi4)
            % 距离DC = L3
            dx2 = xC - L5/2;
            dy2 = yC;
            R2 = sqrt(dx2^2 + dy2^2);
            theta2 = atan2(dy2, dx2);
            A2 = (R2^2 + L4^2 - L3^2) / (2 * L4 * R2);
            
            if abs(A2) <= 1
                delta2 = acos(A2);
                % 两个可能的phi4解
                phi4_1 = theta2 + delta2;
                phi4_2 = theta2 - delta2;
                
                for phi4 = [phi4_1, phi4_2]
                    % 从phi4计算D点
                    xD = L5/2 + L4 * cos(phi4);
                    yD = L4 * sin(phi4);
                    
                    % 计算phi3 (从D到C)
                    phi3 = atan2(yC - yD, xC - xD);
                    
                    % 验证约束: B, C, D构成的三角形
                    % 计算BD长度
                    xBD = xD - xB;
                    yBD = yD - yB;
                    lBD_calc = sqrt(xBD^2 + yBD^2);
                    
                    % 使用余弦定理验证
                    cosBeta = (L2^2 + lBD_calc^2 - L3^2) / (2 * L2 * lBD_calc);
                    
                    if abs(cosBeta) <= 1
                        % 这是一个有效解
                        all_solutions = [all_solutions; phi1, phi2, phi3, phi4];
                    end
                end
            end
        end
    end
    
    % 如果没有找到解，使用数值方法
    if isempty(all_solutions)
        % 使用多个初始值尝试
        initial_guesses = [0, 0, 0, 0;
                           pi/2, pi/2, pi/2, pi/2;
                           -pi/2, -pi/2, -pi/2, -pi/2;
                           pi, pi, pi, pi];
        
        options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
        
        for i = 1:size(initial_guesses, 1)
            try
                solution = fsolve(@(x) inverse_kinematics_equations(x, L0, phi0, L1, L2, L3, L4, L5), ...
                                  initial_guesses(i, :)', options);
                
                % 验证解的有效性
                phi1 = solution(1);
                phi2 = solution(2);
                phi3 = solution(3);
                phi4 = solution(4);
                
                % 计算验证误差
                [L0_calc, phi0_calc] = VMC_position_forward_kinematics(phi1, phi4, L1, L2, L3, L4, L5);
                
                if abs(L0_calc - L0) < 1e-6 && abs(mod(phi0_calc - phi0 + pi, 2*pi) - pi) < 1e-6
                    all_solutions = [all_solutions; phi1, phi2, phi3, phi4];
                end
            catch
                continue;
            end
        end
    end
    
    % 去除重复解（在一定容差内）
    if ~isempty(all_solutions)
        all_solutions = unique(round(all_solutions * 1e6) / 1e6, 'rows');
    end

end