function WBR2_generate_kinematics()
    % 10输入4输出状态空间方程AB矩阵求解，在拟合的平衡点处线性化，而非在0点线性化
    % 脚本会解出A、B矩阵的表达式sym_A和sym_B
    % 写法参考港科开源

    fprintf("\n=== 正在运行" + mfilename + ".m ===\n");
    %% 推导轮腿机器人的运动方程，并进行简化
    % 广义坐标及其导数
    syms ddX_b_h real

    % 物理常数
    syms g real

    % 机体参数
    syms m_b I_b l_b theta_b0 real
    syms theta_b ddtheta_b real      % 机体角度、角速度、角加速度
    syms T_l_to_b real  % 左腿对机体的力和力矩
    syms T_r_to_b real  % 右腿对机体的力和力矩
    syms a_b_h a_b_v real                     % 机体加速度

    % 腿部参数
    syms m_l m_r I_l I_r l_l l_r l_l_d l_r_d theta_l0 theta_r0 real
    syms theta_l dtheta_l ddtheta_l real
    syms theta_r dtheta_r ddtheta_r real
    syms T_wl_to_l real  % 左轮对左腿
    syms T_wr_to_r real  % 右轮对右腿
    syms a_l_h a_l_v a_r_h a_r_v real            % 腿加速度

    % 轮子参数
    syms m_wl m_wr I_wl I_wr real
    syms ddtheta_wl ddtheta_wr real
    syms a_wl_h a_wl_v a_wr_h a_wr_v real  % 轮加速度

    % 几何参数
    syms R R_w I_yaw real

    % Yaw角
    syms ddphi real

    % 最终方程：消去内力，得到五个广义坐标的二阶导数表达式
    % 推导参考 https://github.com/HKUAstar/Balance-Model-2026.git

    % 最终方程1: 整体水平动量
    eq1 = m_b*a_b_h + m_l*a_l_h + m_r*a_r_h + m_wl*a_wl_h + m_wr*a_wr_h ...
        + (T_wl_to_l + T_wr_to_r + I_wl*ddtheta_wl + I_wr*ddtheta_wr)/R;

    % 最终方程2: 机体转动
    eq2 = I_b*ddtheta_b - T_l_to_b - T_r_to_b - m_b*g*l_b*sin(theta_b + theta_b0);

    % 最终方程3: 右腿转动方程
    eq3 = I_r*ddtheta_r - T_wr_to_r + T_r_to_b - m_r*g*l_r_d*sin(theta_r + theta_r0) ...
        - m_wr*a_wr_h*l_r*cos(theta_r) - (T_wr_to_r + I_wr*ddtheta_wr)*l_r*cos(theta_r)/R ...
        - (m_b*a_b_v + m_l*a_l_v + m_r*a_r_v)*l_r*sin(theta_r)/2 ...
        - (m_wr - m_b - m_l - m_r - m_wl)*g*l_r*sin(theta_r)/2;

    % 最终方程4: 左腿转动方程
    eq4 = I_l*ddtheta_l - T_wl_to_l + T_l_to_b - m_l*g*l_l_d*sin(theta_l + theta_l0) ...
        - m_wl*a_wl_h*l_l*cos(theta_l) - (T_wl_to_l + I_wl*ddtheta_wl)*l_l*cos(theta_l)/R ...
        - (m_b*a_b_v + m_l*a_l_v + m_r*a_r_v)*l_l*sin(theta_l)/2 ...
        - (m_wl - m_b - m_l - m_r - m_wr)*g*l_l*sin(theta_l)/2;

    % 最终方程5: Yaw转动方程
    eq5 = I_yaw*ddphi - R_w/R*(T_wl_to_l + I_wl*ddtheta_wl - T_wr_to_r - I_wr*ddtheta_wr);

    fprintf('得到 5 个最终方程（五个广义坐标的二阶导数表达式）\n\n');
    % %% 定义符号变量
    % % ========== 物理参数 ==========
    % syms real                     % 重力加速度
    % syms R R_w real               % 轮半径, 半轮距
    % syms l_l l_r real             % 长度参数
    % % ========== 广义坐标及其导数 (推导 §5.0) ==========
    % % q = [X_b^h, phi, theta_l, theta_r, theta_b]
    % syms ddX_b_h real     % 机体水平位置
    % syms ddphi real           % yaw角
    % syms theta_l dtheta_l ddtheta_l real
    % syms theta_r dtheta_r ddtheta_r real
    % % ========== 轮角速度 ==========
    % syms ddtheta_wl real
    % syms ddtheta_wr real
    % % ========== 加速度变量 (原方程中出现的) ==========
    % syms a_b_h a_b_v real              % 机体加速度
    % syms a_l_h a_l_v a_r_h a_r_v real  % 腿转轴加速度
    % syms a_wl_h a_wl_v a_wr_h a_wr_v real  % 轮加速度

    % %% 方程用广义坐标表示
    % leg_accel_term = (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2;
    % leg_vel_sq_term = (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2;
    % leg_term = leg_accel_term - leg_vel_sq_term;

    % ddtheta_wr_expr = (ddX_b_h + R_w*ddphi)/R - leg_term/R;
    % ddtheta_wl_expr = (ddX_b_h - R_w*ddphi)/R - leg_term/R;

    %% 代换
    kinematics_subs = {
        % 轮加速度 (纯滚动)
        a_wr_h, R*ddtheta_wr;
        a_wl_h, R*ddtheta_wl;
        a_wr_v, sym(0);
        a_wl_v, sym(0);
        
        % 腿转轴加速度
        a_r_h, R*ddtheta_wr + l_r*cos(theta_r)*ddtheta_r - l_r*sin(theta_r)*dtheta_r^2;
        a_r_v, -l_r*sin(theta_r)*ddtheta_r - l_r*cos(theta_r)*dtheta_r^2;
        a_l_h, R*ddtheta_wl + l_l*cos(theta_l)*ddtheta_l - l_l*sin(theta_l)*dtheta_l^2;
        a_l_v, -l_l*sin(theta_l)*ddtheta_l - l_l*cos(theta_l)*dtheta_l^2;
        
        % 机体加速度 (左右腿转轴平均)
        a_b_h, (R*(ddtheta_wr + ddtheta_wl)/2 ...
            + (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2 ...
            - (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2);
        a_b_v, (-(l_r*sin(theta_r)*ddtheta_r + l_l*sin(theta_l)*ddtheta_l)/2 ...
            - (l_r*cos(theta_r)*dtheta_r^2 + l_l*cos(theta_l)*dtheta_l^2)/2);
    };
    eq1_new = subs(eq1, kinematics_subs(:,1), kinematics_subs(:,2));
    eq2_new = subs(eq2, kinematics_subs(:,1), kinematics_subs(:,2));
    eq3_new = subs(eq3, kinematics_subs(:,1), kinematics_subs(:,2));
    eq4_new = subs(eq4, kinematics_subs(:,1), kinematics_subs(:,2));
    eq5_new = subs(eq5, kinematics_subs(:,1), kinematics_subs(:,2));
    eq1_new = simplify(eq1_new);
    eq2_new = simplify(eq2_new);
    eq3_new = simplify(eq3_new);
    eq4_new = simplify(eq4_new);
    eq5_new = simplify(eq5_new);

    %% 轮角加速度代换为广义坐标
    leg_term = (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2 ...
         - (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2;
    ddtheta_wr_sub = (ddX_b_h + R_w*ddphi)/R - leg_term/R;
    ddtheta_wl_sub = (ddX_b_h - R_w*ddphi)/R - leg_term/R;
    wheel_subs = {ddtheta_wr, ddtheta_wr_sub; ddtheta_wl, ddtheta_wl_sub};

    eq1_new = simplify(subs(eq1_new, wheel_subs(:,1), wheel_subs(:,2)));
    eq2_new = simplify(subs(eq2_new, wheel_subs(:,1), wheel_subs(:,2)));
    eq3_new = simplify(subs(eq3_new, wheel_subs(:,1), wheel_subs(:,2)));
    eq4_new = simplify(subs(eq4_new, wheel_subs(:,1), wheel_subs(:,2)));
    eq5_new = simplify(subs(eq5_new, wheel_subs(:,1), wheel_subs(:,2)));

    %% 提取质量矩阵
    ddq = [ddX_b_h; ddphi; ddtheta_l; ddtheta_r; ddtheta_b];
    u = [T_r_to_b; T_l_to_b; T_wr_to_r; T_wl_to_l];
    eqs = {eq1_new, eq2_new, eq3_new, eq4_new, eq5_new}; 
    % M
    M_sym = sym(zeros(5,5));
    for i = 1:5
        for j = 1:5
            M_sym(i,j) = diff(eqs{i}, ddq(j));
        end
    end
    % B
    B_raw = sym(zeros(5,4));
    for i = 1:5
        for j = 1:4
            B_raw(i,j) = diff(eqs{i}, u(j));
        end
    end
    B_sym = -B_raw;
    % g
    g_sym = sym(zeros(5,1));
    for i = 1:5
        g_sym(i) = -(eqs{i} - M_sym(i,:)*ddq - B_raw(i,:)*u);
    end
    save WBR2_MASS.mat;
end
