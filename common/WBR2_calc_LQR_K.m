function [K, OX] = WBR2_calc_LQR_K(Q, R, length_left, length_right, control_freq, STATS_FILE, MASS_MATRIX)
    % 计算指定腿长下的K矩阵（4 * 10），模型：WBR2
    % 因线性化的位置（平衡点）会随着腿长变化，故每次都要重算一遍AB矩阵
    dt = 1 / control_freq;

    %% 参数代换
    [sym_vars, sym_vals] = compute_state(length_left, length_right, STATS_FILE, MASS_MATRIX);
    M_fit = simplify(subs(MASS_MATRIX.M_sym, sym_vars, sym_vals));
    B_fit = simplify(subs(MASS_MATRIX.B_sym, sym_vars, sym_vals));
    g_fit = simplify(subs(MASS_MATRIX.g_sym, sym_vars, sym_vals));

    %% 求解平衡点
    % ========== 广义坐标及其导数 ==========
    syms X_b_h dX_b_h ddX_b_h real
    syms phi dphi ddphi real
    syms theta_b theta_l theta_r real
    syms dtheta_l ddtheta_l real
    syms dtheta_r ddtheta_r real
    syms dtheta_b ddtheta_b real
    % ========== 控制力矩 ==========
    syms T_r_to_b T_l_to_b T_wr_to_r T_wl_to_l real
    
    eq_subs_partial = {
        dtheta_l, 0;
        dtheta_r, 0;
        dtheta_b, 0;
        ddtheta_l, 0;
        ddtheta_r, 0;
        ddtheta_b, 0;
        ddX_b_h, 0;
        ddphi, 0;
        phi, 0;
        dphi, 0;
        X_b_h, 0;
        dX_b_h, 0;
        T_r_to_b, 0;
        T_l_to_b, 0;
        T_wr_to_r, 0;
        T_wl_to_l, 0;
    };
    g_fit_eq = subs(g_fit, eq_subs_partial(:,1), eq_subs_partial(:,2));
    % 符号求解各 theta
    % theta_b 从 g(2)=0
    if g_fit_eq(2) == 0
        theta_b_fit = 0;
    else
        sol_b = solve(g_fit_eq(2) == 0, theta_b, 'Real', true);
        sol_b_vals = double(sol_b);
        [~, idx_b] = min(abs(sol_b_vals));
        theta_b_fit = sol_b_vals(idx_b);
    end
    
    % theta_r 从 g(3)=0
    if g_fit_eq(3) == 0
        theta_r_fit_val = 0;
    else
        sol_r = solve(g_fit_eq(3) == 0, theta_r, 'Real', true);
        sol_r_vals = double(sol_r);
        [~, idx_r] = min(abs(sol_r_vals));
        theta_r_fit_val = sol_r_vals(idx_r);
    end
    
    % theta_l 从 g(4)=0
    if g_fit_eq(4) == 0
        theta_l_fit_val = 0;
    else
        sol_l = solve(g_fit_eq(4) == 0, theta_l, 'Real', true);
        sol_l_vals = double(sol_l);
        [~, idx_l] = min(abs(sol_l_vals));
        theta_l_fit_val = sol_l_vals(idx_l);
    end
    
    theta_star_fit = [theta_l_fit_val; theta_r_fit_val; theta_b_fit];
    
    % 完整平衡点代换 (包括加速度和控制输入)
    eq_subs_fit = {
        theta_l, theta_star_fit(1); theta_r, theta_star_fit(2); theta_b, theta_star_fit(3);
        dtheta_l, 0; dtheta_r, 0; dtheta_b, 0;
        ddtheta_l, 0; ddtheta_r, 0; ddtheta_b, 0;
        ddX_b_h, 0; ddphi, 0;
        phi, 0; dphi, 0; X_b_h, 0; dX_b_h, 0;
        T_r_to_b, 0; T_l_to_b, 0; T_wr_to_r, 0; T_wl_to_l, 0;
    };
    
    M_eq_fit = double(subs(M_fit, eq_subs_fit(:,1), eq_subs_fit(:,2)));
    B_eq_fit = double(subs(B_fit, eq_subs_fit(:,1), eq_subs_fit(:,2)));
    
    dg_dtheta_l_fit = double(subs(diff(g_fit, theta_l), eq_subs_fit(:,1), eq_subs_fit(:,2)));
    dg_dtheta_r_fit = double(subs(diff(g_fit, theta_r), eq_subs_fit(:,1), eq_subs_fit(:,2)));
    dg_dtheta_b_fit = double(subs(diff(g_fit, theta_b), eq_subs_fit(:,1), eq_subs_fit(:,2)));
    
    dg_dtheta_fit = [zeros(5,1), dg_dtheta_l_fit, dg_dtheta_r_fit, dg_dtheta_b_fit];
    A_dyn_fit = M_eq_fit \ dg_dtheta_fit;
    
    A_fit_num = zeros(10, 10);
    B_fit_num = zeros(10, 4);
    
    A_fit_num(1,2) = 1; A_fit_num(3,4) = 1; A_fit_num(5,6) = 1;
    A_fit_num(7,8) = 1; A_fit_num(9,10) = 1;
    
    A_fit_num(2,5) = A_dyn_fit(1,2); A_fit_num(2,7) = A_dyn_fit(1,3); A_fit_num(2,9) = A_dyn_fit(1,4);
    A_fit_num(4,5) = A_dyn_fit(2,2); A_fit_num(4,7) = A_dyn_fit(2,3); A_fit_num(4,9) = A_dyn_fit(2,4);
    A_fit_num(6,5) = A_dyn_fit(3,2); A_fit_num(6,7) = A_dyn_fit(3,3); A_fit_num(6,9) = A_dyn_fit(3,4);
    A_fit_num(8,5) = A_dyn_fit(4,2); A_fit_num(8,7) = A_dyn_fit(4,3); A_fit_num(8,9) = A_dyn_fit(4,4);
    A_fit_num(10,5) = A_dyn_fit(5,2); A_fit_num(10,7) = A_dyn_fit(5,3); A_fit_num(10,9) = A_dyn_fit(5,4);
    
    B_dyn_fit = M_eq_fit \ B_eq_fit;
    B_fit_num(2,:) = B_dyn_fit(1,:); B_fit_num(4,:) = B_dyn_fit(2,:);
    B_fit_num(6,:) = B_dyn_fit(3,:); B_fit_num(8,:) = B_dyn_fit(4,:);
    B_fit_num(10,:) = B_dyn_fit(5,:);

    [Ad, Bd] = c2d(A_fit_num, B_fit_num, dt);
    K = dlqr(Ad, Bd, Q, R, 0); % 状态和输入变量之间不存在耦合关系，不需要最后一项交叉加权矩阵
    OX = theta_star_fit'; % 3个偏移平衡点: 左腿摆角、右腿摆角、底盘俯仰角
end

function [sym_vars, sym_vals] = compute_state(length_left, length_right, STATS_FILE, STATE_SPACE_FILE)
    % 对照https://github.com/HKUAstar/Balance-Model-2026/blob/main/%E6%8E%A8%E5%AF%BC-v2.md的定义
    m_b = STATE_SPACE_FILE.m_b;
    m_l = STATE_SPACE_FILE.m_l;
    m_r = STATE_SPACE_FILE.m_r;
    m_wl = STATE_SPACE_FILE.m_wl;
    m_wr = STATE_SPACE_FILE.m_wr;
    I_b = STATE_SPACE_FILE.I_b;
    I_l = STATE_SPACE_FILE.I_l;
    I_r = STATE_SPACE_FILE.I_r;
    I_wl = STATE_SPACE_FILE.I_wl;
    I_wr = STATE_SPACE_FILE.I_wr;
    I_yaw = STATE_SPACE_FILE.I_yaw;
    l_l = STATE_SPACE_FILE.l_l;
    l_r = STATE_SPACE_FILE.l_r;
    l_l_d = STATE_SPACE_FILE.l_l_d;
    l_r_d = STATE_SPACE_FILE.l_r_d;
    l_b = STATE_SPACE_FILE.l_b;
    R = STATE_SPACE_FILE.R;
    R_w = STATE_SPACE_FILE.R_w;
    g = STATE_SPACE_FILE.g;
    theta_l0 = STATE_SPACE_FILE.theta_l0;
    theta_r0 = STATE_SPACE_FILE.theta_r0;
    theta_b0 = STATE_SPACE_FILE.theta_b0;
    sym_vars = [m_b, m_l, m_r, m_wl, m_wr, ...
                I_b, I_l, I_r, I_wl, I_wr, I_yaw, ...
                l_l, l_r, l_l_d, l_r_d, l_b, R, R_w, g, ...
                theta_l0, theta_r0, theta_b0];

    fit_distance_hip_leg_left = polyval(STATS_FILE.coeff_distance_hip_leg, length_left);
    fit_distance_hip_leg_right = polyval(STATS_FILE.coeff_distance_hip_leg, length_right);
    fit_distance_wheel_leg_left = polyval(STATS_FILE.coeff_distance_wheel_leg, length_left);
    fit_distance_wheel_leg_right = polyval(STATS_FILE.coeff_distance_wheel_leg, length_right);

    m_b_ac = STATS_FILE.mass_body;
    m_l_ac = STATS_FILE.mass_leg_aver;
    m_r_ac = STATS_FILE.mass_leg_aver;
    m_wl_ac = STATS_FILE.mass_wheel_aver;
    m_wr_ac = STATS_FILE.mass_wheel_aver;
    I_b_ac = STATS_FILE.inertia_body;
    I_l_ac = polyval(STATS_FILE.coeff_inertia_leg, length_left);
    I_r_ac = polyval(STATS_FILE.coeff_inertia_leg, length_right);
    I_wl_ac = STATS_FILE.inertia_wheel_aver;
    I_wr_ac = STATS_FILE.inertia_wheel_aver;
    I_yaw_ac = STATS_FILE.inertia_body;
    l_l_ac = length_left;
    l_r_ac = length_right;
    l_l_d_ac = fit_distance_hip_leg_left;
    l_r_d_ac = fit_distance_hip_leg_right;
    l_b_ac = STATS_FILE.distance_body_hip;
    R_ac = STATS_FILE.radius_wheel;
    R_w_ac = STATS_FILE.distance_half_wheel;
    g_ac = STATS_FILE.g_;
    theta_l0_ac = calc_leg_centroid_offset_ang(length_left, fit_distance_wheel_leg_left, l_l_d_ac);
    theta_r0_ac = calc_leg_centroid_offset_ang(length_right, fit_distance_wheel_leg_right, l_r_d_ac);
    theta_b0_ac = STATS_FILE.body_centroid_offset_ang;
    sym_vals = [m_b_ac, m_l_ac, m_r_ac, m_wl_ac, m_wr_ac, I_b_ac, ...
                I_l_ac, I_r_ac, I_wl_ac, I_wr_ac, I_yaw_ac, l_l_ac, ...
                l_r_ac, l_l_d_ac, l_r_d_ac, l_b_ac, R_ac, R_w_ac, g_ac, ...
                theta_l0_ac, theta_r0_ac, theta_b0_ac];
end