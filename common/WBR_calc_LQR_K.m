function [K] = WBR_calc_LQR_K(Q, R, length_left, length_right, control_freq, STATS_FILE, A10_B4_FILE)
    % WBR_calc_LQR_K 计算指定腿长下的K矩阵（4 * 10），模型：WBR
    dt = 1 / control_freq;
    
    [sym_vars, sym_vals] = compute_state(length_left, length_right, STATS_FILE, A10_B4_FILE);
    
    A = A10_B4_FILE.sym_A;
    B = A10_B4_FILE.sym_B;
    
    A_Continuous = subs(A, sym_vars, sym_vals);
    B_Continuous = subs(B, sym_vars, sym_vals);
    
    [Ad,Bd] = c2d(double(A_Continuous), double(B_Continuous), dt);
    
    K = dlqr(Ad, Bd, Q, R, 0); % 状态和输入变量之间不存在耦合关系，不需要最后一项交叉加权矩阵

end

function [sym_vars, sym_vals] = compute_state(length_left, length_right, STATS_FILE, STATE_SPACE_FILE)
    %% 符号常量
    R_w = STATE_SPACE_FILE.R_w;
    R_l = STATE_SPACE_FILE.R_l;
    l_l = STATE_SPACE_FILE.l_l;
    l_r = STATE_SPACE_FILE.l_r;
    l_wl = STATE_SPACE_FILE.l_wl;
    l_wr = STATE_SPACE_FILE.l_wr;
    l_bl = STATE_SPACE_FILE.l_bl;
    l_br = STATE_SPACE_FILE.l_br;
    l_c = STATE_SPACE_FILE.l_c;
    m_w = STATE_SPACE_FILE.m_w;
    m_l = STATE_SPACE_FILE.m_l;
    m_b = STATE_SPACE_FILE.m_b;
    I_w = STATE_SPACE_FILE.I_w;
    I_ll = STATE_SPACE_FILE.I_ll;
    I_lr = STATE_SPACE_FILE.I_lr;
    I_b = STATE_SPACE_FILE.I_b;
    I_z = STATE_SPACE_FILE.I_z;
    g = STATE_SPACE_FILE.g;
    
    %% 实际值
    R_w_ac = STATS_FILE.radius_wheel;
    R_l_ac = STATS_FILE.distance_half_wheel;
    l_l_ac = length_left;
    l_r_ac = length_right;
    l_wl_ac = polyval(STATS_FILE.coeff_distance_wheel_leg, length_left);
    l_wr_ac = polyval(STATS_FILE.coeff_distance_wheel_leg, length_right);
    l_bl_ac =  polyval(STATS_FILE.coeff_distance_hip_leg, length_left);
    l_br_ac =  polyval(STATS_FILE.coeff_distance_hip_leg, length_right);
    l_c_ac = STATS_FILE.distance_body_hip;
    m_w_ac = STATS_FILE.mass_wheel_aver;
    m_l_ac = STATS_FILE.mass_leg_aver;
    m_b_ac = STATS_FILE.mass_body;
    I_w_ac = STATS_FILE.inertia_wheel_aver;
    I_ll_ac = polyval(STATS_FILE.coeff_inertia_leg, length_left);
    I_lr_ac = polyval(STATS_FILE.coeff_inertia_leg, length_right);
    I_b_ac = STATS_FILE.inertia_body;
    I_z_ac = STATS_FILE.inertia_z;
    g_ac = STATS_FILE.g_;
    
    %% compute
    sym_vars = [R_w, R_l, l_l, l_r, l_wl, l_wr, l_bl, l_br, l_c, m_w, m_l, m_b, I_w, I_ll, I_lr, I_b, I_z, g];
    sym_vals = [R_w_ac, R_l_ac, l_l_ac, l_r_ac, l_wl_ac, l_wr_ac, l_bl_ac, l_br_ac, l_c_ac, ...
    m_w_ac, m_l_ac, m_b_ac, I_w_ac, I_ll_ac, I_lr_ac, I_b_ac, I_z_ac, g_ac];
end

