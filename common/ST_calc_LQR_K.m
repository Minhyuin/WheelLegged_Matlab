function [K] = ST_calc_LQR_K(Q, R, length, control_freq, STATS_FILE, A4_B1_FILE)
    % ST_calc_LQR_K 计算指定腿长下的K矩阵（1 * 4），模型：ST
    dt = 1 / control_freq;

    [sym_vars, sym_vals] = compute_state(length, STATS_FILE, A4_B1_FILE);

    A = A4_B1_FILE.sym_A;
    B = A4_B1_FILE.sym_B;

    A_Continuous = subs(A, sym_vars, sym_vals);
    B_Continuous = subs(B, sym_vars, sym_vals);
    
    [Ad,Bd] = c2d(double(A_Continuous), double(B_Continuous), dt);
    
    K = dlqr(Ad, Bd, Q, R, 0); % 状态和输入变量之间不存在耦合关系，不需要最后一项交叉加权矩阵
end

function [sym_vars, sym_vals] = compute_state(length, STATS_FILE, STATE_SPACE_FILE)
    %% 符号常量
    m_w = STATE_SPACE_FILE.m_w;
    m_b = STATE_SPACE_FILE.m_b;
    l = STATE_SPACE_FILE.l;
    g = STATE_SPACE_FILE.g;

    %% 实际值
    m_w_ac = STATS_FILE.mass_wheel_aver;
    m_b_ac = STATS_FILE.mass_body / 2;
    l_ac = polyval(STATS_FILE.coeff_distance_wheel_leg, length);
    g_ac = STATS_FILE.g_;

    %% compute
    sym_vars = [m_w, m_b, l, g];
    sym_vals = [m_w_ac, m_b_ac, l_ac, g_ac];
end