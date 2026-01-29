function [K] = WL_calc_LQR_K(Q, R, length, control_freq, STATS_FILE, A6_B2_FILE)
    %WL_CALC_LQR_K 计算指定腿长下的K矩阵（2 * 6），模型：WL
    %   此处显示详细说明
    dt = 1 / control_freq;
    
    [sym_vars, sym_vals] = compute_state(length, STATS_FILE, A6_B2_FILE);
    
    A = A6_B2_FILE.sym_A;
    B = A6_B2_FILE.sym_B;

    A_Continuous = subs(A, sym_vars, sym_vals);
    B_Continuous = subs(B, sym_vars, sym_vals);
    
    [Ad,Bd] = c2d(double(A_Continuous), double(B_Continuous), dt);
    
    if 0
        if any(isnan(Ad(:))) || any(isinf(Ad(:)))
            disp(Ad);
            error('Ad异常');
        end
        if any(isnan(Bd(:))) || any(isinf(Bd(:)))
            disp(Bd);
            error('Bd异常');
        end
    end
    
    K = dlqr(Ad, Bd, Q, R, 0); % 状态和输入变量之间不存在耦合关系，不需要最后一项交叉加权矩阵

end

function [sym_vars, sym_vals] = compute_state(length, STATS_FILE, STATE_SPACE_FILE)
    %% 符号常量
    R = STATE_SPACE_FILE.R;
    L = STATE_SPACE_FILE.L;
    L_M = STATE_SPACE_FILE.L_M;
    l = STATE_SPACE_FILE.l;
    m_w = STATE_SPACE_FILE.m_w;
    m_p = STATE_SPACE_FILE.m_p;
    M = STATE_SPACE_FILE.M;
    I_w = STATE_SPACE_FILE.I_w;
    I_p = STATE_SPACE_FILE.I_p;
    I_M = STATE_SPACE_FILE.I_M;
    g = STATE_SPACE_FILE.g;

    %% 实际值
    R_ac = STATS_FILE.radius_wheel;
    L_ac = polyval(STATS_FILE.coeff_distance_wheel_leg, length);
    L_M_ac = polyval(STATS_FILE.coeff_distance_hip_leg, length);
    l_ac  = STATS_FILE.distance_body_hip;
    m_w_ac = STATS_FILE.mass_wheel_aver;
    m_p_ac = STATS_FILE.mass_leg_aver;
    M_ac = STATS_FILE.mass_body / 2;
    I_w_ac =  STATS_FILE.inertia_wheel_aver;
    I_p_ac = polyval(STATS_FILE.coeff_inertia_leg, length);
    I_M_ac = STATS_FILE.inertia_body / 2;
    g_ac = STATS_FILE.g_;

    %% compute
    sym_vars = [R, L, L_M, l, m_w, m_p, M, I_w, I_p, I_M, g];
    sym_vals = [R_ac, L_ac, L_M_ac, l_ac, m_w_ac, m_p_ac, M_ac, I_w_ac, I_p_ac, I_M_ac, g_ac];

end
