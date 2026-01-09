function generate_A10_B4()
    % 10输入4输出状态空间方程AB矩阵求解
    % 脚本会解出A、B矩阵的表达式sym_A和sym_B
    % 写法参考港科的开源，核心还是解上交的方程

    fprintf("\n=== 正在运行" + mfilename + ".m ===\n");
    
    %% 定义参数
    % constant ↓
    syms R_w                % 驱动轮半径
    syms R_l                % 驱动轮轮距/2
    syms l_l l_r            % 左右腿长
    syms l_wl l_wr          % 驱动轮质心到左右腿部质心距离
    syms l_bl l_br          % 机体转轴到左腿摆杆质心距离
    syms l_c                % 机体质心到腿部关节中心点距离
    syms m_w m_l m_b        % 驱动轮质量 腿部质量 机体质量
    syms I_w                % 驱动轮转动惯量           (自然坐标系法向)
    syms I_ll I_lr          % 驱动轮左右腿部转动惯量    (自然坐标系法向，实际上会变化)
    syms I_b                % 机体转动惯量             (自然坐标系法向)
    syms I_z                % 机器人z轴转动惯量        (简化为常量)]
    syms g                  % 重力加速度
    % status variable ↓
    syms s                  % 质心位移
    syms d_s                % 质心速度
    syms phi                % 偏航角yaw
    syms d_phi              % 偏航角速度
    syms theta_ll           % 左腿倾斜角             (自然坐标系法向)
    syms d_theta_ll         % 左腿倾斜角速度         (自然坐标系法向)
    syms theta_lr           % 右腿倾斜角             (自然坐标系法向)
    syms d_theta_lr         % 右腿倾斜角速度         (自然坐标系法向)
    syms theta_b            % 机体倾斜角             (自然坐标系法向)
    syms d_theta_b          % 机体倾斜角速度         (自然坐标系法向)
    % independence variable ↓
    syms theta_wl   theta_wr   % 左右驱动轮转角
    syms d_theta_wl  d_theta_wr   
    syms dd_theta_wl dd_theta_wr 
    syms dd_theta_ll dd_theta_lr dd_theta_b
    % control variable ↓
    syms T_wl T_wr T_bl T_br
    
    %% 求解5个广义坐标的二阶导数
    %（https://bbs.robomaster.com/forum.php?mod=viewthread&tid=22756）
    eqn3p11 = (I_w*l_l/R_w+m_w*R_w*l_l+m_l*R_w*l_bl)*dd_theta_wl+(m_l*l_wl*l_bl-I_ll)*dd_theta_ll+(m_l*l_wl+m_b*l_l/2)*g*theta_ll+T_bl-T_wl*(1+l_l/R_w)==0;
    eqn3p12 = (I_w*l_r/R_w+m_w*R_w*l_r+m_l*R_w*l_br)*dd_theta_wr+(m_l*l_wr*l_br-I_lr)*dd_theta_lr+(m_l*l_wr+m_b*l_r/2)*g*theta_lr+T_br-T_wr*(1+l_r/R_w)==0;
    eqn3p13 = -(m_w*R_w*R_w+I_w+m_l*R_w*R_w+m_b*R_w*R_w/2)*dd_theta_wl-(m_w*R_w*R_w+I_w+m_l*R_w*R_w+m_b*R_w*R_w/2)*dd_theta_wr-(m_l*R_w*l_wl+m_b*R_w*l_l/2)*dd_theta_ll-(m_l*R_w*l_wr+m_b*R_w*l_r/2)*dd_theta_lr+T_wl+T_wr==0;
    eqn3p14 = (m_w*R_w*l_c+I_w*l_c/R_w+m_l*R_w*l_c)*dd_theta_wl+(m_w*R_w*l_c+I_w*l_c/R_w+m_l*R_w*l_c)*dd_theta_wr+m_l*l_wl*l_c*dd_theta_ll+m_l*l_wr*l_c*dd_theta_lr-I_b*dd_theta_b+m_b*g*l_c*theta_b-(T_wl+T_wr)*l_c/R_w-(T_bl+T_br)==0;
    eqn3p15 = ((I_z*R_w)/(2*R_l)+I_w*R_l/R_w)*dd_theta_wl-((I_z*R_w)/(2*R_l)+I_w*R_l/R_w)*dd_theta_wr+(I_z*l_l)/(2*R_l)*dd_theta_ll-(I_z*l_r)/(2*R_l)*dd_theta_lr-T_wl*R_l/R_w+T_wr*R_l/R_w==0;
    [dd_theta_wl,dd_theta_wr,dd_theta_ll,dd_theta_lr,dd_theta_b] = solve([eqn3p11,eqn3p12,eqn3p13,eqn3p14,eqn3p15],[dd_theta_wl,dd_theta_wr,dd_theta_ll,dd_theta_lr,dd_theta_b]);
    
    %% 通过计算雅可比矩阵的方法解所需偏导数
    J_A = jacobian([dd_theta_wl,dd_theta_wr,dd_theta_ll,dd_theta_lr,dd_theta_b],[theta_ll,theta_lr,theta_b]);
    J_B = jacobian([dd_theta_wl,dd_theta_wr,dd_theta_ll,dd_theta_lr,dd_theta_b],[T_wl,T_wr,T_bl,T_br]);
    
    %% 状态空间模型
    sym_A = sym('A',[10 10]);
    sym_B = sym('B',[10 4]);
    % 港科+上交的快速AB求解器
    % 填入A数据：a25,a27,a29,a45,a47,a49,a65,a67,a69,a85,a87,a89,a105,a107,a109
    for p = 5:2:9
        A_index = (p - 3)/2;
        sym_A(2,p) = R_w*(J_A(1,A_index) + J_A(2,A_index))/2;
        sym_A(4,p) = (R_w*(- J_A(1,A_index) + J_A(2,A_index)))/(2*R_l) - (l_l*J_A(3,A_index))/(2*R_l) + (l_r*J_A(4,A_index))/(2*R_l);
        for q = 6:2:10
            sym_A(q,p) = J_A(q/2,A_index);
        end
    end
    % A的以下数值为1：a12,a34,a56,a78,a910，其余数值为0
    for r = 1:10
        if rem(r,2) == 0
            sym_A(r,1) = 0; sym_A(r,2) = 0; sym_A(r,3) = 0; sym_A(r,4) = 0; sym_A(r,6) = 0; sym_A(r,8) = 0; sym_A(r,10) = 0;
        else
            sym_A(r,:) = zeros(1,10);
            sym_A(r,r+1) = 1;
        end
    end
    % 填入B数据：b21,b22,b23,b24,b41,b42,b43,b44,b61,b62,b63,b64,b81,b82,b83,b84,b101,b102,b103,b104,
    for h = 1:4
        sym_B(2,h) = R_w*(J_B(1,h) + J_B(2,h))/2;
        sym_B(4,h) = (R_w*(- J_B(1,h) + J_B(2,h)))/(2*R_l) - (l_l*J_B(3,h))/(2*R_l) + (l_r*J_B(4,h))/(2*R_l);
        for f = 6:2:10
            sym_B(f,h) = J_B(f/2,h);
        end
    end
    
    % B的其余数值为0
    for e = 1:2:9
        sym_B(e,:) = zeros(1,4);
    end
    
    %% 保存参数
    save A10_B4.mat
end
