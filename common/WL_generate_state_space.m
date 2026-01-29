function WL_generate_state_space()
    % 6输入2输出状态空间方程AB矩阵求解
    % 脚本会解出A、B矩阵的表达式sym_A和sym_B
    % 根据王工文章公式推导 得到状态空间的符号表达式，部分变量的命名有出入
    % https://zhuanlan.zhihu.com/p/563048952
    %% 一些与whx文章不同定义的解释
    % 我的状态变量顺序
    % 腿相对于竖直方向的夹角 phi
    % 腿相对于竖直方向的角速度 dphi
    % 轮毂位移 s
    % 轮毂速度 ds
    % 底盘仰俯角 q
    % 底盘仰俯角速度 dq

    % 王工的状态变量顺序
    % 腿相对于竖直方向的夹角 theta
    % 腿相对于竖直方向的角速度 dtheta
    % 轮毂位移 s
    % 轮毂速度 ds
    % 底盘仰俯角 phi
    % 底盘仰俯角速度 dphi

    % 这个描述动力学运动学的文件用王工的规则

    fprintf("\n=== 正在运行" + mfilename + ".m ===\n");

    %% 符号变量定义
    syms t
    syms theta(t) x(t) phi(t)
    syms T T_p N P N_M P_M N_f g; %variable
    syms R          % 驱动轮半径
    syms L          % 摆杆重心到驱动轮轴距离
    syms L_M    % 摆杆重心到机体转轴距离
    syms l          % 机体重心到其转轴距离
    syms m_w    % 驱动轮质量
    syms m_p    % 摆杆质量
    syms M      % 机体质量
    syms I_w    % 驱动轮转动惯量
    syms I_p    % 摆杆绕质心转动惯量
    syms I_M    % 机体绕质心转动惯量
    syms d_d_x d_x d_d_theta d_theta d_d_phi d_phi %status variable
    d_syms = [d_d_x,d_x,d_d_theta,d_theta,d_d_phi,d_phi];
    d_syms_raw = [diff(x,t,2),diff(x,t),diff(theta,t,2),diff(theta,t),diff(phi,t,2),diff(phi,t)];
    
    %% 经典力学分析
    % 对驱动轮,化简消得:
    eqn1p3 = diff(x,t,2) == (T-N*R)/(I_w/R+m_w*R);
    eqn1p3 = subs(eqn1p3,d_syms_raw,d_syms);
    % 对摆杆,有:
    eqn1p4 = N - N_M == m_p*diff(x+L*sin(theta),t,2);
    eqn1p4 = subs(eqn1p4,d_syms_raw,d_syms);
    
    eqn1p5 = P - P_M -m_p*g ==  m_p*diff(L*cos(theta),t,2);
    eqn1p5 = subs(eqn1p5,d_syms_raw,d_syms);
    
    eqn1p6 = I_p*diff(theta,t,2) == (P*L+P_M*L_M)*sin(theta)-(N*L+N_M*L_M)*cos(theta)-T+T_p;
    eqn1p6 = subs(eqn1p6,d_syms_raw,d_syms);
    % 对机体,有:
    eqn1p7 = N_M == M*diff(x+(L+L_M)*sin(theta)-l*sin(phi),t,2);
    eqn1p7 = subs(eqn1p7,d_syms_raw,d_syms);                    
    
    eqn1p8 = P_M - M*g == M*diff((L+L_M)*cos(theta)+l*cos(phi),t,2);
    eqn1p8 = subs(eqn1p8,d_syms_raw,d_syms);
    eqn1p9 = I_M*diff(phi,t,2) == T_p + N_M*l*cos(phi)+P_M*l*sin(phi);
    eqn1p9 = subs(eqn1p9,d_syms_raw,d_syms);
    
    [N,P,N_M,P_M,d_d_theta,d_d_x,d_d_phi] = solve(eqn1p7,eqn1p8,eqn1p4,eqn1p5,eqn1p3,eqn1p6,eqn1p9, N,P,N_M,P_M,d_d_theta,d_d_x,d_d_phi);

    %% 状态空间模型
    status_vec = [theta;d_theta;x;d_x;phi;d_phi];
    ctrl_vec = [T;T_p];
    d_status_vec = [d_theta;d_d_theta;d_x;d_d_x;d_phi;d_d_phi];
    
    % 定义平衡点状态
    % exp_status_vec = [0;0;0;diff(x,t);0;0;]; 
    exp_status_vec = [0;0;x;0;0;0;];
    
    % 在平衡点处对其线性化
    sym_A_raw=jacobian(d_status_vec,status_vec);
    sym_B_raw=jacobian(d_status_vec,ctrl_vec);
    sym_A = subs(sym_A_raw,status_vec,exp_status_vec);
    sym_B = subs(sym_B_raw,status_vec,exp_status_vec);

    %% 保存参数
    save A6_B2.mat
end

