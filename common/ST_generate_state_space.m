function WBR_generate_state_space()
    % 4输入1输出状态空间方程AB矩阵求解
    % 脚本会解出A、B矩阵的表达式sym_A和sym_B

    fprintf("\n=== 正在运行" + mfilename + ".m ===\n");

    %% 定义参数
    syms m_w m_b
    syms l
    syms g
    syms sym_A
    syms sym_B

    %% 求解A、B矩阵表达式
    sym_A = [0                    0                     1 0;
             0                    0                     0 1;
             0             m_b*g/m_w                    0 0;
             0          (m_b+m_w)*g/m_w/l               0 0];
    sym_B = [          0;
                       0;
                    1/m_w;
                    1/m_w/l];
    save A4_B1.mat
end
