function [tp, fd] = VMC_torque_inverse_kinematics(tau1, tau4, phi1, phi4, L1, L2, L3, L4, L5)
    %VMC_torque_inverse_kinematics 输入摆杆角度和长度以及用切向力和法向力表示的力矢量，输出并联关节力矩
    %   注意，摆杆长度为正数，摆杆角度零点位置为摆杆竖直，详细见
    %   https://jjyczursn9.feishu.cn/wiki/IW2BwNm7GiU5gZkxX96c0rd5nxe?from=from_copylink
    
    %% temporary I
    sin1 = sin(phi1);
    cos1 = cos(phi1);
    sin4 = sin(phi4);
    cos4 = cos(phi4);
    
    %% calc phi2 and phi3
    xA = -L5 / 2;
    yA = 0;
    xE = L5 / 2;
    yE = 0;
    xB = xA + L1 * cos1;
    yB = yA + L1 * sin1;
    xD = xE + L4 * cos4;
    yD = yE + L4 * sin4;
    
    xBD = xD - xB;
    yBD = yD - yB;
    lBD = sqrt((xBD * xBD) + (yBD * yBD));
    
    cosBeta = ((L2 * L2 + lBD * lBD) - (L3 * L3)) / (2 * L2 * lBD);
    beta = acos(cosBeta);
    alpha = atan2(yBD, xBD);
    
    phi2 = beta + alpha;
    phi3 = pi - beta + alpha;
    
    %% temporary II
    sin2 = sin(phi2);
    cos2 = cos(phi2);
    sin3 = sin(phi3);
    cos3 = cos(phi3);
    
    %% calc phi0 and l0
    xC = L2 * cos2 + L1 * cos1 - L5 / 2;
    yC = L2 * sin2 + L1 * sin1;
    phiC = atan2(yC, xC);
    L0 = sqrt((xC * xC) + (yC * yC));

    %% assemble J Matrix
    sin21 = sin(phi1 - phi2);
    sin43 = sin(phi3 - phi4);
    sin32 = sin(phi2 - phi3);
    J = zeros(2, 2);
    J(1, 1) = +L1 * sin3 * sin21 / sin32;
    J(1, 2) = +L4 * sin2 * sin43 / sin32;
    J(2, 1) = -L1 * cos3 * sin21 / sin32;
    J(2, 2) = -L4 * cos2 * sin43 / sin32; 
    J_T = J';
    J_T_inv = inv(J_T);

    %% assemble R Matrix
    R = [-sin(phiC), cos(phiC); cos(phiC), sin(phiC);];
    R_inv = inv(R);

    %% assemble M Matrix
    M = zeros(2, 2);
    M(1, 1) = 1 / L0;
    M(1, 2) = 0;
    M(2, 1) = 0;
    M(2, 2) = 1;
    M_inv = inv(M);

    output = zeros(2, 1); %#ok<*PREALL>
    input = [tau1; tau4];
    output = M_inv * R_inv * J_T_inv * input; %#ok<*MINV>

    tp = output(1, 1);
    fd = output(2, 1);
end

