function [L0, phi0, phi2, phi3] = VMC_position_inverse_kinematics(phi1, phi4, L1, L2, L3, L4, L5)
    %VMC_position_inverse_kinematics 输入摆杆角度和长度，输出各关节角度
    %   注意，摆杆长度为正数，摆杆角度零点位置为摆杆竖直，详细见
    %   https://jjyczursn9.feishu.cn/wiki/IW2BwNm7GiU5gZkxX96c0rd5nxe?from=from_copylink

    fprintf("\n=== 正在运行" + mfilename + ".m ===\n");
    
    % temporary I
    sin1 = sin(phi1);
    cos1 = cos(phi1);
    sin4 = sin(phi4);
    cos4 = cos(phi4);
    
    % calc phi2 and phi3
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
    
    % temporary II
    sin2 = sin(phi2);
    cos2 = cos(phi2);
    sin3 = sin(phi3);
    cos3 = cos(phi3);
    
    % calc phi0 and l0
    xC = L2 * cos2 + L1 * cos1 - L5 / 2;
    yC = L2 * sin2 + L1 * sin1;
    phiC = atan2(yC, xC);
    phi0 = phiC - pi / 2;
    L0 = sqrt((xC * xC) + (yC * yC));

end

