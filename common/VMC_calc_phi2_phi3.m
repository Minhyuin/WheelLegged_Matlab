function [phi2, phi3] = VMC_calc_phi2_phi3(phi1, phi4, L1, L2, L3, L4, L5)
%CALC_PHI2_PHI3 输入phi1, phi4以及各腿长，输出phi2, phi3
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
end

