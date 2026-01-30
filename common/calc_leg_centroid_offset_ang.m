function [ang] = calc_leg_centroid_offset_ang(length, distance_wheel_leg, distance_hip_leg)
    % calc_leg_centroid_offset_ang 计算腿部质心偏移角度（从腿摆杆方向指向质心，逆时针正方向）
    % 输入参数（全部为非零正整数的标量）：
    %   length: 腿摆杆长度
    %   distance_wheel_leg: 驱动轮质心到腿摆杆质心距离
    %   distance_hip_leg: 机体转轴到腿摆杆质心距离
    % 输出参数：
    %   ang: 所求的质心偏移角度，单位：rad
    %
    % 主要通过几何关系求解，已知length，distance_wheel_leg，distance_hip_leg围成三角形
    % 使用余弦定理求即可

    % 余弦定理：c^2 = a^2 + b^2 - 2*a*b*cos(C)
    % 其中，distance_wheel_leg为c，length为a，distance_hip_leg为b
    % 求解角度C（即质心偏移角度）
    cos_ang = (length^2 + distance_hip_leg^2 - distance_wheel_leg^2) / (2 * length * distance_hip_leg);
    cos_ang = max(-1, min(1, cos_ang));
    ang = acos(cos_ang);
    % 因一般情况下，串联腿质心靠后，故可以认为质心偏移角度为正

end