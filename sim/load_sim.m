function result = load_sim(model_type, show_plot)
% 根据模型选择加载相应的仿真参数
% - model_type 模型选择，0: WBR(全身建模) | 1: WL(单腿建模) | 2: ST(板凳建模)
% - show_plot 是否绘制VMC解算示意图
    result = struct();
    if (model_type == 0)
        result = load('WBR_sim_param.mat');
        result.L1 = result.WBR_SIM.robot_stats.length_thigh;
        result.L2 = result.WBR_SIM.robot_stats.length_calf;
        result.L3 = result.WBR_SIM.robot_stats.length_calf;
        result.L4 = result.WBR_SIM.robot_stats.length_thigh;
        result.L5 = 0.01;
        result.init_phi1 = pi - 0.4;
        result.init_phi4 = 0.4;
        [result.init_L0, result.init_phi0, result.init_phi2, result.init_phi3] = VMC_plot(result.init_phi1, result.init_phi4, result.L1, result.L2, result.L3, result.L4, result.L5, show_plot);
    end
    result.init_vx = 0;
    result.init_z = 1;
    init_leg_length = 0.2;
    init_leg_phi0 = deg2rad(0);
end

function a = deg2rad(b)
    a = b / 180 * pi;
end