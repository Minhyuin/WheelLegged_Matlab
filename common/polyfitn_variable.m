function [Y_sym, Y_sim_param] = polyfitn_variable(X_data, Y_data, rank)
%   根据给定的X和Y序列，对Y进行多项式拟合
    fprintf("\n=== 正在运行" + mfilename + ".m ===\n");
    % 记录开始时间
    fit_start_time = tic;

    %% 提取采样点
    total_samples = numel(Y_data); % 采样点数
    Y_sim_param = zeros(rank+1, rank+1);

    %% 开始拟合
    polyval = polyfitn(X_data, Y_data, rank);
    expon = polyval.ModelTerms;
    coef = polyval.Coefficients;
    nterms = length(coef);
    nvars = size(expon, 2);
    expon_full = zeros(nterms, 2); % 最多处理2个变量
    expon_full(:, 1:min(2, nvars)) = expon(:, 1:min(2, nvars));
    for k = 1:nterms
        exp1 = expon_full(k, 1) + 1;
        exp2 = expon_full(k, 2) + 1;

        if exp1 <= rank + 1 && exp2 <= rank + 1
            Y_sim_param(exp1, exp2) = coef(k);
        else
            error('多项式阶数超过了预分配的数组维度');
        end
    end
    Y_sym = simplify(polyn2sym(polyval));

    %% 显示计算时间
    fit_time = toc(fit_start_time);

    %% 拟合特征和统计信息
    fprintf('多项式拟合完成，共计 %d 个采样点，拟合阶数 %d，计算时间 %.2f 秒。\n', ...
        total_samples, rank, fit_time);
end