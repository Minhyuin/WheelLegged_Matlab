function [K_sym, K_sim_param] = polyfitn_K(X_data, K_data, rank, is_evaluate)
%   根据给定的腿长和K矩阵序列，对K矩阵的每一个元素进行多项式拟合

%   输入参数:
%   X_data/K_data - 采样点，应是一个 （X, K）的张量, 其中K是因变量张量， X是自变量张量
%       格式要求: 
%       sample_idx为采样点数总数
%       K_data(sample_idx, ...) 
%       例如K是个n行m列的二维矩阵，即有n*m个因变量元素需要拟合
%       对于二维矩阵则此格式为
%       K_data = zeros(TOTAL_N_SAMPLES, N_u, N_x)
%       计算时K_data(sample_idx, :, :) = calc_K
%       三维、四维以此类推 
%       X_data(sample_idx, ...) 
%       例如需要进行左右腿的相关拟合，每次采样都有两个独立变量，可理解为1行2列的行向量
%       对于该行向量，格式为
%       X_data = zeros(TOTAL_N_SAMPLES, 2);
%       计算时X_data(sample_idx, :) = [length_left, length_right];
%
%   rank - 拟合阶数，阶数越高，拟合越精准，算力资源消耗也越高 
%   is_evaluate - true or false，决定是否需要输出统计特征
%
%   输出参数:
%   K_sym - 拟合出的表达式
%   K_sim_param - 多项式的所有系数

    fprintf("\n=== 正在运行" + mfilename + ".m ===\n");
    % 记录开始时间
    fit_start_time = tic;

    %% 提取采样点
    dims = size(K_data); % 拟合阶数
    total_samples = dims(1); % 采样点数
    N_u = dims(2);  % K矩阵行数
    N_x = dims(3);  % K矩阵列数
    poly_n = N_u * N_x;
    K_sym = sym(zeros(N_u,N_x)); 
    K_sim_param = zeros(N_u,N_x,rank+1,rank+1);

    %% 用于显示进度
    progress_step = max(1, round(poly_n / 100));
    last_percent = 0;

    %% 开始拟合
    fprintf('K拟合进度:')
    for i = 1 : N_u
        for j = 1 : N_x
            K = squeeze(K_data(:, i, j));
            X = X_data;

            polyval = polyfitn(X, K, rank);
            expon = polyval.ModelTerms;
            coef = polyval.Coefficients;
            nterms = length(coef);
            nvars = size(expon, 2);
            
            for k = 1 : nterms
                exp1 = expon(k, 1) + 1;
                exp2 = expon(k, 2) + 1;
                if exp1 <= rank+1 && exp2 <= rank+1
                    K_sim_param(i, j, exp1, exp2) = coef(k);
                else
                    error('多项式阶数超过了预分配的数组维度');
                end
            end

            K_sym(i, j) = simplify(polyn2sym(polyval));

            % 更新进度显示
            progress_step = progress_step + 1;
            current_percent = floor(progress_step / poly_n * 100);
            if current_percent > last_percent
                fprintf(repmat('#', 1, current_percent - last_percent));
                last_percent = current_percent;
            end
        end
    end

    if last_percent < 100
        fprintf(repmat('#', 1, 100 - last_percent));
    end

    %% 计算耗时
    fit_time = toc(fit_start_time);

    %% 输出拟合特征和统计信息
    fprintf('\n');
    fprintf('拟合完成！\n');
    fprintf('  - 拟合阶数: %d阶\n', rank);
    fprintf('  - 自变量维度: %d (左腿长, 右腿长)\n', size(X_data, 2));
    fprintf('  - 采样点数量: %d\n', total_samples);
    fprintf('  - K矩阵维度: %d×%d\n', N_u, N_x);
    fprintf('  - 待拟合元素总数: %d\n', poly_n);
    fprintf('  - 多项式项数: %d项\n', (rank+1)*(rank+2)/2); % 二元多项式的项数
    fprintf('  - 拟合总耗时: %.2f秒 (约%.2f分钟)\n', fit_time, fit_time/60);
    fprintf('  - 平均每个元素耗时: %.3f毫秒\n', fit_time/poly_n*1000);
    fprintf('\n');
    if (is_evaluate == true)
        evaluate_fit_results(X_data, K_data, K_sym, rank, true);
    end
end