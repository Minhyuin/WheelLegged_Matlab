function cpp_code = matrix_to_cpp(matrix, var_name)
% 将MATLAB矩阵转换为C++风格的一维数组字符串（按行优先顺序展开）
% 输入:
%   matrix: MATLAB矩阵
%   var_name: 变量名
% 输出:
%   cpp_code: C/C++风格的一维组字符串

    [rows, cols] = size(matrix);
    total_elements = rows * cols;

    cpp_code = sprintf('float %s[%d * %d] = {\n    ', var_name, rows, cols);
    
    % 先遍历行，再遍历列
    for i = 1:rows
        for j = 1:cols
            % 格式化数字
            value = matrix(i, j);
            
            % 保留四位小数
            num_str = sprintf('%.4f', value);
            
            cpp_code = [cpp_code, num_str];
            
            % 添加逗号（除了最后一个元素）
            if ~(i == rows && j == cols)
                cpp_code = [cpp_code, ', '];
            end
            
            % 元素换行
            if mod((i-1)*cols + j, cols) == 0 && ~(i == rows && j == cols)
                cpp_code = [cpp_code, '\n    '];
            end
        end
    end
    
    cpp_code = [cpp_code, '\n};\n'];
    
    fprintf(cpp_code);
end