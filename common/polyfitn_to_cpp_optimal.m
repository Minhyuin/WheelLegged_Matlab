function polyfitn_to_cpp_optimal(polyn_matrix, name)
% 直接将多项式转换为C++表达式数组，使用预计算中间变量
% 目前支持rank = 3的操作，当rank增长时需要进行扩充修改
% 代码前需对腿长做预计算
% const float &x1 = _lL;
% const float &x2 = _lR;
% const float x1x1 = x1 * x1;
% const float x2x2 = x2 * x2;
% const float x1x1x1 = x1x1 * x1;
% const float x2x2x2 = x2x2 * x2;
% const float x1x2 = x1 * x2;
% const float x1x2x2 = x1x2 * x2;
% const float x1x1x2 = x1x1 * x2;
%
% 另外，在clang-tidy的环境下，还需要抑制警告
% 在生成的代码头一行加入
% // NOLINTBEGIN(readability-math-missing-parentheses, readability-identifier-naming)
% 在生成的代码尾一行加入
% // NOLINTEND(readability-math-missing-parentheses, readability-identifier-naming)

    fprintf("\n=== 正在运行" + mfilename + ".m ===\n");
    
    %% 开始转换 %%
    
    [rows, cols] = size(polyn_matrix);
    
    fprintf('float %s[%d * %d] = {\n', name, rows, cols);
    
    for i = 1:rows
        for j = 1:cols
            expr = char(polyn_matrix(i, j));
            
            % 将符号表达式转换为数值并格式化
            sym_expr = str2sym(expr);
            expr_decimal = vpa(sym_expr, 8); % 8位精度
            expr_str = char(expr_decimal);
            
            % 进行替换
            % 1. 替换X1, X2为x1, x2
            expr_str = strrep(expr_str, 'X1', 'x1');
            expr_str = strrep(expr_str, 'X2', 'x2');
            
            % 2. 替换幂运算为预计算变量
            expr_str = replace_powers(expr_str);
            
            % 3. 为浮点数添加f后缀
            expr_str = add_f_suffix(expr_str);
            
            % 4. 简化表达式
            expr_str = simplify_expr(expr_str);
            
            fprintf('    %s', expr_str);
            
            if ~(i == rows && j == cols)
                fprintf(',');
            end
            
            if j == cols && i < rows
                fprintf('\n');
            else
                fprintf(' ');
            end
        end
    end
    
    fprintf('\n};\n');
end

function expr = replace_powers(expr_str)
% 替换幂运算为预计算变量
    
    % 替换x1^2 -> x1x1
    expr_str = regexprep(expr_str, 'x1\s*\^\s*2', 'x1x1');
    expr_str = regexprep(expr_str, 'x1\s*\*\s*x1(?!\w)', 'x1x1');
    
    % 替换x2^2 -> x2x2
    expr_str = regexprep(expr_str, 'x2\s*\^\s*2', 'x2x2');
    expr_str = regexprep(expr_str, 'x2\s*\*\s*x2(?!\w)', 'x2x2');
    
    % 替换x1^3 -> x1x1x1
    expr_str = regexprep(expr_str, 'x1\s*\^\s*3', 'x1x1x1');
    expr_str = regexprep(expr_str, 'x1x1\s*\*\s*x1(?!\w)', 'x1x1x1');
    expr_str = regexprep(expr_str, 'x1\s*\*\s*x1x1(?!\w)', 'x1x1x1');
    
    % 替换x2^3 -> x2x2x2
    expr_str = regexprep(expr_str, 'x2\s*\^\s*3', 'x2x2x2');
    expr_str = regexprep(expr_str, 'x2x2\s*\*\s*x2(?!\w)', 'x2x2x2');
    expr_str = regexprep(expr_str, 'x2\s*\*\s*x2x2(?!\w)', 'x2x2x2');
    
    % 替换x1*x2 -> x1x2
    expr_str = regexprep(expr_str, 'x1\s*\*\s*x2', 'x1x2');
    expr_str = regexprep(expr_str, 'x2\s*\*\s*x1', 'x1x2');
    
    % 替换x1^2*x2 -> x1x1x2
    expr_str = regexprep(expr_str, 'x1\s*\^\s*2\s*\*\s*x2', 'x1x1x2');
    expr_str = regexprep(expr_str, 'x1x1\s*\*\s*x2', 'x1x1x2');
    
    % 替换x1*x2^2 -> x1x2x2
    expr_str = regexprep(expr_str, 'x1\s*\*\s*x2\s*\^\s*2', 'x1x2x2');
    expr_str = regexprep(expr_str, 'x1\s*\*\s*x2x2', 'x1x2x2');
    expr_str = regexprep(expr_str, 'x1x2\s*\*\s*x2', 'x1x2x2');
    
    expr = expr_str;
end

function expr = add_f_suffix(expr_str)
% 为浮点数添加f后缀
    
    % 匹配浮点数（包括负数和科学计数法）
    pattern = '(\-?\d+\.\d+(?:[eE][\+\-]?\d+)?)';
    
    % 查找所有浮点数
    [matches, tokens] = regexp(expr_str, pattern, 'match', 'tokens');
    
    for k = 1:length(matches)
        % 检查是否已经有f后缀
        if ~endsWith(matches{k}, 'f')
            % 添加f后缀
            new_num = [matches{k}, 'f'];
            expr_str = strrep(expr_str, matches{k}, new_num);
        end
    end
    
    expr = expr_str;
end

function expr = simplify_expr(expr_str)
% 简化表达式
    
    % 移除空格
    expr_str = strrep(expr_str, ' ', '');
    
    % 简化1*乘法
    expr_str = regexprep(expr_str, '1\*', '');
    expr_str = regexprep(expr_str, '\*1(?![0-9])', '');
    
    % 简化+和-运算符
    expr_str = strrep(expr_str, '+-', '-');
    expr_str = strrep(expr_str, '-+', '-');
    
    % 如果表达式以+开头，移除+
    if expr_str(1) == '+'
        expr_str = expr_str(2:end);
    end
    
    expr = expr_str;
end