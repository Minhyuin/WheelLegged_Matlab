function polyfitn_to_cpp(polyn_matrix,name)
% 直接将多项式转换为C++表达式数组，未使用预计算中间变量
% 但可以支持任意的rank

    fprintf("\n=== 正在运行" + mfilename + ".m ===\n");
    
    %% 开始转换 %%
    fprintf('float %s[%d * %d] = {\n(',name,size(polyn_matrix, 1),size(polyn_matrix, 2));
    for i = 1:size(polyn_matrix, 1)
        %fprintf('    {');
        for j = 1:size(polyn_matrix, 2)
            expr = char(polyn_matrix(i, j));
            
            % 使用 str2sym 将字符串转换为符号表达式
            sym_expr = str2sym(expr);
            
            % 使用 vpa 函数将符号表达式转换为指定精度的小数形式
            expr_decimal = vpa(sym_expr, 10); % 这里的 6 是精度，您可以根据需要调整
            
            % 将符号表达式转换为字符串形式
            expr_str = char(expr_decimal);
            
            % 使用 regexprep 进行替换
            new_expr = regexprep(expr_str, '(\w+)\^(\d+)', 'powf($1, $2)');
            
            % 将小数格式化为字符串，并手动添加 'f'
            new_expr_with_f = regexprep(new_expr, '(\d+\.\d+)', '$1f');
    
            add1 = strrep(new_expr_with_f, ' + ', ') + (');
    
            add2 = strrep(add1, ' - ', ') - (');
    
            replaceX = strrep(add2, 'X', 'x');
            
            fprintf('%s', replaceX);
            if j < size(polyn_matrix, 2)
                fprintf('), (');
            end
        end
        %fprintf('}');
        if i < size(polyn_matrix, 1)
            fprintf('),');
        end
        if i == size(polyn_matrix, 1) && j == size(polyn_matrix, 2)
            fprintf('\n');
        else
            fprintf('\n(');
        end
    end
    fprintf(')};\n');

end