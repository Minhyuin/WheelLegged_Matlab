function get_matrix_param_in_C_code(K_LQR,name)

fprintf('float %s[%d*%d] = {\n',name,size(K_LQR, 1),size(K_LQR, 2));
for i = 1:size(K_LQR, 1)
    fprintf('    {');
    for j = 1:size(K_LQR, 2)
        expr = char(K_LQR(i, j));
        
        % 使用 str2sym 将字符串转换为符号表达式
        sym_expr = str2sym(expr);
        
        % 使用 vpa 函数将符号表达式转换为指定精度的小数形式
        expr_decimal = vpa(sym_expr, 20); % 这里的 6 是精度，您可以根据需要调整
        
        % 将符号表达式转换为字符串形式
        expr_str = char(expr_decimal);
        
        % 使用 regexprep 进行替换
        new_expr = regexprep(expr_str, '(\w+)\^(\d+)', 'pow$2f($1)');
        
        % 将小数格式化为字符串，并手动添加 'f'
        new_expr_with_f = regexprep(new_expr, '(\d+\.\d+)', '$1f');
        
        fprintf('%s', new_expr_with_f);
        if j < size(K_LQR, 2)
            fprintf(', ');
        end
    end
    fprintf('}');
    if i < size(K_LQR, 1)
        fprintf(',');
    end
    fprintf('\n');
end
fprintf('};\n');

end