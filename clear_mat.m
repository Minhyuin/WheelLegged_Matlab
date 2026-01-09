% QUICK_CLEAR_MAT 快速清除根目录下的.mat文件并打印文件名
clear; clc;

fprintf('快速清除.mat文件脚本\n');
fprintf('当前目录: %s\n\n', pwd);

% 获取所有.mat文件
mat_files = dir('*.mat');

if isempty(mat_files)
    fprintf('未找到.mat文件。\n');
    return;
end

fprintf('找到 %d 个.mat文件:\n', length(mat_files));
fprintf('----------------------------------------\n');

for i = 1:length(mat_files)
    try
        filename = mat_files(i).name;
        delete(filename);
        fprintf('✓ 已删除: %s\n', filename);
    catch ME
        fprintf('✗ 删除失败: %s\n', filename);
    end
end

fprintf('\n操作完成。\n');