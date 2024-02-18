function [outputimg] = truncateslice(inputimg,width)

%channel =2*width, 5, 3
% 获取原始数组的 中间位置切面 的起始和结束索引
middle_index = ceil(size(inputimg, 3) / 2);
start_index = middle_index - width+1;
end_index = middle_index + width;

% 从原始数组中选择连续的二维slices，并生成新的三维矩阵
outputimg = inputimg(:, :, start_index:end_index);


end

