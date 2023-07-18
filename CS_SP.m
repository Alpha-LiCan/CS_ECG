
function hat_s = CS_SP(y, T, K)
    % 输入参数说明：
    % y：观测向量
    % T：字典矩阵
    % K：稀疏度

    % 初始化稀疏系数
    N = size(T, 2);
    hat_s = zeros(N, 1);

    % 初始化残差
    res = y;

    % 迭代追踪
    for k = 1:K
        % 计算投影系数
        proj = T' * res;
        
        % 找出具有最大投影系数的原子索引
        [~, idx] = max(abs(proj));

        % 更新稀疏系数
        hat_s(idx) = hat_s(idx) + proj(idx);
        
        % 更新残差
        res = y - T * hat_s;
    end
end
