function M_approx= OPTSPACE(M, r)
% OPTSPACE 算法
% 输入：
%     M: 待恢复的矩阵the observed matrix
%     r: 希望恢复出来的矩阵的秩
% 输出：
%     M_approx: 恢复矩阵
%% trimming
[m, n] = size(M);
E = logical(M);
M_E = M;

num = sum(sum(M_E)); % 矩阵中的总degree

avg_row = 2 * num/ m;
avg_col = 2 * num/ n;

M_row = sum(M_E, 2);
M_col = sum(M_E, 1);

index_row = M_row > avg_row;
index_col = M_col > avg_col;

M_trim = M;
M_trim(index_row,:) = 0;
M_trim(:, index_col) = 0;
%%
[X, S, Y] = svds(M_trim, r);

X = sqrt(m) * X;
Y = sqrt(n) * Y;

%S_0 = S(1:r, 1:r)/eps;  %运算中不需要
%%
MaxIter = 1000;
tol = 1e-4;
tau = 1e-3; % 步长初值

for k = 1: MaxIter
    S = optimize_F_S(X, Y, M_E, E);

    grad_X = (E.* (X * S * Y' - M_E)) * Y * S';
    grad_Y = (E.* (X * S * Y' - M_E))' * X * S;
    
    t = tau;
    
    for i = 1: 50
        if(cost_F(X - t*grad_X, Y - t*grad_Y,S,M_E, E) - cost_F(X, Y,S,M_E, E) < ...
                t * 0.5 * (norm(grad_X, 'fro')^2 + norm(grad_Y, 'fro')^2))
            break;
        end
        t = t/2;
    end
    
    X = X - t*grad_X;
    Y = Y - t*grad_Y;
    
    eps = norm(E.*(M_E - X*S*Y'), 'fro') / norm(M_E, 'fro');
    % fprintf('%dth iteration: %e \n', k, eps);
    
    if (eps < tol)
        % fprintf('convergent');
        break;
    end
end

M_approx = X*S*Y';
end
%%
function S = optimize_F_S(X, Y, M_E, E)
% X, Y: the input value of variable in the function F
% M_E： the observed matrix
% E: the mask matrix
% S：the matrix S that minimize the function F given value X, Y
% 想着化成A*S_ij = b的形式求解

[~, k] = size(X);
b = X'* M_E * Y;
b = b(:); %拉成列向量

A = zeros(k*k, k*k);
for i = 1: k
    for j = 1:k
        index = (j-1) * k + i;
        coef = X' * (X(:,i) * Y(:, j)' .* E) * Y;
        A(:, index) = coef(:);
    end
end
S = A\b;

S = reshape(S, k, k);
end
%%
function F = cost_F(X, Y, S, M_E, E)
%计算损失函数
F = sum( sum( ( (X*S*Y' - M_E).*E ).^2 ) )/2 ;
end