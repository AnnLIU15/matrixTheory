function [ X_opt, iter ] = admmapAXB( A, B, X, M, known, para)    

max_iter = para.admmap_iter;
tol = para.admmap_tol;
beta = para.admmap_beta;
kappa = para.admmap_kappa;
rho0 = para.admmap_rho;
betamax = para.admmap_betamax;

AB = A' * B;
Y = Aopr(X);
W = X;
M_fro = norm(M, 'fro');
obj_val = zeros(max_iter, 1);

[m, n] = size(X);
C = zeros(2*m, 2*n);
C(m+1:2*m, n+1:2*n) = X;

for k = 1 : max_iter
    % X update
    last_X = X;
    temp = W - 1/beta * Aadj(Y);
    [U, sigma, V] = svd(temp);
    X = U * (sign(sigma) * max(abs(sigma) - 1/beta, 0)) * V';

    % W update
    last_W = W;
    W = 1/(2*beta) * (beta * (M - X) - (AB + Y(1:m, 1:n) + Y(m+1:2*m, n+1:2*n))) .* known...
        + X + 1/beta * (AB + Aadj(Y));

    % Y update
    Y = Y + beta * (Aopr(X)+Bopr(W, known)-C);
    
    % rho update
    temp = beta * max(norm(X - last_X, 'fro'), norm(W - last_W, 'fro')) / norm(C, 'fro');
    if temp < kappa
        rho = rho0;
    else
        rho = 1;
    end
    
    % beta update
    beta = min(betamax, rho*beta);

    delta = norm(X - last_X, 'fro') / M_fro;

    if delta < tol
        fprintf('    converged at\n');
        fprintf('    iter %d, ||X_k+1-X_k||_F/||M||_F=%.4f\n', k, delta);
        break ;
    end
    
    obj_val(k) = -trace(A*W*B') + trace(Y*(Aopr(X)+Bopr(W, known)-C)')...
        + nuclear_norm(X) + beta / 2 * norm(Aopr(X)+Bopr(W, known)-C,'fro')^2;

    if k > 1 && abs(obj_val(k) - obj_val(k-1)) < tol
        fprintf('    converged at\n');
        fprintf('    iter %d, obj value=%.4f\n', k, obj_val(k));        
        break;
    end
    if k == max_iter
       fprintf('beta = %.4f\n', beta);
end

X_opt = X;
iter = k;

end