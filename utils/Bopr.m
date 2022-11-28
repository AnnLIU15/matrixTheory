function Y = Bopr(X, k)
    [m, n] = size(X);
    Y = zeros(2*m, 2*n);
    Y(1:m, 1:n) = -X;
    Y(m+1:2*m, n+1:2*n) = X .* k;
end