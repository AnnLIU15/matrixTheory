function Y = Aopr(X)
    [m, n] = size(X);
    Y = zeros(2*m, 2*n);
    Y(1:m, 1:n) = X;
end