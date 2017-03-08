function v = vpolyval(p, x)
%VPOLYVAL Vector Polynomial Evaluator
%   Because the standard polyval doesn't support vector polynomials :(
    x = x(:);
    v = zeros(size(p, 1), size(x, 1));
    if isscalar(x)
        v(:) = p * (x .^ (size(p, 2)-1:-1:0) )';
    else
        for vi = 1:size(v, 1)
            v(vi, :) = polyval(p(vi, :), x);
        end
    end
end