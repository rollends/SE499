function dpp = polydiff( pp, order )
%POLYDIFF Differentiate Polynomial
%   Returns the coefficients of the polynomial equal to the order'th derivative of the
%   polynomial passed in.
    assert( order >= 0, 'Only non-zero positive integer derivative orders are permitted.');
    assert( isnumeric(order), 'Only non-zero positive integer derivative orders are permitted.');
    
    N = size(pp, 2) - 1;
    Nf = N - order;
    
    if order > N + 1
        % 0 identically.
        dpp = 0 * pp;
        return;
    end
    
    % Zero out the coefficients that are eliminated and shift.
    pp(:,(N+1-(order-1)):(N+1)) = 0;
    pp = circshift( pp, order, 2 );
    
    % Construct the power-multiplier that comes from differentiating powers of x.
    powers = Nf:-1:0;
    c = horzcat(    zeros(1, order), ...
                    factorial(order + powers) ./ factorial(powers) );
    
    dpp = pp .* repmat(c, size(pp, 1), 1);
end

