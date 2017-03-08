function O = gramschmidt( E )
%GRAMSCHMIDT Performs Gram-Schmidt Orthoganilization 
%   Generates an orthonormal basis out of the columns of E
    
    O = E;
    O(:, 1) = O(:, 1) / norm(O(:, 1));
    for c = 2:size(O, 2)
        rbase = repmat(O(:, c), 1, c - 1);
        proj = sum( O(:,1:(c-1)) .* rbase, 1 );
        O(:, c) = O(:, c) - sum( proj .* rbase, 2 );
        O(:, c) = O(:,c) / norm(O(:,c));
    end
end

