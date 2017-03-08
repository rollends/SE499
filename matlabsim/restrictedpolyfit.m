function poly = restrictedpolyfit( order, delta, controlPoints, varargin )
%RESTRICTEDPOLYFIT Restricted Polynomial Fit
%   Generates a polynomial of the order requested such that it fits closely to the control
%   points within a 1-norm error specified by delta. This function assumes a linear change
%   between control points and attempts to fit a polynomial as closely to these lines.
    ip = inputParser;
    addRequired(ip, 'order', @isnumeric);
    addRequired(ip, 'delta', @isnumeric);
    addRequired(ip, 'controlPoints', @isnumeric);
    addParameter(ip, 'upperBound', NaN);
    addParameter(ip, 'lowerBound', NaN);
    parse(ip, order, delta, controlPoints, varargin{:});
    
    order = ip.Results.order;
    delta = ip.Results.delta;
    controlPoints = ip.Results.controlPoints;
    ub = [];
    lb = [];
    
    if ~isnan(ip.Results.upperBound)
        ub = ip.Results.upperBound * ones(1, (order+1)*2);
    end
    if ~isnan(ip.Results.lowerBound)
        lb = ip.Results.lowerBound * ones(1, (order+1)*2);
    end

    % Resample between control points linearly to help us generate a polynomial that 
    % tightly fits a line between waypoints.
    Y = generateIntermediatePoints( controlPoints );
    
    % Construct standard Least Squares matrix (i.e. C such that x = pinv(C)*d is the least
    % squares solution to the standard poly fit problem)
    [C, d] = constructLSQMatrix( Y, order );
    
    % Construct (In)Equality constraint matrices
    [A, b, Aeq, beq] = constructConstraints( controlPoints, order, delta );
    
    % Optimization Options (new Interior Point algorithm)
    option = optimoptions(  @lsqlin, ...
                            'Algorithm', 'interior-point', ...
                            'MaxIterations', 1000, ...
                            'ConstraintTolerance', delta * 1e-1);
    
    % Least Squares with Linear Constraints (Equality and Inequality)
    [p2,~] = lsqlin(double(C), double(d), ... Least Squares System
                    double(A), double(b), ... Linear Constraints (1-Norm Ball around Control Points)
                    double(Aeq), double(beq), ... Equality Constraints
                    lb, ... Lower Bound
                    ub, ... Upper Bound
                	[], ... Initializer (doesn't work for the more recent algos)
                    option);
    
    % Reshape coefficient vector into 'pp' form.
    poly = transpose(reshape(p2, order + 1, 2));
end

function Yc = generateIntermediatePoints( controls )
    Yc = controls(:,1);
    SampleRate = 20;
    for i = 1:(size(controls, 2)-1)
        samples = vertcat( linspace(controls(1,i),controls(1,i+1),SampleRate),...
                           linspace(controls(2,i),controls(2,i+1),SampleRate) );
        Yc = horzcat(Yc, samples(:,2:end));
    end
end

function [A, b, Aeq, beq] = constructConstraints( Yc, N, delta )
    A = zeros(2 * numel(Yc), (N+1) * 2);
    b = zeros(2 * numel(Yc), 1);

    index = 1;
    
    % Go through every coordinate
    for xi = 1:size(Yc, 1)
        % Go through every point
        for pi = 1:size(Yc, 2)
            base = (xi-1)*(N+1);
            lambda = (pi - 1) / (size(Yc,2) - 1);
            polyLambda = lambda .^ (N:-1:0);

            % Enforce upper bound
            A(index, (base+1):(base+(N+1))) = polyLambda;
            b(index) = Yc(xi, pi) + delta;

            % Next row of C
            index = index + 1;

            % Enforce lower bound
            A(index, (base+1):(base+(N+1))) = -polyLambda;
            b(index) = -Yc(xi, pi) + delta;

            % Next row of C
            index = index + 1;
        end
    end
    
    Aeq = zeros(4, (N+1) * 2);
    beq = zeros(4, 1);

    Aeq(1, 1:(N+1)) = circshift(eye(1, N+1), -1);
    Aeq(2, (N+2):end) = circshift(eye(1, N+1), -1);
    beq(1:2) = Yc(:,1);

    Aeq(3, 1:(N+1)) = ones(1, N+1);
    Aeq(4, (N+2):end) = ones(1, N+1);
    beq(3:4) = Yc(:,end);
end

function [C, d] = constructLSQMatrix( Y, N )
    C = zeros(numel(Y), (N+1) * 2);
    d = zeros(numel(Y), 1);

    index = 1;
    % Go through every coordinate
    for xi = 1:size(Y, 1)
        % Go through every point
        for pi = 1:size(Y, 2)
            base = (xi-1)*(N+1);
            lambda = (pi - 1) / (size(Y,2) - 1);
            polyLambda = lambda .^ (N:-1:0);

            C(index, (base+1):(base+(N+1))) = polyLambda;
            d(index) = Y(xi, pi);

            % Next row of C
            index = index + 1;
        end    
    end
end