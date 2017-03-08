function controller = polysfcontrol( pp, varargin )
%POLYSFCONTROL Polynomial Path Following Controller
%   This path following controller uses Seret-Frenet Frames. The
%   controller, after cancelling out non-linearities, implements
%   an LQR optimal controller.
    assert( 2 == size(pp, 1), '%s', ...
            'The polynomial (pp) is assumed to result in a 2-vector.' );
    assert( 3 <= size(pp, 2), '%s', ...
            'The polynomial (pp) should atleast be of 3rd order.' );
    controller = @(x) polycontroller(pp, x);
end

function u = polycontroller(pp, x)
    % Columnize
    x = x(:);

    % Find the closest parameter eta.
    de = 1e-1;
    eta = besteta(pp, de, x);
    
    % Normal/Tangent vectors and perform gram-schmidt
    e = polypathbasis(pp, eta);
    
    % Find coordinate in linearized space.
    vxi = xi(pp, x, e, eta);

    % Optimal Gain (calculated via LQR with Q = I, R = 1)
    K = [1 sqrt(3)];
         
    % Calculate v optimal control
    v = -dot(K, vxi); 

    % Calculate 2nd Order Lie Derivatives (that is LieFF and LieGF)
    [LFF, LGF] = lie2(pp, x, e, eta);
    
    % Calculate the control signal required for the non-linear system. Should think of a
    % smart way to deal with catastrophic cancelling here.
    u = (v - LFF) / LGF;
    
    % Make sure to fail if we ended up generating a really large signal,
    % maybe lieGF was too small or lieFF is too big.
    assert( norm(u, Inf) < 1e10, ...
            'Control Signal has become too large to handle. \n\t LieGF = %g. \n\t LieFF = %g', ...
            LGF,...
            LFF );
end

function e = polypathbasis(pp, eta)
    % Tangent Vector
    t = vpolyval(polydiff(pp, 1), eta);
    
    % Normalize
    t = t / norm(t);
    
    % Rotate vector to form full orthonormal basis.
    t2 = [0 -1; 1 0] * t;
    
    % Basis
    e = [t, t2];
end

function [f, df] = flow(x)
    f = [cos(x(3)); sin(x(3)); 0];
    df = [  0, 0, -sin(x(3));
            0, 0, cos(x(3));
            0, 0, 0 ];
end

function [h, dh] = obs(x)
    dh = [1 0 0; 0 1 0];
    h = dh * x;
end

function xi = xi(pp, x, e, eta)
    [f, ~] = flow(x);
    [h, dh] = obs(x);
    g = [0 ; 0 ; 1];
    
    % Polynomial value and derivatives at closest point.
    s = vpolyval(pp, eta);
    ds = vpolyval(polydiff(pp, 1), eta);
    dds = vpolyval(polydiff(pp, 2), eta);
    
    % Helper derivatives
    dwdh = e(:,1)' / norm(ds);
    dsdh = e(:,1) * e(:,1)';
    dhdt = dh * f;
    
    % Basis Vector Derivatives
    R = [0 -1; 1 0];
    de1dl = dds / norm(ds) - dot(dds, ds)/dot(ds,ds) * e(:,1);
    de1dx = de1dl * dwdh * dh;
    de2dl = R * de1dl;
    de2dx = de2dl * dwdh * dh;

    xi1 = dot(h - s, e(:, 2));
    xi2 = dot((eye(2) - dsdh) * dhdt, e(:,2)) + dot(h - s, de2dx * f);
    
    xi = vertcat(xi1, xi2);
end

function [lff, lgf] = lie2(pp, x, e, eta)
    [f, df] = flow(x);
    [h, dh] = obs(x);
    g = [0 ; 0 ; 1];
    
    % Polynomial value and derivatives at closest point.
    s = vpolyval(pp, eta);
    ds = vpolyval(polydiff(pp, 1), eta);
    dds = vpolyval(polydiff(pp, 2), eta);
    ddds = vpolyval(polydiff(pp, 3), eta);
    
    % Helper derivatives
    dwdh = e(:,1)' / norm(ds);
    dsdh = e(:,1) * e(:,1)';
    dhdt = dh * f;
    
    % Basis Vector Derivatives
    R = [0 -1; 1 0];
    de1dl = dds / norm(ds) - dot(dds, ds)/dot(ds,ds) * e(:,1);
    de1dx = de1dl * dwdh * dh;
    de2dl = R * de1dl;
    de2dx = de2dl * dwdh * dh;

    % dXi2dx = d/dx [ xi2 ]
    % xi2 = dot((eye(2) - dsdh) * dhdt, e(:,2)) +
    %       dot(h - s, de2dx * f);
    dXi2dx = zeros(1, length(x));
    
    for xi = 1:length(dXi2dx)
        % Part A : d/dx [ dot((eye(2) - dsdh) * dhdt, e(:,2)) ]
        
            % dot((eye(2) - dsdh) * dhdt, d/dx [ e(:,2) ] )
            A1 = dot((eye(2) - dsdh) * dhdt, de2dx(:, xi));

            % dot( d/dx [ (eye(2) - dsdh) * dhdt ], e(:,2) )
            A2 = dot(-(de1dx(:,xi)*e(:,1)' + e(:,1)*de1dx(:,xi)') * dhdt + (eye(2) - dsdh) * dh * df(:, xi), e(:,2));
        
        % Part B : d/dx [ dot( h - s, de2dx * f ) ]
        
            % dot( d/dx [ h - s ], de2dx * f )
            B1 = dot( dh(:,xi) - dsdh * dh(:,xi), de2dx * f );

            % Helper Second Derivatives!
            ddwdhdx = ( norm(ds) * de1dx(:, xi)' - e(:,1)' / norm(ds) * dot(ds, dds * dwdh * dh(:, xi)) ) / dot(ds, ds);
            dde1dldx = ( norm(ds) * ddds * dwdh * dh(:,xi) - dds / dot(ds,ds) * dot(ds, dds * dwdh * dh(:,xi)) ) / dot(ds, ds) - ...
                       dot(dds,ds)/dot(ds,ds) * de1dx(:,xi) - ... dot(dds, ds)/dot(ds,ds) * d/dx [e(:,1)]
                       (    dot(ds,ds) * (dot(ddds * dwdh * dh(:,xi), ds) + dot(dds, dds * dwdh * dh(:,xi))) - ... d/dx [ dot(dds, ds)/dot(ds,ds) ] e(:,1)
                            2 * dot(dds, ds) * dot(ds * dwdh * dh(:,xi), ds) ) ...
                        / (dot(ds,ds)^2) * e(:,1);
            dde1ddx = ( dde1dldx * dwdh + de1dl * ddwdhdx ) * dh;
            dde2ddx = R * dde1ddx;

            % dot( h - s, d/dx[ de2dx * f ] )
            B2 = dot( h - s, dde2ddx * f + de2dx * df(:, xi) );

        dXi2dx(xi) = A1 + A2 + B1 + B2;
    end
    
    lff = dXi2dx * f;
    lgf = dXi2dx * g;
end

function s = saturate(v, low, up)
    s = min(max(v, low), up);
end

function best = besteta(coefs, ~, x)
    global eta;

    % Find the closest position on the WHOLE polynomial.
    if eta < 0
        % Initial pass, attempt a full estimate by a broad search.
        n = 1000;
        lowBound = 0;
        upBound = 1;
        dI = 1/n;
%     else
%         % Attempt a small pertubation of eta, and see if we can find a
%         % better estimate of eta.
%         n = 100;
%         lowBound = eta - de;
%         upBound = eta + de;
%         dI = (2*de)/n;
        samplePoints = vpolyval(coefs, linspace(lowBound, upBound, n));
        error = samplePoints - repmat(x(1:2), 1, n);
        [~, I] = min(sum(error .^ 2, 1));
        eta = saturate(lowBound + (I-1) * dI, lowBound, upBound);    
    end
    
    % Use previous estimate to refine our search for minimal distance point.
    options = optimset('MaxIter', 20);
    h = obs(x);
    eta = fminsearch(@(l) norm(vpolyval(coefs, l) - h), eta, options);
    best = eta;
end
