classdef FramePathFollower < SpliningRobot
    %FRAMEPATHFOLLOWER Path follower using tangent/normal vector projection.
    %   
    
    properties
        basis;
    end
    
    properties
        s;
        ds;
        dds;
        ddds;
        curveSpeed;
        
        dwdh;
        dsdh;
        dedl;
    end
    
    methods
        function obj = FramePathFollower(pp)
            obj = obj@SpliningRobot(pp);
            %obj.gain = 10 * obj.gain;
        end
    end
        
    methods(Access = protected)       
        function xi = transverseState(obj)
            obj.cacheCurveState();
            obj.cacheFrameState();
            
            e = obj.basis;
            
            % System 
            [f, ~] = obj.flow();
            [h, dh] = obj.obs();

            % Polynomial value and derivatives at closest point.
            s = obj.s;
            
            % Helper derivatives
            xi1 = dot(h - s, e(:, 2));
            %xi2 = dot((eye(2) - dsdh) * dhdt, e(:,2)) + dot(h - s, de2dx * f);
            xi2 = dot(dh * f, e(:,2));

            xi = vertcat(xi1, xi2);
        end

        function [lff, lgf] = lie2(obj)  
            % System
            [f, df] = obj.flow();
            [~, dh] = obj.obs();
            g = [0 ; 0 ; 1];

            % Polynomial value and derivatives at closest point.
            s = obj.s; % vpolyval(obj.poly, obj.eta);
            ds = obj.ds; % vpolyval(polydiff(obj.poly, 1), obj.eta);
            dds = obj.dds; % vpolyval(polydiff(obj.poly, 2), obj.eta);
            ddds = obj.ddds; % vpolyval(polydiff(obj.poly, 3), obj.eta);
            cSpeed = obj.curveSpeed;
            
            % Helper derivatives
            e = obj.basis;
            dwdh = obj.dwdh;

            % Basis Vector Derivatives
            R = [0 -1; 1 0];
            de2dx = obj.dedl(:,2) * dwdh * dh;

            lff = dot( de2dx * f, dh * f );
            lgf = dot( e(:,2), dh * df * g );
            
            % dXi2dx = d/dx [ xi2 ]
            % xi2 = dot((eye(2) - dsdh) * dhdt, e(:,2)) +
            %       dot(h - s, de2dx * f);
%             dXi2dx = zeros(1, length(obj.x));
% 
%             for xi = 1:length(dXi2dx)
%                 % Part A : d/dx [ dot((eye(2) - dsdh) * dhdt, e(:,2)) ]
% 
%                     % dot((eye(2) - dsdh) * dhdt, d/dx [ e(:,2) ] )
%                     A1 = dot((eye(2) - dsdh) * dhdt, de2dx(:, xi));
% 
%                     % dot( d/dx [ (eye(2) - dsdh) * dhdt ], e(:,2) )
%                     A2 = dot(-(de1dx(:,xi)*e(:,1)' + e(:,1)*de1dx(:,xi)') * dhdt + (eye(2) - dsdh) * dh * df(:, xi), e(:,2));
% 
%                 % Part B : d/dx [ dot( h - s, de2dx * f ) ]
% 
%                     % dot( d/dx [ h - s ], de2dx * f )
%                     B1 = dot( dh(:,xi) - dsdh * dh(:,xi), de2dx * f );
% 
%                     % Helper Second Derivatives!
%                     ddwdhdx = ( cSpeed * de1dx(:, xi)' - e(:,1)' / cSpeed * dot(ds, dds * dwdh * dh(:, xi)) ) / (cSpeed^2);
%                     dde1dldx = ( cSpeed * ddds * dwdh * dh(:,xi) - dds / (cSpeed^2) * dot(ds, dds * dwdh * dh(:,xi)) ) / (cSpeed^2) - ...
%                                dot(dds,ds)/(cSpeed^2) * de1dx(:,xi) - ... dot(dds, ds)/dot(ds,ds) * d/dx [e(:,1)]
%                                (    (cSpeed^2) * (dot(ddds * dwdh * dh(:,xi), ds) + dot(dds, dds * dwdh * dh(:,xi))) - ... d/dx [ dot(dds, ds)/dot(ds,ds) ] e(:,1)
%                                     2 * dot(dds, ds) * dot(ds * dwdh * dh(:,xi), ds) ) ...
%                                 / (cSpeed^4) * e(:,1);
%                     dde1ddx = ( dde1dldx * dwdh + obj.dedl(:,1) * ddwdhdx ) * dh;
%                     dde2ddx = R * dde1ddx;
% 
%                     % dot( h - s, d/dx[ de2dx * f ] )
%                     B2 = dot( h - s, dde2ddx * f + de2dx * df(:, xi) );
% 
%                 dXi2dx(xi) = A1 + A2 + B1 + B2;
%             end
% 
%             lff = dXi2dx * f;
%             lgf = dXi2dx * g;
        end
    end
    
    methods(Access = private)
        function cacheCurveState(obj)
            obj.s = vpolyval(obj.poly, obj.eta);
            obj.ds = vpolyval(polydiff(obj.poly, 1), obj.eta);
            obj.dds = vpolyval(polydiff(obj.poly, 2), obj.eta);
            obj.ddds = vpolyval(polydiff(obj.poly, 3), obj.eta);
            obj.curveSpeed = norm(obj.ds);
        end
        
        function cacheFrameState(obj)
            obj.basis = obj.polypathbasis();
            e = obj.basis;
            
            % Helper derivatives
            obj.dwdh = e(:,1)' / obj.curveSpeed;
            obj.dsdh = e(:,1) * e(:,1)';

            % Basis Vector Derivatives
            R = [0 -1; 1 0];
            de1dl = obj.dds / obj.curveSpeed - dot(obj.dds, obj.ds) / (obj.curveSpeed^2) * e(:,1);
            de2dl = R * de1dl;
            
            obj.dedl = [de1dl, de2dl];
        end
        
        function e = polypathbasis(obj)
            % Tangent Vector
            t = vpolyval(polydiff(obj.poly, 1), obj.eta);

            % Normalize
            t = t / norm(t);

            % Rotate vector to form full orthonormal basis.
            t2 = [0 -1; 1 0] * t;

            % Basis
            e = [t, t2];
        end
    end
    
end

