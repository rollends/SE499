classdef SylvesterPathFollower < SpliningRobot
    %SYLVESTERPATHFOLLOWER Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = SylvesterPathFollower(pp)
            obj = obj@SpliningRobot(pp);
        end
    end
    
    methods(Access = protected)
        function xi = transverseState(obj)
            xi =[ obj.splinelevelset; obj.lieF ];
        end

        function [lff, lgf] = lie2(obj)
            J = obj.splinejacobian();
            W = obj.splinehessian();
            [f, df] = obj.flow();
            [~, dh] = obj.obs();
            
            lff = (dh * f)'*W*(dh * f);
            lgf = dot(J, dh * df(:,3));
        end
    end
    
    methods(Access = private)
        
        function v = splinehessian(obj)
            cell = num2cell(obj.poly);
            [A1,A2,A3,A4,A5,A6] = cell{1,:};
            [B1,B2,B3,B4,B5,B6] = cell{2,:};
            x1 = obj.x(1);
            x2 = obj.x(2);

            v = hessianquintic(A1,A2,A3,A4,A5,A6,B1,B2,B3,B4,B5,B6,x1,x2);
        end

        function v = splinejacobian(obj)
            cell = num2cell(obj.poly);
            [A1,A2,A3,A4,A5,A6] = cell{1,:};
            [B1,B2,B3,B4,B5,B6] = cell{2,:};
            x1 = obj.x(1);
            x2 = obj.x(2);

            v = jacobianquintic(A1,A2,A3,A4,A5,A6,B1,B2,B3,B4,B5,B6,x1,x2);
        end
        
        function lf = lieF(obj)
            J = obj.splinejacobian;
            f = obj.flow;
            [~, dh] = obj.obs;
            lf = dot(J, dh * f);
        end

        function v = splinelevelset(obj)
            cell = num2cell(obj.poly);
            [A1,A2,A3,A4,A5,A6] = cell{1,:};
            [B1,B2,B3,B4,B5,B6] = cell{2,:};
            x1 = obj.x(1);
            x2 = obj.x(2);

            if abs(A1) < 1e-12 && abs(B1) < 1e-12
                v = levelquartic(A2,A3,A4,A5,A6,B2,B3,B4,B5,B6,x1,x2);
            else
                v = levelquintic(A1,A2,A3,A4,A5,A6,B1,B2,B3,B4,B5,B6,x1,x2);
            end
        end
    end
end

