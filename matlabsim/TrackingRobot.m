classdef TrackingRobot < handle
    %TRACKINGROBOT Implements a simple look-ahead tracking controller.
    %   Detailed explanation goes here
    
    properties        
        t       = 0;
        x       = zeros(3, 1);
        poly    = zeros(2, 1);
        s       = zeros(2, 1);
        ds      = zeros(2, 1);
        eta     = 0;
        
        gain    = [-10, -1; 1, -10];
        
        solveT  = [];
        solveS  = [];
        solveX  = [];
    end

    properties(Access = private)
        firstRun    = true;
        BlockSize   = 500;
        NDataStored = 0;
    end
    
    methods(Static, Access = protected)
        function s = saturate(v, low, up)
            s = min(max(v, low), up);
        end
    end
    
    methods
        function obj = TrackingRobot(pp)
            obj.poly = pp;
        end
        
        function v = linearController(obj)
            % Calculate v
            v = dot(obj.gain, obj.x(1:2) - obj.s);
        end
        
        function dx = stepRobot( obj, t, x )
            % Save state
            obj.x = x(:);
            obj.t = t;
            
            % Find closest point
            obj.s = vpolyval(obj.poly, t);
            obj.ds = vpolyval(polydiff(obj.poly, 1), t);
            
            % Control System
            f = obj.flow(1);
            [~,dh] = obj.obs();
            g = [0; 0; 1];
            
            l = 0.05;
            tau = dh * f;
            R = [0 -1; 1 0];
            beta = [tau, l * R * tau];
            
            ya = x(1:2) + l * tau;
            
            if abs(t - 1) < 1e-3
                % Reached end point...don't bother integrating.
                dx = 0 * x;
                return;
            end
            
            % Evaluate flow of system
            u = obj.ds + obj.gain * (ya - obj.s);
            control = beta \ u;
            
            dx = obj.flow(control(1)) + g * control(2);
            
            % Save solution
            obj.storeDataRow();
        end
    end
    
    methods(Access = protected)
        function [f, df] = flow(obj, v)
            xv = obj.x;
            f = v * [cos(xv(3)); sin(xv(3)); 0];
            df = [  0, 0, - v * sin(xv(3));
                    0, 0, v * cos(xv(3));
                    0, 0, 0 ];
        end

        function [h, dh] = obs(obj)
            dh = [1 0 0; 0 1 0];
            h = dh * obj.x;
        end
    end
    
    methods(Access = private)
        function storeDataRow(obj)
            if size(obj.solveX, 2) > obj.NDataStored + 1
                obj.solveX(:, (obj.NDataStored+1):(obj.NDataStored+obj.BlockSize)) = 0;
                obj.solveS(:, (obj.NDataStored+1):(obj.NDataStored+obj.BlockSize)) = 0;
            end
            
            f = obj.NDataStored + 1;
            
            obj.solveT(:, f) = obj.t;
            obj.solveX(:, f) = obj.x;
            obj.solveS(:, f) = obj.s;
            
            obj.NDataStored = f;
        end
    end
end

