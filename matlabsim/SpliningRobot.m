classdef(Abstract) SpliningRobot < handle
    %SPLININGROBOT Spline Following Robot Simulator
    %   Implements essential functions to simulate spline following techniques. Must be
    %   extended with a specific implementation of a path follower.
    
    properties        
        t       = 0;
        x       = zeros(3, 1);
        poly    = zeros(2, 1);
        u       = 0;
        
        eta     = zeros(1, 1);
        xi      = zeros(2, 1);
        
        lieFF   = 0;
        lieGF   = 0;
        
        gain    = [-1/0.1^2, -sqrt(2) / 0.1];
        s       = zeros(2, 1);
        
        solveT  = [];
        solveXi = [];
        solveEta= [];
        solveX  = [];
        solveS  = [];
        solveU  = [];
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
        function obj = SpliningRobot(pp)
            obj.poly = pp;
        end
        
        function v = linearController(obj)
            % Calculate v
            v = dot(obj.gain, obj.xi);
        end
        
        function dx = stepRobot( obj, t, x )
            % Save state
            obj.x = x(:);
            obj.t = t;
            
            % Find closest point
            obj.eta = obj.findClosestPoint();
            obj.s = vpolyval(obj.poly, obj.eta);
            
            if abs(obj.eta - 1) < 1e-3
                % Reached end point...don't bother integrating.
                dx = 0 * x;
                return;
            end
            
            % Calculate state along transverse directions.
            obj.xi = obj.transverseState();
            
            % We can cache the lie derivatives here.
            [obj.lieFF, obj.lieGF] = obj.lie2();
            
            % Control System
            f = obj.flow(); 
            g = [0; 0; 1];
            
            % Evaluate flow of system
            v = obj.linearController();
            obj.u = (v - obj.lieFF) / obj.lieGF;
            dx = f + g * obj.u;
            
            % Save solution
            obj.storeDataRow();
        end

    end
    
    methods(Access = protected)
        function [f, df] = flow(obj)
            xv = obj.x;
            v = 1;
            f = v * [cos(xv(3)); sin(xv(3)); 0];
            df = [  0, 0, - v * sin(xv(3));
                    0, 0, v * cos(xv(3));
                    0, 0, 0 ];
        end

        function [h, dh] = obs(obj)
            dh = [1 0 0; 0 1 0];
            h = dh * obj.x;
        end
        
        function eta = findClosestPoint(obj)
            h = obj.obs();
            
            if obj.firstRun
                obj.firstRun = false;
                
                n = 1000;
                lowBound = 0;
                upBound = 1;
                dI = 1/n;
                
                samplePoints = vpolyval(obj.poly, linspace(lowBound, upBound, n));
                error = samplePoints - repmat(h, 1, n);
                [~, I] = min(sum(error .^ 2, 1));
                obj.eta = SpliningRobot.saturate(lowBound + (I-1) * dI, lowBound, upBound);    
            end
            
            % Use previous estimate to refine our search for minimal distance point.
            % options = optimset('MaxIter', 100);
            % eta = fminsearch(@(l) norm(vpolyval(obj.poly, l) - h)^2, obj.eta, options);
            dpoly = polydiff(obj.poly, 1);

            eta = obj.eta;
            ferror = @(l) norm(h - vpolyval(obj.poly, l))^2;
            fgrad = @(l) -2*dot(h - vpolyval(obj.poly, l), vpolyval(dpoly, l));

            alph = min(1, 1 / norm(vpolyval(dpoly, eta)));
            
            for i = 1:500 	
                estimate = eta - alph * fgrad(eta);
                improvement = ferror(estimate) - ferror(eta);
                if( improvement < 0 )
                    eta = estimate;
                    alph = alph * 1.2;
                else
                    alph = alph * 0.7;
                    continue;
                end
                
                if norm(fgrad(eta)) <= 1e-6
                    break
                end
            end
        end
    end

    methods(Access = private)
        function storeDataRow(obj)
            if size(obj.solveXi, 2) > obj.NDataStored + 1
                obj.solveXi(:, (obj.NDataStored+1):(obj.NDataStored+obj.BlockSize)) = 0;
                obj.solveEta(:, (obj.NDataStored+1):(obj.NDataStored+obj.BlockSize)) = 0;
                obj.solveX(:, (obj.NDataStored+1):(obj.NDataStored+obj.BlockSize)) = 0;
                obj.solveS(:, (obj.NDataStored+1):(obj.NDataStored+obj.BlockSize)) = 0;
                obj.solveU((obj.NDataStored+1):(obj.NDataStored+obj.BlockSize)) = 0;
            end
            
            f = obj.NDataStored + 1;
            
            obj.solveT(:, f) = obj.t;
            obj.solveX(:, f) = obj.x;
            obj.solveEta(:, f) = obj.eta;
            obj.solveXi(:, f) = obj.xi;
            obj.solveS(:, f) = obj.s;
            obj.solveU(f) = obj.u;
            
            obj.NDataStored = f;
        end
    end
    
    methods(Abstract, Access = protected)
        % Self-explanatory.
        xi              = transverseState(obj)
        
        % Second lie derivatives.
        [lieFF, lieGF]  = lie2(obj)
    end
    
end

