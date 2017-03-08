function Path = rrt( x0, domain, isInGoal, isColliding )
%RRTSTAR Summary of this function goes here
%   Detailed explanation goes here
    goalReached = false;
    
    if ~isrow(x0)
        x0 = x0';
    end
    
    % Store point in set...and construct empty edge set
    V = x0(1:2);
    E = [];
    
    lambda = 0.5;
    
    while ~goalReached
        % Pick random point in domain
        xRandom = ((domain(:,2) - domain(:,1)) .* rand(2, 1) + domain(:,1))';
        
        % Find nearest point in set.
        dist = sum((V - repmat(xRandom, size(V, 1), 1)) .^ 2, 2);
        [~,sortInd] = sort(dist);
        
        I = -1;
        for ind = sortInd'
            % Find edge going towards the nearest point (should only be one!)
            if ind ~= 1
                [~, backedge] = findbackedge(E, ind);
                vector1 = V(backedge(2),:) - V(backedge(1),:);
                vector1 = vector1 / norm(vector1);
            else
                vector1 = [cos(x0(3)), sin(x0(3))];
            end
            xNearest = V(ind, :);
            vector2 = (xRandom - xNearest);
            vector2 = vector2 / norm(vector2);
            
            % Prune this vertex if it changes the angle too large.
            if dot(vector1, vector2) > 0 && ...
               abs(acos(dot(vector1, vector2))) < pi/8
                % Found an index with minimum distance that satisfies curve
                % constraints
                I = ind;
                break;
            end
        end
        if I < 0
            continue
        end
        
        %
        xCandidate = xRandom;
        sampleLine = vertcat(   linspace(xNearest(1), xRandom(1), 100),...
                                linspace(xNearest(2), xRandom(2), 100));
        if isColliding(sampleLine)
            continue
        end
        
        % Add vertex and edge to graph (V,E)
        V = vertcat(V, xCandidate);
        E = vertcat(E, [I, size(V, 1)]);
        
        % Check end condition.
        goalReached = isInGoal(xCandidate);
    end
    
    % If we reached the goal, the last point must be candidate that
    % triggered the break. Use back-edges to construct the path from the
    % root to the leaf node.
    I = size(V, 1);
    Path = xCandidate;
    while I ~= 1
        [I, ~] = findbackedge(E, I);
        Path = vertcat(V(I,:), Path);
    end
end

function [I, backedge] = findbackedge(edgelist, vertex)
    backedge = edgelist(edgelist(:,2) == vertex, :);
    I = backedge(1);
end