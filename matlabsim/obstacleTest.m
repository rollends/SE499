function collision = obstacleTest( obstacles, points )
%OBSTACLETEST Summary of this function goes here
%   Detailed explanation goes here
    collision = false;
    
    for p = points
        if any(sum((obstacles - repmat(p', size(obstacles,1), 1)) .^ 2, 2) <= 4)
            collision = true;
            return;
        end
    end
end

