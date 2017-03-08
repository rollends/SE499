function stop = plotdriver( t, y, varargin )
%PLOTDRIVER Summary of this function goes here
%   Detailed explanation goes here
    stop = 0;
    
    if strcmp(varargin{1}, 'done')
        return;
    elseif strcmp(varargin{1}, 'init')
        
    end
    
    if size(y, 2) == 2
        plot(y(1,end), y(2,end), '.b');
    else
        plot(y(1), y(2), '.b');
    end
    drawnow;
end

