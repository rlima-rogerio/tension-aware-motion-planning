function retVal = isFree(q)
% function isFree(q)
% A function that checks if a chosen point inside a C-space is (retVal):
%   0 (false):  belongs an obstacle
%   1 (true):   is a free C-space
% 	q is the configuration
% retVal 

% Global variables
global vectorObs;

xq = q(1);
yq = q(2);

% Initialize C-space as FREE
in = true; 

for k=1:length(vectorObs)
    xobs = vectorObs{k}.Vertices(:,1);
    yobs = vectorObs{k}.Vertices(:,2);
	[inside,onedge] = inpolygon(xq,yq,xobs,yobs);
    
    % If an obstacle is detected, exit the loop and returns FALSE
    if (inside && ~onedge)
        in = false;
        break;
    end
end


retVal = in;
