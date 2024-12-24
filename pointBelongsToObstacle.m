%% This function has inputs a point and a cell array where each cell represents 
%  an obstacle as a polyshape object. This function returns the object
%  containing the point.
function [obs, idx] = pointBelongsToObstacle(point, allObs)

for m=1:length(allObs) % First obstacle is the environment boundary
    minimums(m) = min(vecnorm([allObs{m}.Vertices - point]'));
end
[min_,idx] = min(minimums);
obs = allObs(idx);

end