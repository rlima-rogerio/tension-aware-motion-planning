%% Computes the distance given the graph nodes
%
% INPUT:
%   G:      Graph in the 2D coordinate domain
%   nodeID: Array of node IDs connecting the start and goal nodes
%
% OUTPUT:
%   tetherLength: Length of the tether using Euclidean distance from 's' to 'g'
function tetherLength = calcTetherLength(G, nodeIDs)

coord2D = zeros(2, length(nodeIDs));
for nodeIdx = 1 : length(nodeIDs)
    coord2D(:, nodeIdx) = G.coord(nodeIDs(nodeIdx));
end

tetherLength = 0;
for coordIdx = 1 : size(coord2D,2) - 1
    x1 = coord2D(1, coordIdx);
    x2 = coord2D(1, coordIdx + 1);
    y1 = coord2D(2, coordIdx);
    y2 = coord2D(2, coordIdx + 1);
    tetherLength = tetherLength + norm([x2-x1, y2-y1]);
end
