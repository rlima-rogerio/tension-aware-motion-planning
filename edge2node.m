%% This function converts a visibility list into a node by grouping 
%  pairwise visible vertex (from the visibility table).
% INPUT:
%   vl: visilibity list in a table format
% OUTPUT:
%   nodes: cell array containing pairs of visible vertices as a single node
%
% Example:
% 
%   Vertex #    |   Visibility List
%   ----------------------------------
%       vs      |   1
%       1*      |   [2, 4]
%       2       |   [1, 3]
%       3       |   [2, 5]
%       4       |   [1, 5]
%       5^      |   [3, 4] 
%
% Legend:
%   vs: Virtual start point
%    *: Start
%    ^: Goal
function nodes = edge2node(vl)
% # of vertices in the environment, including the starting and goal points
num_vertices = size(vl,1); 
nodes = [];
neighborList = [];
nn = 1;
    for k = 1 : num_vertices
        neighborList = find(vl(k,:));
        for neighbor = neighborList
            nodes{nn} = [k,neighbor];
            nn = nn + 1;
        end
    end
end
