%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   A* modified into Dijkstra (Robotics Toobox) by making h(n) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [path,cost] = dijkstra(g, vstart, vgoal)
%PGraph.Astar path finding
%
% PATH = G.Astar(V1, V2) is the lowest cost path from vertex V1 to
% vertex V2.  PATH is a list of vertices starting with V1 and ending
% V2.
%
% [PATH,C] = G.Astar(V1, V2) as above but also returns the total cost
% of traversing PATH.
%
% Notes::
% - Uses the efficient A* search algorithm.
% - The heuristic is the distance function selected in the constructor, it
%   must be  admissible, meaning that it never overestimates the actual
%   cost to get to the nearest goal node.
%
% References::
% - Correction to "A Formal Basis for the Heuristic Determination of Minimum Cost Paths".
%   Hart, P. E.; Nilsson, N. J.; Raphael, B.
%   SIGART Newsletter 37: 28-29, 1972.
%
% See also PGraph.goal, PGraph.path.


directed = true;
%opt.directed = false;

%opt = tb_optparse(opt, varargin);

% The set of vertices already evaluated.
closedSet = [];

% The set of tentative vertices to be evaluated, initially containing the start node
openSet = [vstart];
came_from(vstart) = 0;    % The map of navigated vertices.

g_score(vstart) = 0;    % Cost from start along best known path.
h_score(vstart) = 0;    % Reduces to Dijkstra
% Estimated total cost from start to goal through y.
f_score(vstart) = g_score(vstart) + h_score(vstart);

while ~isempty(openSet)
    % current := the node in openSet having the lowest f_score[] value
    [~,k] = min(f_score(openSet));
    vcurrent = openSet(k);

    if vcurrent == vgoal
        % we have arrived!
        path = [];
        p = vgoal;
        while true
            path = [p path];
            p = came_from(p);
            if p == 0
                break;
            end
        end
        if nargout > 1
            cost = f_score(vgoal);
        end
        return
    end

    %remove current from openSet
    openSet = setdiff(openSet, vcurrent);
    %add current to closedSet
    closedSet = union(closedSet, vcurrent);

    if directed
        [neighbours,costs] = g.neighbours_out(vcurrent);

    else
        [neighbours,costs] = g.neighbours(vcurrent);
    end

    for k=1:length(neighbours)

        neighbour = neighbours(k);

        if ismember(neighbour, closedSet)
            continue;
        end
        tentative_g_score = g_score(vcurrent) + costs(k);

        if ~ismember(neighbour, openSet)
            %add neighbor to openSet
            openSet = union(openSet, neighbour);
            h_score(neighbour) = 0;%g.distance(neighbour, vgoal);
            tentative_is_better = true;
        elseif tentative_g_score < g_score(neighbour)
            tentative_is_better = true;
        else
            tentative_is_better = false;
        end
        if tentative_is_better
            % we found an edge that belongs to the path

            % came_from is an array that records the path taken
            %  - length g.n
            %  -  came_from(A) = B means a path segment from A -> B
            %                         if came_from(neighbour) > 0
            %                             disp('problem')
            %                         end
            came_from(neighbour) = vcurrent;
            g_score(neighbour) = tentative_g_score;
            f_score(neighbour) = g_score(neighbour) + h_score(neighbour);
        end
    end
end
path = [];
cost = Inf;
end