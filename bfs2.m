%% Breath-First-Search (BFS) implementation
function [path]= bfs2(g, goal)

C = [];                         % Closed set
O = CQueue();                   % Open set
v = 1;                          % Start node
p = [1 0];                      % Array to save back-pointers

%O = fifo('enqueue',O,v,0);      % Add v to set O
O.push(v)                   % Add v to set O
while (~O.isempty())
    u = O.pop();                % Remove node 'u' from 'O'
    if (isempty(find(C == u, 1)))  % If 'u' is not an element of 'C'
        C = [C u];              % Add 'u' to 'C'
        nb = g.neighbours_out(u);   % Get all neighbors of 'u'
                                % Add all nb of u that r not in C to O
        %nb = unique(nb);
        for n=1:length(nb)
            isIn = find(C == nb(n), 1);
            if(isempty(isIn))
                % O = fifo('enqueue',O,nb(n),u);
                O.push(nb(n));
                p = [p; nb(n) u];
                g.setvdata(nb(n), u);    % Saves the back-pointer
            end % end-if
        end %end-for
    end % end-if
end % end-while


visited_nodes = [p(:,1)];
[b,m1,n1] = unique(visited_nodes,'first');
[c1,d1] =sort(m1);
% c1 = c1(1:end-1); % Remove last line from the list
gs = p(c1,:); % Graph searched: [node, weight, back-pointer]

bp = goal;
path = [goal];
idx = find(gs(:,1) == goal);
k = 1;
while(bp ~= 1) % Node #1 is the starting point
    bp = gs(idx,2);
    path = [path bp];
    idx = find(gs(:,1) == bp);
    
    if (k > g.n)        % If all nodes were visited and did not reach the
        path = -1;      % start point (node 1), then, a valid path there is
        break;          % no exists!
    end                 
    k = k + 1;
end

path = path(end:-1:1);
