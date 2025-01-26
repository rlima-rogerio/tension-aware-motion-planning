function list = partitionLine(p1,p2,len)
% function partitionedLine (p1,p2,len)
% A function that partitions a line into segments of length 'len' and
% returns all the partitioned segments in a 'nx4' matrix representing each
% coordinate (x,y) of every endpoint
%   INPUT
%       p1:     1st endpoint (x1,y1) of the original line
%       p2:     2nd endpoint (x2,y2) of the original line
%       len:    the length of the new partitioned segments
%   OUTPUT
%       pp:     (nx4) matrix representing all endpoints of the partitioned
%               line. Ex. (xp1,yp1,xp2,yp2) for every line of the matrix.
tol = 1e-3;

x = [p1(1) p2(1)];
y = [p1(2) p2(2)];

if (x(1) ~= x(2))
    ps1(1) = x(find(x==min(x),1));
    ps1(2) = y(find(x==min(x),1));
    ps2(1) = x(find(x==max(x),1));
    ps2(2) = y(find(x==max(x),1));
else
    ps1(1) = x(find(y==min(y),1));
    ps1(2) = y(find(y==min(y),1));
    ps2(1) = x(find(y==max(y),1));
    ps2(2) = y(find(y==max(y),1));
end

span_x = ps2(1) - ps1(1);
span_y = ps2(2) - ps1(2);

dx = ps2(1) - ps1(1);
dy = ps2(2) - ps1(2);

theta = atan2(dy,dx);

% plot([ps1(1) ps2(1)], [ps1(2) ps2(2)],'k'); hold on; grid on;

idx = 1;
paux = ps1;
if (norm([dx dy]) < len) % Partitioning is impossible!)
    pp = [p1 p2];
elseif (p1(1) == p2(1)) % same x
    aux_y = [min(y):len:max(y)]';
    aux_x = p1(1)*ones(length(aux_y),1);
    dummy = zeros(size(aux_x));
    pp = [dummy(2:end), dummy(2:end), aux_x(2:end), aux_y(2:end)];
elseif (p1(2) == p2(2)) % same y
    aux_x = [min(x):len:max(x)]';
    aux_y = p1(2)*ones(length(aux_x),1);
    dummy = zeros(size(aux_x));
    pp = [dummy(2:end), dummy(2:end), aux_x(2:end), aux_y(2:end)];
else
    while(1)
        pp(idx,1) = paux(1); % initial x
        pp(idx,2) = paux(2); % initial y
        pp(idx,3) = pp(idx,1) + len*cos(theta);
        pp(idx,4) = pp(idx,2) + len*sin(theta);

        paux = pp(idx,3:4);
        
        % Checks if the 2nd endpoint (right-most) extrapolates the original
        % limit. If so, terminates the partitioning and put the original
        % endpoint as the last coordinate.
        if (span_x >= span_y)
            if ( (pp(idx,3) > ps2(1)) )     %|| (pp(idx,4) > ps2(2) ) 
                pp(idx,3) = ps2(1);
                pp(idx,4) = ps2(2);
                break;
            end
        else
            if ( pp(idx,4) > ps2(2) ) 
                pp(idx,3) = ps2(1);
                pp(idx,4) = ps2(2);
                break;
            end
        end
        idx = idx + 1;
    end % while
end % else

list = [ps1(1) ps1(2);pp(:,3:4)];