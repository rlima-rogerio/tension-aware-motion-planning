%% This function has four input parameters
% INPUT:
%   pt1: the first triangle vertex
%   pt2: the second triangle vertex
%   pt3: the third triangle vertex
%   obs: a polyshape object shape describing the obstacle
% OUTPUT:
%   isTriangleValid: a boolean variable indicating wether the triangle formed by
%            (pt1,pt2,pt3) intersects the object 'obs'. 
function isTriangleValid = checkTriangulation(pt1,pt2,pt3,obs)

pts = [pt1; pt2; pt3]';
ps = polyshape(pts(1,:),pts(2,:));

% Computes intersection between two 2-D areas
int = intersect(ps,obs);

% % DEBUG
% figure
% plot(obs); hold on;
% plot(ps);
% plot(pt2(1),pt2(2),'ko'); % plot point B


if (isempty(int.Vertices))
    isTriangleValid = false;
else
    isTriangleValid = true;
end

end