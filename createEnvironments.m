%% Create sets of enviroments
close all
clear
clc

SAVE_OBSTACLES = false;
PLOT_OBSTACLE_ID = true;
setID = 0;

%% Transformations
% Rotation
R = @(theta)([cos(theta) -sin(theta)    0; ...
              sin(theta) cos(theta)     0; ...
              0          0              1]);

% Translation
T = @(d)([1 0 d(1); ...
          0 1 d(2); ...
          0 0 1]);

% Scale
S = @(s)([s 0 0; 0 s 0; 0 0 1]);

%% Shapes
% Square (L = 1) or rectangle (L ~= 1)
Rectangle =@(L)([0 L L 0; ...
                 0 0 1 1; ...
                 1 1 1 1]);

% Hexagon
c60 = cos(pi/3); s60 = sin(pi/3);
Hexagon =([0, 1, 1+c60, 1,      0       -c60; ...
           0, 0, s60,   2*s60,  2*s60   s60; ...
           1, 1, 1,     1,      1,      1]);



%% Set 1
load Boundary.mat

setID = setID + 1;
obsIdx = 1;

obsIdx = obsIdx + 1;
rot = pi/2;     trans = [5.9, 6];     scale = 0.45;       rectWidth = 12;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = pi/2;     trans = [4.1, 3];     scale = 0.45;       rectWidth = 10;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = 0;     trans = [6.5, 2];     scale = 0.45;       rectWidth = 6;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% Plot obstacles
figure;
for k = 1 : length(vectorObs)
    plot(vectorObs{k}); hold on;
    
    if PLOT_OBSTACLE_ID
        if k > 1
            % Computes the centroid to plot the obstacle numbers
            [xc, yc] = centroid(vectorObs{k});
            text(xc, yc, sprintf('%d', k - 1))
        end
    end
end
pbaspect([1 1 1]);

if SAVE_OBSTACLES
    setName = sprintf('Obstacles_%dSet%d', length(vectorObs) - 1, setID);
    %setName = sprintf('Obstacles_Set%d', setID);
    save(setName,'vectorObs');
end

clearvars vectorObs 

%% Set 2
load Boundary.mat

setID = setID + 1;
obsIdx = 1;

obsIdx = obsIdx + 1;
rot = pi/2;     trans = [1.5, 7];     scale = 1.0;    rectWidth = 1.5;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = pi/2;     trans = [1.5, 2.8];     scale = 1.0;    rectWidth = 2.5;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = pi/2;     trans = [3.2, 1.5];     scale = 0.8;       rectWidth = 2;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = pi/3;     trans = [3.5, 7];     scale = 0.5;       rectWidth = 5;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = pi/7;     trans = [5.5, 5.3];     scale = 0.5;       rectWidth = 11;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = pi/2;     trans = [8, 2.5];     scale = 0.5;       rectWidth = 7;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = 0;     trans = [8.5, 5.3];     scale = 0.5;       rectWidth = 4;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = 0;     trans = [5.5, 3];     scale = 1.;    rectWidth = 2;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = pi/2;     trans = [7.5, 8];     scale = 0.6;
shape = Hexagon;
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% Plot obstacles
figure;
for k = 1 : length(vectorObs)
    plot(vectorObs{k}); hold on;
    
    if PLOT_OBSTACLE_ID
        if k > 1
            % Computes the centroid to plot the obstacle numbers
            [xc, yc] = centroid(vectorObs{k});
            text(xc, yc, sprintf('%d', k - 1), "HorizontalAlignment", "center")
        end
    end
end
pbaspect([1 1 1]);

if SAVE_OBSTACLES
    setName = sprintf('Obstacles_%dSet%d', length(vectorObs) - 1, setID);
    save(setName,'vectorObs');
end

%% Set 3
load Boundary.mat

setID = setID + 1;
obsIdx = 1;

obsIdx = obsIdx + 1;
rot = pi/4;     trans = [2., 5.];     scale = 0.3;       rectWidth = 8;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = pi/4;     trans = [4.2, 7.2];     scale = 0.3;       rectWidth = 8;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = pi/4;     trans = [3, 3];     scale = 0.3;       rectWidth = 12;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = pi/4;     trans = [6, 6];     scale = 0.3;       rectWidth = 12;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = pi/4;     trans = [5.8, 2.8];     scale = 0.3;       rectWidth = 8;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = pi/4;     trans = [8, 5];     scale = 0.3;       rectWidth = 8;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% Top-right rectangular obstacle
obsIdx = obsIdx + 1;
rot = 0;     trans = [8.5, 8.6];     scale = 0.5;       rectWidth = 2;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% Plot obstacles
figure;
for k = 1 : length(vectorObs)
    plot(vectorObs{k}); hold on;
    
    if PLOT_OBSTACLE_ID
        if k > 1
            % Computes the centroid to plot the obstacle numbers
            [xc, yc] = centroid(vectorObs{k});
            text(xc, yc, sprintf('%d', k - 1))
        end
    end
end
pbaspect([1 1 1]);

if SAVE_OBSTACLES
    setName = sprintf('Obstacles_%dSet%d', length(vectorObs) - 1, setID);
    save(setName,'vectorObs');
end


%% Set 4
load Boundary.mat

setID = setID + 1;
obsIdx = 1;

obsIdx = obsIdx + 1;
rot = pi/2;     trans = [5, 5.];     scale = 3;       rectWidth = 2;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% Plot obstacles
figure;
for k = 1 : length(vectorObs)
    plot(vectorObs{k}); hold on;
    
    if PLOT_OBSTACLE_ID
        if k > 1
            % Computes the centroid to plot the obstacle numbers
            [xc, yc] = centroid(vectorObs{k});
            text(xc, yc, sprintf('%d', k - 1))
        end
    end
end
pbaspect([1 1 1]);

if SAVE_OBSTACLES
    setName = sprintf('Obstacles_%dSet%d', length(vectorObs) - 1, setID);
    save(setName,'vectorObs');
end

%% Set 5
load Boundary.mat

setID = setID + 1;
obsIdx = 1;

obsIdx = obsIdx + 1;
rot = pi/2;     trans = [5, 5.];     scale = 1;       rectWidth = 2;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = 0;     trans = [5.5, 3];     scale = 1.8;    rectWidth = 1.5;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = pi/2;     trans = [5.5, 8];     scale = 1;
shape = Hexagon;
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% Plot obstacles
figure;
for k = 1 : length(vectorObs)
    plot(vectorObs{k}); hold on;
    
    if PLOT_OBSTACLE_ID
        if k > 1
            % Computes the centroid to plot the obstacle numbers
            [xc, yc] = centroid(vectorObs{k});
            text(xc, yc, sprintf('%d', k - 1))
        end
    end
end
pbaspect([1 1 1]);

if SAVE_OBSTACLES
    setName = sprintf('Obstacles_%dSet%d', length(vectorObs) - 1, setID);
    save(setName,'vectorObs');
end
%% Set 6(R)
load BoundaryW20H10.mat

setID = setID + 1;
obsIdx = 1;

% 01
obsIdx = obsIdx + 1;
rot = pi/2;     trans = [1.5, 6.8];     scale = 0.7;    rectWidth = 2;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 02
obsIdx = obsIdx + 1;
rot = pi/2;     trans = [1.5, 2.5];     scale = 1.0;    rectWidth = 2.5;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 03
obsIdx = obsIdx + 1;
rot = pi/2;     trans = [4, 1.5];     scale = 0.8;       rectWidth = 2;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 04
obsIdx = obsIdx + 1;
rot = pi/2;     trans = [15., 5.0];     scale = 0.5;       rectWidth = 10;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 05
obsIdx = obsIdx + 1;
rot = pi/2;     trans = [4, 6];     scale = 0.5;       rectWidth = 11;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 06
obsIdx = obsIdx + 1;
rot = pi/2;     trans = [10.2, 6];     scale = 0.5;       rectWidth = 7;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 07
obsIdx = obsIdx + 1;
rot = pi/2;     trans = [17.0, 7.8];     scale = 0.3;       rectWidth = 8;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 08
obsIdx = obsIdx + 1;
rot = 0;     trans = [17.0, 5.3];     scale = 0.5;       rectWidth = 4;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 09
obsIdx = obsIdx + 1;
rot = 0;     trans = [12.5, 4];     scale = 1.;    rectWidth = 2;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 10
obsIdx = obsIdx + 1;
rot = pi/2;     trans = [7.0, 8.5];     scale = 0.5;
shape = Hexagon;
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 11
obsIdx = obsIdx + 1;
rot = 0;     trans = [1.5, 5.0];     scale = 0.5;       rectWidth = 4;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');


% 12
obsIdx = obsIdx + 1;
rot = pi/2;     trans = [1.5, 8.5];     scale = 0.8;    rectWidth = 1;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 13
obsIdx = obsIdx + 1;
rot = pi/2;     trans = [5.5, 4.5];     scale = 0.5;
shape = Hexagon;
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 14
obsIdx = obsIdx + 1;
rot = pi/3;     trans = [10, 1.5];     scale = 0.7;       rectWidth = 3;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 15
obsIdx = obsIdx + 1;
rot = pi/3;     trans = [13, 1.5];     scale = 0.7;       rectWidth = 3;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 16
obsIdx = obsIdx + 1;
rot = pi/3;     trans = [16, 1.5];     scale = 0.7;       rectWidth = 3;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 17
obsIdx = obsIdx + 1;
rot = pi/2;     trans = [6.0, 6.5];     scale = 0.8;    rectWidth = 1.3;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 18
obsIdx = obsIdx + 1;
rot = 0;     trans = [9.5, 9];     scale = 1.;    rectWidth = 1.8;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 19
obsIdx = obsIdx + 1;
rot = pi/2;     trans = [18.5, 2.5];     scale = 0.7;       rectWidth = 5;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 20
obsIdx = obsIdx + 1;
rot = 0;     trans = [8.2, 6.0];     scale = 0.7;       rectWidth = 3;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 21
obsIdx = obsIdx + 1;
rot = pi/4;     trans = [7.3, 2.5];     scale = 2.5;    rectWidth = 1.0;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 22
obsIdx = obsIdx + 1;
rot = -pi/4;     trans = [13, 6];     scale = 0.6;       rectWidth = 5;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 23
obsIdx = obsIdx + 1;
rot = pi/4;     trans = [12, 8.5];     scale = 0.6;       rectWidth = 4;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 24
obsIdx = obsIdx + 1;
rot = 0;     trans = [14.7, 8.5];     scale = 0.5;       rectWidth = 6;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% 25
obsIdx = obsIdx + 1;
rot = 0;     trans = [18.8, 7.7];     scale = 0.7;
shape = Hexagon;
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

% Plot obstacles
figure;
for k = 1 : length(vectorObs)
    plot(vectorObs{k}); hold on;
    
    if PLOT_OBSTACLE_ID
        if k > 1
            % Computes the centroid to plot the obstacle numbers
            [xc, yc] = centroid(vectorObs{k});
            text(xc, yc, sprintf('%d', k - 1), "HorizontalAlignment", "center","FontSize",6)
        end
    end
end
pbaspect([1 1 1]);
% axis equal;

if SAVE_OBSTACLES
    setName = sprintf('Obstacles_%dSet%d', length(vectorObs) - 1, setID);
    save(setName,'vectorObs');
end