%% Create sets of enviroments
close all
clear
clc

SAVE_OBSTACLES = true;
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
load boundaryBenchmark1.mat

setID = setID + 1;
obsIdx = 1;

obsIdx = obsIdx + 1;
rot = pi/2;     trans = [2.5, 2.0];     scale = 0.25;       rectWidth = 12;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = 0;     trans = [4.2, 4.2];     scale = 0.7;       rectWidth = 1;
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
pbaspect([1 1 1]); axis equal;

if SAVE_OBSTACLES
    setName = sprintf('Obstacles_Benchmark%d', setID);
    %setName = sprintf('Obstacles_Set%d', setID);
    save(setName,'vectorObs');
end

clearvars vectorObs 

%% Set 2
load boundaryBenchmark2.mat

setID = setID + 1;
obsIdx = 1;

obsIdx = obsIdx + 1;
rot = pi/2;     trans = [7.5, 2.0];     scale = 0.25;       rectWidth = 12;
shape = Rectangle(rectWidth);
obs =  T(trans) * S(scale) * R(rot) * T(-mean(shape,2)) * shape;
vectorObs{obsIdx} = polyshape(obs(1:2, :)');

obsIdx = obsIdx + 1;
rot = 0;     trans = [14.2, 4.2];     scale = 0.7;       rectWidth = 1;
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
            text(xc, yc, sprintf('%d', k - 1), "HorizontalAlignment", "center")
        end
    end
end
pbaspect([1 1 1]); axis equal;

if SAVE_OBSTACLES
    setName = sprintf('Obstacles_Benchmark%d', setID);
    %setName = sprintf('Obstacles_Set%d', setID);
    save(setName,'vectorObs');
end
