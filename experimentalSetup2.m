%% Create an environment from the setup

% clc
% clear
% close all

% SETUP_1 = 0;  
% SETUP_2 = 1;   
% % Define setup environment: {SETUP_1 | SETUP_2}
% ENVIRONMENT_SETUP = SETUP_1;

SAVE_EXPERIMENTAL_OBSTACLES =  true;

% 1240 x 2035mm (H x W): Dimension of the physical setup
% 3044 x 1848px (H x W): Dimension of the image
% Rescale the image to be in a 1-to-1 ratio
Hbox = 1240;
Wbox = 2035;


% Get the relative folder path to access the sub-folder 'figs'
cd(strcat('figs',filesep)); fp_figs = [pwd filesep]; cd ..;

vertex = 1.0e+03 *[
    0.1067    0.1482
    0.2403    0.9151
    0.3705    0.2819
    0.5500    0.6600 % 0.5370    0.6559
    0.7364    1.0453
    0.8314    0.5317
    0.9475    0.3417
    1.1269    0.8870
    1.2993    0.2819
    1.4963    0.5397
    1.6827    0.8201
    1.7566    0.4050
    1.8973    0.1728
    1.8938    1.0558
];

figure;
if (EXPERIMENT_SETUP == SETUP_1)
    img = imread(strcat(fp_figs,"setup.jpg"));
    
    img_size = size(img);
    ratioW = Wbox/img_size(2);
    ratioH = Hbox/img_size(1);
    ratio = mean([ratioH, ratioW]);
    img = imresize(img, ratio);

    imshow(img), hold on;
    
    % Boundaries
    env_height_in = size(img,1); 
    env_width_in = size(img,2); 
    wall_thickness = 0.02*(env_height_in + env_width_in)/2;  % Thickness
    env_width_out = env_height_in + wall_thickness;
    env_length_out = env_width_in + wall_thickness;
    
    x_inPlot = [0 env_height_in env_height_in 0];
    y_inPlot = [0 0 env_width_in env_width_in];
    x_outPlot = [-wall_thickness env_width_out env_width_out -wall_thickness];
    y_outPlot = [-wall_thickness -wall_thickness env_length_out env_length_out];
    
    x_in = x_inPlot/1e3;
    y_in = y_inPlot/1e3;
    x_out = x_outPlot/1e3;
    y_out = y_outPlot/1e3;
    
    env_boundariesPlot = polyshape({y_inPlot,y_outPlot}, {x_inPlot,x_outPlot});
    plot(env_boundariesPlot); axis tight;

    % Takes the input data and draws the limit of the environment    
    env_boundaries = polyshape({y_in,y_out}, {x_in,x_out});
    ind = 1;
    vectorObs{ind} = env_boundaries;
    ind = ind + 1;
    
    plot(env_boundaries); hold on;
    % pbaspect([1 1 1]);
    
    
    
    % Obstacle 1
    xyObsPlot = [vertex(1:3,:)];
    plot(polyshape(xyObsPlot(:,1), xyObsPlot(:,2)), FaceColor=[0.5 0.5 0.5], FaceAlpha=mu(ind));
    % Computes the centroid to plot the obstacle numbers
    [xc, yc] = centroid(polyshape(xyObsPlot(:,1), xyObsPlot(:,2))); text(xc, yc, sprintf('%d', ind-1), "HorizontalAlignment","center")
    xyObsPlot(:,2) = size(img,1) - xyObsPlot(:,2); % Take y's complement
    xyObs = xyObsPlot/1e3;
    obs = polyshape(xyObs(:,1), xyObs(:,2));
    vectorObs{ind} = obs;
    ind = ind + 1;
    
    % Obstacle 2
    xyObsPlot = [vertex(9,:); vertex(12:13,:)];
    plot(polyshape(xyObsPlot(:,1), xyObsPlot(:,2)), FaceColor=[0.5 0.5 0.5], FaceAlpha=mu(ind));
    % Computes the centroid to plot the obstacle numbers
    [xc, yc] = centroid(polyshape(xyObsPlot(:,1), xyObsPlot(:,2))); text(xc, yc, sprintf('%d', ind-1), "HorizontalAlignment","center")
    xyObsPlot(:,2) = size(img,1) - xyObsPlot(:,2); % Take y's complement
    xyObs = xyObsPlot/1e3;
    obs = polyshape(xyObs(:,1), xyObs(:,2));
    vectorObs{ind} = obs;
    ind = ind + 1;
    
    % Obstacle 3
    xyObsPlot = [vertex(6:7,:); vertex(10,:)];
    plot(polyshape(xyObsPlot(:,1), xyObsPlot(:,2)), FaceColor=[0.5 0.5 0.5], FaceAlpha=mu(ind));
    % Computes the centroid to plot the obstacle numbers
    [xc, yc] = centroid(polyshape(xyObsPlot(:,1), xyObsPlot(:,2))); text(xc, yc, sprintf('%d', ind-1), "HorizontalAlignment","center")
    xyObsPlot(:,2) = size(img,1) - xyObsPlot(:,2); % Take y's complement
    xyObs = xyObsPlot/1e3;
    obs = polyshape(xyObs(:,1), xyObs(:,2));
    vectorObs{ind} = obs;
    ind = ind + 1;
    
    % Obstacle 4
    xyObsPlot = [vertex(5,:); vertex(8,:); vertex(11,:); vertex(14,:)];
    plot(polyshape(xyObsPlot(:,1), xyObsPlot(:,2)), FaceColor=[0.5 0.5 0.5], FaceAlpha=mu(ind));
    % Computes the centroid to plot the obstacle numbers
    [xc, yc] = centroid(polyshape(xyObsPlot(:,1), xyObsPlot(:,2))); text(xc, yc, sprintf('%d', ind-1), "HorizontalAlignment","center")
    xyObsPlot(:,2) = size(img,1) - xyObsPlot(:,2); % Take y's complement
    xyObs = xyObsPlot/1e3;
    obs = polyshape(xyObs(:,1), xyObs(:,2));
    vectorObs{ind} = obs;
    ind = ind + 1;
    
    if(SAVE_EXPERIMENTAL_OBSTACLES)
        save('Obstacles_4setup1','vectorObs')
    end

elseif(EXPERIMENT_SETUP == SETUP_2)
    img = imread(strcat(fp_figs,"setup.jpg"));

    img_size = size(img);
    ratioW = Wbox/img_size(2);
    ratioH = Hbox/img_size(1);
    ratio = mean([ratioH, ratioW]);
    img = imresize(img, ratio);

    imshow(img), hold on;
    
    % Boundaries
    env_height_in = size(img,1); 
    env_width_in = size(img,2); 
    wall_thickness = 0.02*(env_height_in + env_width_in)/2;  % Thickness
    env_width_out = env_height_in + wall_thickness;
    env_length_out = env_width_in + wall_thickness;
    
    x_inPlot = [0 env_height_in env_height_in 0];
    y_inPlot = [0 0 env_width_in env_width_in];
    x_outPlot = [-wall_thickness env_width_out env_width_out -wall_thickness];
    y_outPlot = [-wall_thickness -wall_thickness env_length_out env_length_out];
    
    x_in = x_inPlot/1e3;
    y_in = y_inPlot/1e3;
    x_out = x_outPlot/1e3;
    y_out = y_outPlot/1e3;
    
    env_boundariesPlot = polyshape({y_inPlot,y_outPlot}, {x_inPlot,x_outPlot});
    plot(env_boundariesPlot); axis tight;

    % Takes the input data and draws the limit of the environment    
    env_boundaries = polyshape({y_in,y_out}, {x_in,x_out});
    ind = 1;
    vectorObs{ind} = env_boundaries;
    ind = ind + 1;
    
    plot(env_boundaries); hold on;    
    
    % Obstacle 1
    xyObsPlot = [vertex(1:2,:); vertex(3,:)];
    plot(polyshape(xyObsPlot(:,1), xyObsPlot(:,2)), FaceColor=[0.5 0.5 0.5], FaceAlpha=mu(ind));
    % Computes the centroid to plot the obstacle numbers
    [xc, yc] = centroid(polyshape(xyObsPlot(:,1), xyObsPlot(:,2))); text(xc, yc, sprintf('%d', ind-1), "HorizontalAlignment","center")
    xyObsPlot(:,2) = size(img,1) - xyObsPlot(:,2); % Take y's complement
    xyObs = xyObsPlot/1e3;
    obs = polyshape(xyObs(:,1), xyObs(:,2));
    vectorObs{ind} = obs;
    ind = ind + 1;
    
    % Obstacle 2
    xyObsPlot = [vertex(4,:); vertex(5,:); vertex(8,:)];
    plot(polyshape(xyObsPlot(:,1), xyObsPlot(:,2)), FaceColor=[0.5 0.5 0.5], FaceAlpha=mu(ind));
    % Computes the centroid to plot the obstacle numbers
    [xc, yc] = centroid(polyshape(xyObsPlot(:,1), xyObsPlot(:,2))); text(xc, yc, sprintf('%d', ind-1), "HorizontalAlignment","center")
    xyObsPlot(:,2) = size(img,1) - xyObsPlot(:,2); % Take y's complement
    xyObs = xyObsPlot/1e3;
    obs = polyshape(xyObs(:,1), xyObs(:,2));
    vectorObs{ind} = obs;
    ind = ind + 1;
    
    % Obstacle 3
    xyObsPlot = [vertex(6:7,:); vertex(10,:)];
    plot(polyshape(xyObsPlot(:,1), xyObsPlot(:,2)), FaceColor=[0.5 0.5 0.5], FaceAlpha=mu(ind));
    % Computes the centroid to plot the obstacle numbers
    [xc, yc] = centroid(polyshape(xyObsPlot(:,1), xyObsPlot(:,2))); text(xc, yc, sprintf('%d', ind-1), "HorizontalAlignment","center")
    xyObsPlot(:,2) = size(img,1) - xyObsPlot(:,2); % Take y's complement
    xyObs = xyObsPlot/1e3;
    obs = polyshape(xyObs(:,1), xyObs(:,2));
    vectorObs{ind} = obs;
    ind = ind + 1;
    
    % Obstacle 4
    xyObsPlot = [vertex(11,:); vertex(12,:); vertex(9,:); vertex(13,:); vertex(14,:)];
    plot(polyshape(xyObsPlot(:,1), xyObsPlot(:,2)), FaceColor=[0.5 0.5 0.5], FaceAlpha=mu(ind));
    % Computes the centroid to plot the obstacle numbers
    [xc, yc] = centroid(polyshape(xyObsPlot(:,1), xyObsPlot(:,2))); text(xc, yc, sprintf('%d', ind-1), "HorizontalAlignment","center")
    xyObsPlot(:,2) = size(img,1) - xyObsPlot(:,2); % Take y's complement
    xyObs = xyObsPlot/1e3;
    obs = polyshape(xyObs(:,1), xyObs(:,2));
    vectorObs{ind} = obs;
    ind = ind + 1;
    
    if(SAVE_EXPERIMENTAL_OBSTACLES)
        save('Obstacles_4setup2','vectorObs')
    end
end

% Dijkstra on G
pathG_Dks_Map = [];
for k=1 : length(pathG_Dks)
    xy = G.coord(pathG_Dks(k));
    pathG_Dks_Map = [pathG_Dks_Map; 1e3*xy(1), Hbox-1e3*xy(2)];
end
plot(pathG_Dks_Map(:,1), pathG_Dks_Map(:,2), MarkerFaceColor="yellow", MarkerEdgeColor="black", Marker="o", MarkerSize=8, LineWidth=2, Color="red")

% BFS on GE
pathGE_BFS_Map = [];
for k=1 : length(pathG_BFS)
    xy = G.coord(pathG_BFS(k));
    pathGE_BFS_Map = [pathGE_BFS_Map; 1e3*xy(1), Hbox-1e3*xy(2)];
end
plot(pathGE_BFS_Map(:,1), pathGE_BFS_Map(:,2), MarkerFaceColor="yellow", MarkerEdgeColor="black", Marker="o", MarkerSize=8, LineWidth=2, Color="blue")

% Dijkstra on GE (Tension-aware)
pathGE_Dks_Map = [];
for k=1 : length(pathGE_Dks)
    xy = G.coord(pathGE_Dks(k));
    pathGE_Dks_Map = [pathGE_Dks_Map; 1e3*xy(1), Hbox-1e3*xy(2)];
end
plot(pathGE_Dks_Map(:,1), pathGE_Dks_Map(:,2), MarkerFaceColor="yellow", MarkerEdgeColor="black", Marker="o", MarkerSize=8, LineWidth=2, Color="green")
plot(pathGE_Dks_Map(1,1), pathGE_Dks_Map(1,2), MarkerFaceColor="red", MarkerEdgeColor="black", Marker="o", MarkerSize=8)
plot(pathGE_Dks_Map(end,1), pathGE_Dks_Map(end,2), MarkerFaceColor="green", MarkerEdgeColor="black", Marker="o", MarkerSize=8)

if SAVE_FILE  
    filename = strcat(sprintf('result_paths_img_setup%d',EXPERIMENT_SETUP),'.eps');
    filepath = strcat(fp_figs, filename);
    saveas(gcf,filepath,'epsc');
end