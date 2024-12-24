%% MINIMUM TENSION MOTION PLANNER
% AUTHOR:           Rogerio Lima
% DATE & PLACE:     November/2023 @ Morgantown/USA
% VERSION:          1
% DESCRIPTION:      This version searches the tension-aware graph (Gt) with 
%                   Dijkstra, and searches G using Dijkstra and BFS.
% =========================================================================

%% Pre set-up
clc
clear
close all


ENABLE_RANDOM_FRICTION = false;
VERTEX = 0;     % Each vertex may have different coeff. of friction
OBSTACLE = 1;   % Vertices of a single obstacle have the same friction
ENVIRONMENT = 2;% All vertices (all obstacles) have the same friction

% Define friction mode: {VERTEX | OBSTACLE | ENVIRONMENT}
FRICTION_MODE = OBSTACLE;

% This affects the script "experimentalSetup.m", which plots the image setup 
% along with the paths. It assumes the path is run first in this script.
SETUP_0 = 0; % Simulation only (SIM_n, n = 1...6)
SETUP_1 = 1; % Experiment 1
SETUP_2 = 2; % Experiment 2

% Define setup environment: {SETUP_0 |SETUP_1 | SETUP_2}
EXPERIMENT_SETUP = SETUP_0;

% Simulation environment (Valid only with SETUP_0)
SIM_0 = 0;
SIM_1 = 1;
SIM_2 = 2;
SIM_3 = 3;
SIM_4 = 4;
SIM_5 = 5;
SIM_6 = 6; % Choose environment and start/goal position manually

% Only takes effect if "EXPERIMENT_SETUP = SETUP_0"
SIMULATION_SETUP = SIM_1;

% Print experimental results
PRINT_EXPERIMENTAL_RESULTS = true;

% Start/Goal positions: false:hard-coded | true:UI
UI_INTERFACE = false;

% Save file: 0:disable | 1:enable
SAVE_FILE = false;

% Debug mode: 0:disable | 1:text | 2:plot
DEBUG = 0;

% Pause mode: false:disable | true:enable
PAUSE_ON = false;

% Save prints to a "log.txt" file
DATALOGGER_ON_DEBUG = true;

if (DATALOGGER_ON_DEBUG && DEBUG > 0)
    fid = fopen('log.txt', 'w');
end

% Get the relative folder path to access the sub-folder 'figs'
cd(strcat('figs',filesep)); fp_figs = [pwd filesep]; cd ..;
% opengl('save','hardware') 
fs = 16;                % font size
afs = 14;               % axis font size
lw = 1.0;               % line width
dt = 0.8;               % offset distance for the numbers (vertices ids)

% Warning off
warning('off', 'MATLAB:polyshape:boundary3Points');

% Global variables
global vectorObs;
lastNode = 0;

% Local Variables
rho = 0.01; % Linear mass density of the tether, [kg/m]
g = 9.81;   % Gravitational acceleration, [m/s^2]
mu_floor = 0.6;         % Floor's coefficient of friction, []
C = [0 mu_floor 1];     % Term that adapts the linear /friction term to the cost function. {C = mu_g 2D planner; C = 1 for 3D planner)
T0 = 1;
                        % 1: Capstan mode only
COST_FUNCTION_MODE = 2; % 2: Capstan + linear friction for ground vehicles
                        % 3: Capstan + linear friction for aerial vehicles
%% Real-world experiment
if (EXPERIMENT_SETUP ~= 0)
    results = printExperimentalResults();
end

%% Simulation
% Create a graph
G = PGraph();
Gt = PGraph();
rng(5)

% Load obstacles
if (EXPERIMENT_SETUP == SETUP_1)
    obsFilename = 'Obstacles_4setup1';
    ps = [0.04, 1.08];  % Start position:   Obstacles_4setup1*
    pg = [1.93, 0.23];  % Goal position:    Obstacles_4setup1* 1.95, 0.25
elseif (EXPERIMENT_SETUP == SETUP_2)
    obsFilename = 'Obstacles_4setup2';
    ps = [0.05, 0.17];  % Start position:  Obstacles_4setup2* [0.07, 0.07];
    pg = [1.98, 1.10];  % Goal position:   Obstacles_4setup2* [1.98, 1.10]
else
    if (SIMULATION_SETUP == SIM_0)
        obsFilename = 'Obstacles_3';
        mu = [1, 0.9, 0.2, 0.9];
        ps = [0.5, 4.63];
        pg = [4.65, 0.30];
    elseif (SIMULATION_SETUP == SIM_1)
        obsFilename = 'Obstacles_3Set1';
        mu = [1, 0.4, 0.2, 0.4];
        ps = [1.50, 5.00];
        pg = [8.50, 5.00];
    elseif (SIMULATION_SETUP == SIM_2)
        obsFilename = 'Obstacles_9Set2';
        mu =  [1, 0.4, 0.2, 0.4, 0.2, 0.4, 0.4, 0.1, 0.4, 0.1];
        ps = [1.00, 9.00];
        pg = [9.00, 1.50];
    elseif (SIMULATION_SETUP == SIM_3)
        obsFilename = 'Obstacles_Benchmark1';
        mu = [1, 0.5, 0.1];
        ps = [0.20, 4.80];
        pg = [4.80, 0.20];
    elseif (SIMULATION_SETUP == SIM_4)
        obsFilename = 'Obstacles_Benchmark2';
        mu = [1, 0.5, 0.1];
        ps = [0.20, 4.80];
        pg = [14.80, 0.20];
    elseif (SIMULATION_SETUP == SIM_5)      % 1 obstacle
        obsFilename = 'Obstacles_1Set4';
        ps = [1.50, 5.00];
        pg = [8.50, 5.00];
        mu = rand(1, 5);    % Heterogeneous obstacles
    elseif (SIMULATION_SETUP == SIM_6)      % Choose environments
        obsFilename = 'Obstacles_5';
        mu = rand(1, length(vectorObs));    
        ps = [0.4758, 4.6020];              % Choose start position
        pg = [4.7177, 0.2755];              % Choose goal position
    end
    

end
load(obsFilename)


% Array of friction coefficients generated randomly
if (ENABLE_RANDOM_FRICTION)
    if (FRICTION_MODE == OBSTACLE)
        mu = rand(1, length(vectorObs));    % Heterogeneous obstacles
        % For testing special cases depending on the # of obstacles
        mu2 = [1      0.5   0.1];
        mu3 = [1      0.05  0.9     0.1];
        mu4 = [1      0.6   0.6     0.4       0.6];
        mu5 = [1      0.1   0.9     0.1       0.1     0.1];
        mu7 = [1      0.1   0.9     0.1       0.9     0.5     0.9     0.9];
        mu9=  [1      0.1   0.1     0.1       0.5     0.9     0.9     0.1       0.9     0.1];
        mu10= [1      0.1   0.8     0.1       0.5     0.9     0.9     0.1       0.8     0.1     0.1];
    %     mu = mu2
        
    elseif (FRICTION_MODE == ENVIRONMENT)
        mu = ones(1, length(vectorObs));   % Homogeneous obstacles
    end
else
    % Experiment setup 
    if (EXPERIMENT_SETUP ~= 0)
        muL = results.mu.plastic_cylinder;  % Lower mu
        muH = results.mu.cardboard;         % Higher mu
        mu = [1 muH muH muL muH]; % Original
    end
end


% Plot the environment boundaries
plot(vectorObs{1}); hold on;

% Skip 1st obstacle (environment boundaries)
for m = 2 : length(vectorObs) 
    aux = vectorObs{m}.Vertices; 
    
    if (FRICTION_MODE == VERTEX)
        plot(vectorObs{m}, FaceColor=[0.5 0.5 0.5], FaceAlpha=0.7); hold on;
    else
        plot(vectorObs{m}, FaceColor=[0.5 0.5 0.5], FaceAlpha=mu(m)); hold on;
    end

    % Computes the centroid to plot the obstacle numbers
    [xc, yc] = centroid(vectorObs{m});
    text(xc, yc, sprintf('%d', m-1), "HorizontalAlignment","center")
end
pbaspect([1 1 1]);

if UI_INTERFACE
    % % Warnings the user to pick up the starting point
    f = msgbox('Please, select the Start Position...', 'Start','warn');
    uiwait(f);
    % Start position
    [x_start,y_start] = ginput(1);
    ps = [x_start, y_start];

    % Warnings the user to pick up the goal position
    f = msgbox('Please, select the Goal Position...', 'Goal','warn');
    uiwait(f);
    % Goal position
    [x_goal,y_goal] = ginput(1);
    pg = [x_goal, y_goal];
end    


sorted_nodes = ps;
idx = 2;
% Skip 1st obstacle (environment boundaries)
for m=2:length(vectorObs) 
    aux = vectorObs{m}.Vertices;    
    for n=1:size(aux,1)
        sorted_nodes(idx,:)=aux(n,:);
        idx = idx + 1;
    end
end
% Goal position is the last position of all nodes.
sorted_nodes(idx, :) = pg;

% Define friction at vertex level
if (FRICTION_MODE == VERTEX)
    mu = rand(1, length(sorted_nodes));
end

start = 1;
goal = size(sorted_nodes, 1);
% Virtual nodes used in the graph where the nodes are formed by edges. A
% pair of connecting vertices in graph G is a single node in graph Gt.
virtual_start = start - 1;
virtual_goal = goal + 1; 

% DEBUG
if (DEBUG == 2)
    for n=1:size(sorted_nodes, 1)
        plot(sorted_nodes(n,1), sorted_nodes(n,2), 'o'); hold on
        if PAUSE_ON
            pause(3);
        end
    end
end

if (EXPERIMENT_SETUP ~= 0)
    PARTITION_LENGTH = 0.002;
else
    PARTITION_LENGTH = 0.01;
end
connectIdx = 1;
num_nodes = size(sorted_nodes, 1);
visibilityTable = false(num_nodes);            % Visibility table
for k=1:num_nodes - 1
    for j=k+1:num_nodes
        pp = partitionLine(sorted_nodes(k,:),...
            sorted_nodes(j,:), PARTITION_LENGTH);
        flagFreeNode = true;
        for m=1:size(pp,1)
            if (~isFree(pp(m,:)))
                flagFreeNode = false;   % Needed just one false
            end
        end

        if (flagFreeNode)
            visibilityTable(k,j) = true;
            visibilityTable(j,k) = true;
            kj(connectIdx,1:2) = [k,j]; % Pair of valid (connecting) nodes
            jk(connectIdx,1:2) = [j,k]; % Bi-direction visibility
            connectIdx = connectIdx + 1;
        end
    end
end

visIdx = [kj; jk];

% Virtual start node: [0,1] (Should come before all nodes)
firstNodeEdgeID = 1;
Gt.add_node([virtual_start start],'name', firstNodeEdgeID);

% Dimensions of visibility table
sizeVisTable = size(visibilityTable,1);

% DEBUG LEVEL 1
if (DEBUG >= 1)
    % Prints: nodeEdgeIndex: [vertex_1, vertex2] --> nodeID
    fprintf('%d: [%d,%d] --> %d\r\n',1, virtual_start, 1, matrix2digit(1,1, sizeVisTable ));
    if (DATALOGGER_ON_DEBUG)
        fprintf(fid,'%d: [%d,%d] --> %d\r\n',1, virtual_start, 1, matrix2digit(1,1, sizeVisTable ));
    end
end
% DEBUG LEVEL 2
if (DEBUG == 2)
    figure;
    Gt.plot();
    lmt = 0.3;
    axis equal; xlim([-lmt virtual_goal+lmt]); ylim([-lmt virtual_goal+lmt]);
end


% Encodes the visibility table in a coordinate pair where each coordinate
% represents the vertex ID. For example, if the second and fourth vertices
% are visible to each other, then there will be a pair (2,4) that indicates
% its visibility.
visiblePairs = edge2node(visibilityTable);

vertexIdxID = [1, 1];
numNodesEdge = length(visiblePairs);
for k = 2 : numNodesEdge+1
    % Add nodes to the graph
    Gt.add_node(visiblePairs{k-1}, 'name', matrix2digit(visiblePairs{k-1}(1), visiblePairs{k-1}(2), sizeVisTable));
    d = matrix2digit(visiblePairs{k-1}(1),visiblePairs{k-1}(2), sizeVisTable);
    vertexIdxID(k,:) = [k,d];
    
    % DEBUG LEVEL 1
    if (DEBUG >= 1)
        % Prints: nodeEdgeIndex: [vertex_1, vertex2] --> nodeID
        fprintf('%d: [%d,%d] --> %d\r\n', k, visiblePairs{k-1}(1), visiblePairs{k-1}(2), matrix2digit(visiblePairs{k-1}(1), visiblePairs{k-1}(2), sizeVisTable));
        fprintf('  [nodeIdx, nodeID] = [%d,%d] \r\n\n', k, d);
        if (DATALOGGER_ON_DEBUG)
            fprintf(fid, '%d: [%d,%d] --> %d\r\n', k, visiblePairs{k-1}(1), visiblePairs{k-1}(2), matrix2digit(visiblePairs{k-1}(1), visiblePairs{k-1}(2), sizeVisTable));
            fprintf(fid, '  [nodeIdx, nodeID] = [%d,%d] \r\n\n', k, d);
        end
    end
    % DEBUG LEVEL 2
    if (DEBUG == 2)
        Gt.plot()
    end

end

% Virtual goal node: [goal, goal+1]
% lastNodeEdgeID is encoded to be after the last possible node.
lastNodeEdgeID = sizeVisTable^2 + 1;
Gt.add_node([goal virtual_goal], 'name', lastNodeEdgeID);
lastNodeEdgeIndex = Gt.n;

% DEBUG LEVEL 1
if (DEBUG == 1)
    % Prints: nodeEdgeIndex: [vertex_1, vertex2] --> nodeID
    fprintf('%d: [%d,%d] --> %d\r\n',k+1, goal,virtual_goal, lastNodeEdgeID);
end

all_nodes = [visIdx(:,1);visIdx(:,2)];  % All connected nodes in an array
[b,m1,n1] = unique(all_nodes,'first');  % Disregard repeated nodes
[c1,d1] =sort(m1);                      % and take just one
b = b(d1);

% Add nodes to the graph G
for k=1:size(sorted_nodes,1)
    G.add_node(sorted_nodes(k,:));
end

% Create a Boolean matrix initialized with all entries as 'false' to allow
% creating unique edges in graph G
visibilityCheck = false(size(visibilityTable));

% Containers map for tracking the edge IDs
edgeID = containers.Map;    

% Add first and last nodes to the list
kkext = [virtual_start, start 
         visIdx
         goal, virtual_goal];

% Iterates over all pair of connecting vertices
iter = 0;
numConnectingVertices = size(kkext,1);
for cvItr = 1 : numConnectingVertices
    iter = iter + 1;
    nextPairs = find(kkext(:,1) == kkext(cvItr,2));
    for npItr = nextPairs'

        pt1 = kkext(cvItr, 1);
        pt2 = kkext(cvItr, 2);
        pt3 = kkext(npItr, 2);

        % Case for the virtual start and goal points
        if (pt1 == virtual_start || pt3 == virtual_goal)
            % Add two nodes and an egde between them for each triad of valid vertices.
            nodeIdx1 = Gt.lookup(sprintf('%d', matrix2digit(pt1, pt2, sizeVisTable)));
            nodeIdx2 = Gt.lookup(sprintf('%d', matrix2digit(pt2, pt3, sizeVisTable)));

            % --- Changing in the cost function
            if (pt1 == virtual_start)
                vtxb = G.coord(pt2);
                vtxc = G.coord(pt3);
                side_a = hypot(vtxb(1) - vtxc(1), vtxb(2) - vtxc(2));
                Le = side_a;            % side_a = BC,
            end
            if (pt3 == virtual_goal)
                vtxa = G.coord(pt1);
                vtxb = G.coord(pt2); 
                side_c = hypot(vtxa(1) - vtxb(1), vtxa(2) - vtxb(2));
                Le = side_c;            % side_c = AB
            end
            linearFriction = (rho * g * Le/2);
            % --- Changing in the cost function

            Gt.add_edge(nodeIdx1, nodeIdx2);
            costFunction = C(COST_FUNCTION_MODE) * linearFriction;        
            Gt.setcost(Gt.ne, costFunction);

            key = sprintf("%d.%d", nodeIdx1, nodeIdx2);

            % DEBUG LEVEL 1
            if (DEBUG >= 1)
                % Prints: nodeEdgeIndex: [vertex_1, vertex2] --> nodeID
                fprintf('%d: [%d,%d] --> %d\r\n', iter, pt1, pt2, matrix2digit(pt1, pt2, sizeVisTable ));
                fprintf('%d: [%d,%d] --> %d\r\n', iter, pt2, pt3, matrix2digit(pt2, pt3, sizeVisTable ));
                fprintf('%d: cost = %1.2f \r\n\n', iter, costFunction); 
                if (DATALOGGER_ON_DEBUG)
                    fprintf(fid, '%d: [%d,%d] --> %d\r\n', iter, pt1, pt2, matrix2digit(pt1, pt2, sizeVisTable ));
                    fprintf(fid, '%d: [%d,%d] --> %d\r\n', iter, pt2, pt3, matrix2digit(pt2, pt3, sizeVisTable ));
                    fprintf(fid, '%d: cost = %1.2f \r\n\n', iter, costFunction); 
                end
            end
        % General case for all other nodes
        else
            % Triangle sides
            vtxa = G.coord(pt1);
            vtxb = G.coord(pt2);
            vtxc = G.coord(pt3);
            side_a = hypot(vtxb(1) - vtxc(1), vtxb(2) - vtxc(2));
            side_b = hypot(vtxc(1) - vtxa(1), vtxc(2) - vtxa(2));
            side_c = hypot(vtxa(1) - vtxb(1), vtxa(2) - vtxb(2));
            angB = computeAngleB(side_a, side_b, side_c);
    
            [obs, obsIdx] = pointBelongsToObstacle(sorted_nodes(pt2,:), vectorObs);
    
            % Add two nodes and an egde between them for each triad of valid vertices.
            nodeIdx1 = Gt.lookup(sprintf('%d', matrix2digit(pt1, pt2, sizeVisTable)));
            nodeIdx2 = Gt.lookup(sprintf('%d', matrix2digit(pt2, pt3, sizeVisTable)));

            key = sprintf("%d.%d", nodeIdx1, nodeIdx2);
            %Gt.add_edge(nodeIdx1, nodeIdx2);

            isTriangleValid = checkTriangulation(sorted_nodes(pt1,:), sorted_nodes(pt2,:), sorted_nodes(pt3,:), obs{1});
    
            if (isTriangleValid)
                theta = pi - angB;
                Le = side_a + side_c; % side_a = BC, side_c = AB
                linearFriction = (rho * g * Le/2);
                % Capstan equation as cost function
                if ((FRICTION_MODE == OBSTACLE) || (FRICTION_MODE == ENVIRONMENT))
                    costFunction = T0 * (exp(mu(obsIdx) * theta) - 1) + C(COST_FUNCTION_MODE) * linearFriction;
                else
                    costFunction = T0 * (exp(mu(pt2) * theta) - 1) + C(COST_FUNCTION_MODE) * linearFriction;
                end

                % Avoids edge multiplicity in Gt (creates only one edge)
                %if (~edgeID.isKey(key))
                Gt.add_edge(nodeIdx1, nodeIdx2);
                Gt.setcost(Gt.ne, costFunction);
                %end

                % Avoids edge multiplicity (creates only one edge)
                if(~visibilityCheck(pt1, pt2))
                    G.add_edge(pt1, pt2);
                    cf = G.distance(pt1, pt2);
                    G.setcost(G.ne, cf);

                    visibilityCheck(pt1, pt2) = true;
                end
                
                % Avoids edge multiplicity (creates only one edge)
                if(~visibilityCheck(pt2, pt3))
                    G.add_edge(pt2, pt3);
                    cf = G.distance(pt2, pt3);
                    G.setcost(G.ne, cf);

                    visibilityCheck(pt2, pt3) = true;
                end
            else
                costFunction = 1e9;
            end
            key = sprintf("%d.%d", nodeIdx1, nodeIdx2);

            % DEBUG LEVEL 1
            if (DEBUG >= 1)
                % Prints: nodeEdgeIndex: [vertex_1, vertex2] --> nodeID
                fprintf('%d: [%d,%d] --> %d\r\n', iter, pt1, pt2, matrix2digit(pt1, pt2, sizeVisTable ));
                fprintf('%d: [%d,%d] --> %d\r\n', iter, pt2, pt3, matrix2digit(pt2, pt3, sizeVisTable ));
                fprintf('%d: cost = %1.2f \r\n\n', iter, costFunction);
                if (DATALOGGER_ON_DEBUG)
                    fprintf(fid, '%d: [%d,%d] --> %d\r\n', iter, pt1, pt2, matrix2digit(pt1, pt2, sizeVisTable ));
                    fprintf(fid, '%d: [%d,%d] --> %d\r\n', iter, pt2, pt3, matrix2digit(pt2, pt3, sizeVisTable ));
                    fprintf(fid, '%d: cost = %1.2f \r\n\n', iter, costFunction);
                end
            end
        end
        % Records the edgeID in function of pair (nodeID1,nodeID2)
        edgeID(key) = Gt.ne;
        
        % DEBUG LEVEL 1
        if (DEBUG >= 1)
            %{
            nodeID1 = Gt.lookup(sprintf('%d', matrix2digit(pt1, pt2, sizeVisTable)))
            nodeID2 = Gt.lookup(sprintf('%d', matrix2digit(pt2, pt3, sizeVisTable)))
            retrievedEdgeIDFromVertices = edgeID( sprintf("%d.%d", nodeID1, nodeID2))
            Gt.cost(retrievedEdgeIDFromVertices)
            %}
        end
    end
end

% Graph search on Visibility Graph G: BFS
pathG_BFS = bfs2(G, G.n);
tetherLengthGE_BFS = calcTetherLength(G, pathG_BFS);

% Graph search on Visibility Graph G: Dijkstra (has physical meaning)
[pathG_Dks, costDks] = dijkstra(G, start, G.n);
tetherLengthG_Dks = calcTetherLength(G, pathG_Dks);

% Graph search on Visibility Graph Gt: Tension-aware graph (Dijkstra on Gt)
[pathGE_Dks_, costGE_Dks] = dijkstra(Gt, start, Gt.n);
pathGE_Dks = [];
for k = 2 : length (pathGE_Dks_) - 1
    pathGE_Dks = [pathGE_Dks, Gt.coord(pathGE_Dks_(k))'];
end
pathGE_Dks = unique(pathGE_Dks, 'stable'); % Remove repeated entries
tetherLengthGE_Dks = calcTetherLength(G, pathGE_Dks);


% Uncomment to plot the visibility graph
% G.plot();
G.highlight_path(pathG_Dks, 'EdgeColor', 'r', 'NodeSize', 7);   % Dijkstra (A* with h(n) = 0)
G.highlight_path(pathGE_Dks, 'EdgeColor', 'g', 'NodeSize', 7);  % Dijkstra - Node-edge domain
G.highlight_path(pathG_BFS, 'EdgeColor', 'b', 'NodeSize', 7); % Path from A* motion planning
plot(ps(1), ps(2), 'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k' ,'MarkerSize', 7);
plot(pg(1), pg(2), 'Marker','o','MarkerFaceColor','g','MarkerEdgeColor','k' ,'MarkerSize', 7);
box; axis tight; axis equal;
% grid off;

%set(gcf, 'Position', [50, 50, 800, 500]); % [left bottom width height]
if SAVE_FILE  
    filename = strcat(lower(obsFilename),'.eps');
    filepath = strcat(fp_figs, filename);
    saveas(gcf,filepath,'epsc');
end


% Graph Gt
figure
Gt.plot();
Gt.highlight_path(pathGE_Dks_,'EdgeColor','g', 'NodeSize', 7);
lmt = 0.3;
axis equal; xlim([-lmt virtual_goal+lmt]); ylim([-lmt virtual_goal+lmt]);
box;

% Compute friction 
[totalFrictionG_Dks, tensionPerSegment_Dks] = calcFrictionFromPath(Gt, pathG_Dks, sizeVisTable, edgeID);
[totalFrictionGE_Dks, tensionPerSegment_TAw] = calcFrictionFromPath(Gt, pathGE_Dks, sizeVisTable, edgeID);
[totalFrictionGE_BFS, tensionPerSegment_BFS] = calcFrictionFromPath(Gt, pathG_BFS, sizeVisTable, edgeID);

if (DATALOGGER_ON_DEBUG && DEBUG > 0)
    fclose(fid);
end


%%
% clc
fprintf('\nSIMULATION RESULTS\n')
fprintf('Graph Search\t| Cost\t| Length\t| Legend\n')
fprintf('----------------+-------+-----------+---------\n')
fprintf('BFS\t\t\t\t| %1.2f \t| %1.2fm \t| %s\n', totalFrictionGE_BFS, tetherLengthGE_BFS, 'blue')
fprintf('Dijkstra\t\t| %1.2f\t| %1.2fm \t| %s\n', totalFrictionG_Dks, tetherLengthG_Dks, 'red')
fprintf('Tension-aware\t| %1.2f\t| %1.2fm \t| %s\n', totalFrictionGE_Dks, tetherLengthGE_Dks, 'green')

% DEBUG LEVEL 1
if (DEBUG == 1)
    % These  values should be the same!
    fprintf('\n');
    fprintf('These cost function values (or friction) should be the same: [%1.5f | %1.5f]\n\n', ...
        costGE_Dks, totalFrictionGE_Dks);
end

%% Plot simulation results on top of real-world experiments
if (EXPERIMENT_SETUP ~= 0)
    % Call script "experimentSetup.m" to plot the image and solutions
    experimentalSetup2
end

function results = printExperimentalResults()
    % Friction estimation
    T0 = 0.5;    % Mass load weight [lbs]
    wrapAngle = pi; % Radians
    
    % White plastic cylinder covered with copper tape on surface
    T_plastic_cylinder = 0.90; % Pounds
    mu_plastic_cylinder = (log(T_plastic_cylinder) - log(T0))/wrapAngle;
    T_plastic_cylinder_check = T0 * exp(wrapAngle * mu_plastic_cylinder);
    
    % Cardboard cylinder 
    T_cardboard = 1.36; % Pounds
    mu_cardboard = (log(T_cardboard) - log(T0))/wrapAngle;
    T_cardboard_check = T0 * exp(wrapAngle * mu_cardboard);

    % Measures (trials split in rows)
    T_BFS = [1.10 1.14 1.13 1.14 1.25
             1.11 1.08 1.08 1.05 1.08];
    T_Dks = [1.07 1.10 1.08 1.05 1.03       % 1.37 1.30 1.48 1.35 1.33
             0.99 1.01 1.02 1.04 1.00];
    T_TAw = [0.82 0.84 0.88 0.91 0.97
             0.99 1.01 1.02 1.04 1.00];
    
    % \tau = \delta T = T_0 * (e^{\mu*\theta} - 1) = T - T_0
    deltaT_BFS = T_BFS - T0;
    deltaT_Dks = T_Dks - T0;
    deltaT_TAw = T_TAw - T0;

%     T_BFS = T_BFS/T0;   % Normalization
%     T_Dks = T_Dks/T0;   % Normalization
%     T_TAw = T_TAw/T0;   % Normalization

    % Experiment #1
    % Mean
    Tavg_BFS1 = mean(deltaT_BFS(1,:));
    Tavg_Dks1 = mean(deltaT_Dks(1,:));
    Tavg_TAw1 = mean(deltaT_TAw(1,:));
    % Standard Deviation
    Tstd_BFS1 = std(deltaT_BFS(1,:));
    Tstd_Dks1 = std(deltaT_Dks(1,:));
    Tstd_TAw1 = std(deltaT_TAw(1,:));
    % Length
    L_BFS1 = 2.95;
    L_Dks1 = 2.37;
    L_TAw1 = 2.39;
    
    % Experiment #2
    % Mean
    Tavg_BFS2 = mean(deltaT_BFS(2,:));
    Tavg_Dks2 = mean(deltaT_Dks(2,:));
    Tavg_TAw2 = mean(deltaT_TAw(2,:));
    % Standard Deviation
    Tstd_BFS2 = std(deltaT_BFS(2,:));
    Tstd_Dks2 = std(deltaT_Dks(2,:));
    Tstd_TAw2 = std(deltaT_TAw(2,:));
    % Length
    L_BFS2 = 3.13;
    L_Dks2 = 2.45;
    L_TAw2 = 2.45;
    
    fprintf('EXPERIMENTAL RESULTS\n')
    fprintf('Graph Search\t| Cost (std)\t| Length\t| Legend\n')
    fprintf('----------------+---------------+-----------+---------\n')
    fprintf('BFS\t\t\t\t| %1.2f (%1.2f) \t| %1.2fm \t| %s\n', Tavg_BFS1, Tstd_BFS1, L_BFS1, 'blue')
    fprintf('Dijkstra\t\t| %1.2f (%1.2f)\t| %1.2fm \t| %s\n', Tavg_Dks1, Tstd_Dks1, L_Dks1, 'red')
    fprintf('Tension-aware\t| %1.2f (%1.2f)\t| %1.2fm \t| %s\n', Tavg_TAw1, Tstd_TAw1, L_TAw1, 'green')
    fprintf('----------------+---------------+-----------+---------\n')
    fprintf('BFS\t\t\t\t| %1.2f (%1.2f) \t| %1.2fm \t| %s\n', Tavg_BFS2, Tstd_BFS2, L_BFS2, 'blue')
    fprintf('Dijkstra\t\t| %1.2f (%1.2f)\t| %1.2fm \t| %s\n', Tavg_Dks2, Tstd_Dks2, L_Dks2, 'red')
    fprintf('Tension-aware\t| %1.2f (%1.2f)\t| %1.2fm \t| %s\n', Tavg_TAw2, Tstd_TAw2, L_TAw2, 'green')
    
    results.mu.plastic_cylinder = mu_plastic_cylinder;
    results.mu.cardboard = mu_cardboard;
end
