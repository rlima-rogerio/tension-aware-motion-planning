% =========================================================================
% WHO:      Rogerio Lima
% WHEN      August/2020 
% WHERE:    WVU - Morgantown
% =========================================================================

%% Pre set-up
clc
clear all
close all

% Global variables
global vectorObs;

%% Code

% Pops up a window querying the dimensions of a rectangular environment
prompt = {'Enter width size (m):','Enter height size (m):'};
dlgtitle = 'Rectangular Environment';
dims = [1 50];
definput = {'5','5'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
env_width_in = str2double(answer{1}); 
env_length_in = str2double(answer{2}); 

wall_thickness = 0.02*(env_width_in + env_length_in)/2;  % Thickness
env_width_out = env_width_in + wall_thickness;
env_length_out = env_length_in + wall_thickness;
x_in = [0 env_width_in env_width_in 0];
y_in = [0 0 env_length_in env_length_in];
x_out = [-wall_thickness env_width_out env_width_out -wall_thickness];
y_out = [-wall_thickness -wall_thickness env_length_out env_length_out];


% Takes the input data and draws the limit of the environment
env_boundaries = polyshape({x_in,x_out},{y_in,y_out});
plot(env_boundaries); hold on;
pbaspect([1 1 1]);

% Stores the C-space boundary
vectorObs{1} = env_boundaries;

% Initialize the state of the FSM
stateObstacles = 'queryObsVertices';
ind = 2;

while (true)
    switch stateObstacles
        case 'getObsVertices' 
            [xobs,yobs] = ginput;
            stateObstacles = 'checkObsVertices';
            
        case 'checkObsVertices'
            if (length(xobs) > 2)
                stateObstacles = 'createPolygon';
            else
                stateObstacles = 'invalidObs';
            end
            
        case 'createPolygon'
            obs = polyshape(xobs, yobs);
            vectorObs{ind} = obs;
            ind = ind + 1;
            stateObstacles = 'plotObsVertices';
            
        case 'saveObsVertices'
            save('Obstacles','vectorObs');
            stateObstacles = 'getObsVertices';
            % Finishes the obstacle creation
            break;
            
        case 'plotObsVertices'
            plot(obs);
            stateObstacles = 'queryObsVertices';
            
        case 'queryObsVertices'
            createNewObs = questdlg('Add a new obstacle?', ...
                                    'Obstacle creation', ...
                              'Yes','No','opts');
            if (strcmp(createNewObs,'Yes'))
                stateObstacles = 'getObsVertices';
            else
                stateObstacles = 'saveObsVertices';
            end
            
       case 'invalidObs'
            f = msgbox('Please, try again...', 'Warning','warn');
            uiwait(f);
            stateObstacles = 'getObsVertices';
        otherwise
            f = msgbox('Unexpected error...', 'Error','error');
            uiwait(f);
            stateObstacles = 'getObsVertices';
    end % switch-case
end % while


% Initialize the state of the FSM for C-space checking
stateCheckCSpace = 'queryChecking';

while (true)
    switch stateCheckCSpace
        case 'queryChecking' 
            createNewObs = questdlg('Do you want to check the C-space?', ...
                                    'C-space checking', ...
                              'Yes','No','opts');
            if (strcmp(createNewObs,'Yes'))
                stateCheckCSpace = 'getConfigurationQ';
            else
                stateCheckCSpace = 'exitChecking';
            end
            
        case 'getConfigurationQ'
            q = ginput(1);
            stateCheckCSpace = 'checkIsFree';
            
        case 'checkIsFree'
            ret = isFree(q);
            if(ret)
                stateCheckCSpace = 'spaceFree';     % true
            else
                stateCheckCSpace = 'spaceNotFree';  % false
            end
            
        case 'spaceFree'
            f = msgbox('C-space is Free!', 'Result','warn');
            uiwait(f);
            stateCheckCSpace = 'queryChecking';            
            
        case 'spaceNotFree'
            f = msgbox('It is inside of an obstacle...', 'Result','error');
            uiwait(f);
            stateCheckCSpace = 'queryChecking';

        case 'exitChecking'
            break;        
            
        otherwise
            f = msgbox('Unexpected error...', 'Error','error');
            uiwait(f);
            stateCheckCSpace = 'queryChecking';
    end % switch-case
end % while


%{
Refs:
https://www.mathworks.com/help/matlab/ref/msgbox.html
https://www.mathworks.com/help/matlab/ref/polyshape.html
https://www.mathworks.com/help/matlab/ref/inpolygon.html
%}

