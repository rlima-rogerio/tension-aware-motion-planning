clc
clear
close all

% Save file: 0:disable | 1:enable
SAVE_FILE = true;

% Get the relative folder path to access the sub-folder 'figs'
cd(strcat('figs',filesep)); fp_figs = [pwd filesep]; cd ..;

% Load saved figures
c=hgload(strcat(fp_figs,'benchmark1.fig'));
k=hgload(strcat(fp_figs,'benchmark2.fig'));

% Prepare subplots
figure
h(1)=subplot(2,1,1); box;
h(2)=subplot(2,1,2); box;

% Paste figures on the subplots
copyobj(allchild(get(c,'CurrentAxes')),h(1)); axis equal;
copyobj(allchild(get(k,'CurrentAxes')),h(2)); axis equal;

axis(h(1), 'equal');
axis(h(2), 'equal');

set(h(1), 'xlim', [-2, 16]);
set(h(2), 'xlim', [-2, 16]);

% Add title
t(1) = title(h(1),'(a)');
t(2) = title(h(2),'(b)');

% Add legends
% l(1)=legend(h(1),'LegendForFirstFigure')
% l(2)=legend(h(2),'LegendForSecondFigure')

%set(gcf, 'Position', [50, 50, 800, 500]); % [left bottom width height]
if SAVE_FILE  
    filename = 'result-benchmark.eps';
    filepath = strcat(fp_figs, filename);
    saveas(gcf,filepath,'epsc');
end