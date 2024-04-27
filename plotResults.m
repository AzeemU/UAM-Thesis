%close all; clear; clc
function [] = plotResults(Link_Flow, Community_Demand, fairType, fairVal, budget, Vertiport_Status)

%Default budget
if nargin < 5
    budget = 17;
    Vertiport_Status = ones(17,1);
end

set(0,'defaultTextInterpreter','latex');


%% Load Data

% Basic network information: nodes, edges, communities
Data_Locations = readtable('AustinUAMNetwork.xlsx');
Data_Edges = readtable('AustinUAMNetworkEdges.xlsx');
S = load('AustinZipCodes.mat').S;

crmname = 'tv';

% Use alpha = 1 results as the colormap baseline
Link_Flow_1 = readmatrix([crmname, '_link_alpha_1_delta_0.9.csv']);
% Route_Flow_1 = readmatrix('route_flow_alpha_1.csv');
Community_Demand_1 = readmatrix([crmname, '_community_alpha_1_delta_0.9.csv']);

% The current alpha data - use if want to run this as not a function
% Alpha = 2; % select the alpha value
% DeltaThresh = 1000;
% budget = 15; %max 17

% switch fairType
%     case 'alpha'
%         %For alpha fairness (if main.m was run)
%         Link_Flow = readmatrix([crmname, '_link_alpha_',num2str(fairVal),'_delta_0.9.csv']);
%         % Route_Flow = readmatrix(['route_flow_alpha_',num2str(Alpha),'.csv']);
%         Community_Demand = readmatrix([crmname,'_community_alpha_',num2str(fairVal),'_delta_0.9.csv']);
%     case 'threshold'
%         %For threshold fairness (if thresholdmain was run)
%         %Shouldn't really be used bc vertiport selection with a value of 17
%         %does the same thing
%         Link_Flow = readmatrix([crmname, '_link_thresh_',num2str(fairVal),'_delta_0.9.csv']);
%         % Route_Flow = readmatrix(['route_flow_alpha_',num2str(Alpha),'.csv']);
%         Community_Demand = readmatrix([crmname,'_community_thresh_',num2str(fairVal),'_delta_0.9.csv']);
%     case 'netdesign'
%         %For vertiportselectionmain.m (threshold + vertiport selection)
%         Link_Flow = readmatrix([crmname, '_link_thresh_',num2str(fairVal),'_selection_', num2str(budget), '_delta_0.9.csv']);
%         % Route_Flow = readmatrix(['route_flow_alpha_',num2str(Alpha),'.csv']);
%         Community_Demand = readmatrix([crmname,'_community_thresh_',num2str(fairVal),'_selection_', num2str(budget), '_delta_0.9.csv']);
%     otherwise
%         print('Not a valid case. Use alpha, threshold, or netdesign.')
% end


% normalize data for color coding
CD_Max = max(Community_Demand_1);
LF_Max = max(Link_Flow_1);

Community_Demand = Community_Demand/CD_Max;
Link_Flow = Link_Flow/LF_Max;

fig = figure('Position', [50 50 950 1000]);
t = tiledlayout(1,1);


%% Draw Zipcode Areas
ax1 = geoaxes(t);

Lon_Upper = -97.5704;
Lon_Lower = -98.0054;
Lat_Upper = 30.5716;
Lat_Lower = 30.1544;

% colormap for communities
cmap1 = flipud(pink);

for i = 1:height(S)

    Boundaries_X = S(i).X;
    Boundaries_Y = S(i).Y;
    
    shape = geopolyshape(Boundaries_Y,Boundaries_X);
    pg = geoplot(shape,'LineWidth',0.5); hold on
    color = 'k';
    pg.FaceColor = cmap1(min(256,max(floor(Community_Demand(i)*256),1)),:);
    pg.FaceAlpha = 1;
    pg.EdgeColor = color;

%     AvgPosition = mean(S(i).BoundingBox);
%     text(AvgPosition(2),AvgPosition(1),num2str(i),'Color','b','FontSize',14);
    %hold on
end

%% Draw UAM Network
% colormap for link flow
cmap2 = redblue(256);
%cmap2 = flipud(colormap(parula));

N_Locations = height(Data_Locations);
N_Routes = height(Data_Edges);

for i = 1:N_Routes/2

    Sub_Table = Data_Edges(i,:);

    Airport1 = Sub_Table.Vertiport1;
    Airport2 = Sub_Table.Vertiport2;
    Capacity = Sub_Table.Capacity;

    Location_Airport1 = Data_Locations(strcmp(Data_Locations.Name, Airport1),:);
    Location_Airport2 = Data_Locations(strcmp(Data_Locations.Name, Airport2),:);

    Airport1Lat = Location_Airport1.Latitude;
    Airport1Lon = Location_Airport1.Longitude;
    Airport2Lat = Location_Airport2.Latitude;
    Airport2Lon = Location_Airport2.Longitude;

    geoplot([Airport1Lat Airport2Lat],[Airport1Lon Airport2Lon],'Color',[cmap2(min(256,ceil(Link_Flow(i)*256+0.001)),:),0.8],'LineWidth', Capacity*1);
    hold on

%     AvgPosition = [(Airport1Lat+Airport2Lat)/2, (Airport1Lon+Airport2Lon)/2-0.006];
%     text(AvgPosition(1),AvgPosition(2),num2str(i),'Color','r','FontSize',18); hold on

end

%% Draw Vertiports
for i = 1:N_Locations

    Sub_Table = Data_Locations(i,:);

    Latitude = Sub_Table.Latitude;
    Longitude = Sub_Table.Longitude;
    Name = Sub_Table.Name;
    Capacity = Sub_Table.Capacity;
    if round(Vertiport_Status(i)) == 1
        portcolor = [82/256, 204/256, 85/256];
    elseif round(Vertiport_Status(i)) == 0
        portcolor = [227/256 108/256 5/256];
    end

    pg = geoplot(Latitude,Longitude,'o', 'Color', portcolor, 'MarkerFaceColor',portcolor, 'MarkerSize', Capacity*1.6); hold on
%     text(Latitude-0.025,Longitude-0.03,Name,'Color','r','FontSize',14); hold on
%     text(Latitude,Longitude-0.005,num2str(i),'Color','k','FontSize',16); hold on

end

geobasemap streets-light

geolimits([Lat_Lower Lat_Upper],[Lon_Lower Lon_Upper]);
set(gca, 'FontName', 'Times', 'FontSize', 22);
if strcmp(fairType, 'alpha')
    title(['$\alpha =~$', num2str(fairVal)], 'FontSize', 22);
elseif strcmp (fairType, 'netdesign') || strcmp(fairType, 'thresh')
    title(['$\Delta =~$', num2str(fairVal)], 'FontSize', 22);
    subtitle(sprintf('Vertiports: %d\nGini Coeff: %.3f', budget, ginicoeff(abs(Community_Demand))), 'FontSize', 18);
end

%% Colorbars
colormap(ax1, cmap1)
cb1 = colorbar(ax1);
cb1.Layout.Tile = 'south';
cb1.Label.String = 'Community Demand Served';
clim([0 CD_Max]);

ax2 = axes(t);
set(gca, 'FontName', 'Times', 'FontSize', 22);
colormap(ax2, cmap2)
ax2.Visible = 'off';
cb2 = colorbar(ax2);
clim([0 LF_Max]);
cb2.Layout.Tile = 'east';
cb2.Label.String = 'Link Traffic Volume';
cb2.TickLabelsMode = 'manual';
cb2.Ticks = [0 LF_Max];
cb2.TickLabels{1} = 'Low';
cb2.TickLabels{2} = 'High';

%fontname(fig,'Arial')
end