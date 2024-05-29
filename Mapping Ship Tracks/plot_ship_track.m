% This code makes a plot of the ships from SBCEX 2022 in reference to the
% shipping lanes and line arrays that were deployed. 

% The base map and code comes from Sripps. 
% The additional code was written by Alexandra M. Hopps-McDaniel

clear all; close all;
currentFilePath = mfilename('fullpath');

[currentFolderPath, ~, ~] = fileparts(currentFilePath);

% Add path to distinguishable colors
colors_path = fullfile(currentFolderPath,'distinguishable_colors');
addpath(colors_path) 

% Add path to bathy
bathy_path = fullfile(currentFolderPath,'bathy');
addpath(bathy_path) 

% Add path to distinguishable colors
colors_path = fullfile(currentFolderPath,'Ship coordinates spreadsheets');
addpath(colors_path) 

SBCEXP22_map_shallowWaterB

% Number of colors you want
numColors = 10;


% Get distinguishable colors
colors = distinguishable_colors(numColors);


%title({'SBCEXP22 Deployment Locations with SIEM ARISTOTLE Ship Track', ...

%%
ship_coordinates = readtable('CARMEN_SOUTH.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;


plot(-lon,lat,'color',colors(2,:),'linewidth',1.3,'DisplayName','CARMEN')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'>','color',[0 0 0],'markerfacecolor',colors(2,:),'DisplayName','')

%%
ship_coordinates = readtable('MSC_CANCUN_SOUTH.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(1,:),'linewidth',1.3,'DisplayName','MSC CANCUN')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'>','color',[0 0 0],'markerfacecolor',colors(1,:),'DisplayName','')

%%
ship_coordinates = readtable('ALS_APOLLO_SOUTH.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;

set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(3,:),'linewidth',1.3,'DisplayName','ALS APOLLO')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'>','color',[0 0 0],'markerfacecolor',colors(3,:),'DisplayName','')

%%
ship_coordinates = readtable('PHOENIX_ADMIRAL_SOUTH1.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(4,:),'linewidth',1.3,'DisplayName','PHOENIX ADMIRAL 1')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'>','color',[0 0 0],'markerfacecolor',colors(4,:),'DisplayName','')

%%
ship_coordinates = readtable('PHOENIX_ADMIRAL_SOUTH2.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(5,:),'linewidth',1.3,'DisplayName','PHOENIX ADMIRAL 2')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'>','color',[0 0 0],'markerfacecolor',colors(5,:),'DisplayName','')

%%
ship_coordinates = readtable('MAERSK_KENSINGTON_SOUTH.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(6,:),'linewidth',1.3,'DisplayName','MAERSK KENSINGTON')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'>','color',[0 0 0],'markerfacecolor',colors(6,:),'DisplayName','')

%%
ship_coordinates = readtable('GSL_DOROTHEA_SOUTH.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(7,:),'linewidth',1.3,'DisplayName','GSL DOROTHEA')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'>','color',[0 0 0],'markerfacecolor',colors(7,:),'DisplayName','')

%%
ship_coordinates = readtable('DEE4_FIG_SOUTH.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;

set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(8,:),'linewidth',1.3,'DisplayName','DEE4 FIG')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'>','color',[0 0 0],'markerfacecolor',colors(8,:),'DisplayName','')
%%
ship_coordinates = readtable('ATLANTIC_SAIL_SOUTH.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(9,:),'linewidth',1.3,'DisplayName','ATLANTIC SAIL SOUTH')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'>','color',[0 0 0],'markerfacecolor',colors(9,:),'DisplayName','')

%%
ship_coordinates = readtable('ATLANTIC_SAIL_NORTH.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(9,:),'linewidth',1.3,'LineStyle','--', 'DisplayName','ATLANTIC SAIL NORTH')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'<','color',[0 0 0],'markerfacecolor',colors(9,:),'DisplayName','')



%%
ship_coordinates = readtable('HISTRIA_GIADA_SOUTH.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(10,:),'linewidth',1.3,'DisplayName','HISTRIA GIADA SOUTH')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'>','color',[0 0 0],'markerfacecolor',colors(10,:),'DisplayName','')

%%
ship_coordinates = readtable('HISTRIA_GIADA_NORTH.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(10,:),'linewidth',1.3,'LineStyle','--', 'DisplayName','HISTRIA GIADA NORTH')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'<','color',[0 0 0],'markerfacecolor',colors(10,:),'DisplayName','')

%%
ship_coordinates = readtable('NORDBAY_NORTH.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(1,:),'linewidth',1.3,'DisplayName','NORDBAY','LineStyle','--')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'<','color',[0 0 0],'markerfacecolor',colors(1,:),'DisplayName','')

%%
ship_coordinates = readtable('MSC_KATRINA_NORTH.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(2,:),'linewidth',1.3,'LineStyle','--','DisplayName','MSC KATRINA')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'<','color',[0 0 0],'markerfacecolor',colors(2,:),'DisplayName','')

%%
ship_coordinates = readtable('MAERSK_KLEVEN_NORTH.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(7,:),'linewidth',1.3,'LineStyle','--','DisplayName','MAERSK KLEVEN')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'<','color',[0 0 0],'markerfacecolor',colors(7,:),'DisplayName','')

%%
ship_coordinates = readtable('SIEM_ARISTOTLE_OTHER.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(4,:),'linewidth',1.3,'LineStyle',':', 'DisplayName','SIEM ARISTOTLE')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'>','color',[0 0 0],'markerfacecolor',colors(4,:),'DisplayName','')

%%
ship_coordinates = readtable('GLEN_CANYON_OTHER.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(5,:),'linewidth',1.3,'LineStyle',':','DisplayName','GLEN CANYON')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'>','color',[0 0 0],'markerfacecolor',colors(5,:),'DisplayName','')

%%
ship_coordinates = readtable('TORRENTE_OTHER.csv');

lat = ship_coordinates.LAT;
lon = ship_coordinates.LONG;
set(hL,'AutoUpdate','on')
plot(-lon,lat,'color',colors(8,:),'linewidth',1.3,'LineStyle',':','DisplayName','TORRENTE')
hL = legend;
set(hL,'AutoUpdate','off')
plot(-lon(end),lat(end),'<','color',[0 0 0],'markerfacecolor',colors(8,:),'DisplayName','')


%%
name='SBCEX 2022 Ship Tracks.png';
path=fullfile(currentFolderPath,name);
save
print('-dpng',path)

