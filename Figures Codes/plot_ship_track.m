%
%  genereate a SBCEXP22 shallow water map,and overlay the Marina Dia
%  ship track on top.

clear all;
close all;
addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\sbc-vla-data\SBCEX2022\PLOT_SHIP_TRACKS\');
addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\sbc-vla-data\SBCEX2022\PLOT_SHIP_TRACKS\bathy');  % to plot base map
addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Notes\SBCEX Spreadsheets\mapping_spreadsheets\NORTH\');
addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Notes\SBCEX Spreadsheets\mapping_spreadsheets\SOUTH\');
addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Notes\SBCEX Spreadsheets\mapping_spreadsheets\OTHER\');
SBCEXP22_map_shallowWaterB


%addpath(C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\sbc-vla-data\SBCEX2022\PLOT_SHIP_TRACKS\distinguishable_colors
% Number of colors you want
numColors = 10;
addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\sbc-vla-data\SBCEX2022\PLOT_SHIP_TRACKS\distinguishable_colors\');

% Get distinguishable colors
colors = distinguishable_colors(numColors);


%title({'SBCEXP22 Deployment Locations with SIEM ARISTOTLE Ship Track', ...
%'Shown is MSC CANCUN Ship Track from JD145 19:50:55 - JD145 20:25:30'})


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

print('-dpng','C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\matlab-bborca\Inversion\Inversion_4000_iteration_scans\other_figs\Ship Tracks\Ship_Track.png')

