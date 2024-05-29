%% JASA Figures 

% This code is used to create the figures used in the JASA "Deep sediment 
% heterogeneity inferred using very low-frequency features from merchant
% ships" by Alexandra M. Hopps-McDaniel et. al.

% The figures are as follows
%   - beta_1,2 values for all ships
%   - DL1 Sound Speeds with error bars for all ships
%   - DL1 Sound Speeds on a map of the NEMP (need to fix a bug)
%   - DL1 Sound Speeds with source and reciever prior model for the top
%     layers of the NEMP

%% B12 singular points

SHIPS = {'CARMEN'	'MSC CANCUN'	'ALS APOLLO'	'PHOENIX ADMIRAL 1'	'PHOENIX ADMIRAL 2'	'MAERSK KENSINGTON'	'GSL DOROTHEA'	'DEE4 FIG'	'MSC DON GIOVANNI'	'WAN HAI 301'	'SEAMAX GREENWICH'	'EVER FASHION'	'MAERSK KINLOSS'	'ADVENTURE OF THE SEAS' 'ATLANTIC SAIL'	'HISTRIA GIADA'	'NORDBAY'	'MSC KATRINA'	'MAERSK KLEVEN'	'AURIGA LEADER'	'GSL NICOLETTA'	'SIEM ARISTOTLE'	'GLEN CANYON'	'TORRENTE'};
indices = 1:numel(SHIPS);

% 
b12=[23.4	23.8	23.7	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
24.4	24.6	23.9	24.1	24.7	24.7	23.9	24.3	24.3	24.3	24.5	24.2	24.2	24.6	24.8	24.3	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
25.7	26.2	25.8	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	27.2	NaN	27.7	27.2	NaN	NaN	NaN	NaN	NaN
NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	27.5	26.9	27.5	27	27.8	26.9	27.3	NaN	NaN	NaN
NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	27	27.3	NaN	NaN	NaN	NaN	NaN
NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	22.3	NaN	NaN
NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	23.5	23	23.3
NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	24.6	NaN	NaN];

figure;
shift_amount=0.1;

% southern shippign lane
scatter(indices-shift_amount, b12(1,:),'filled','MarkerFaceColor',[0.6350 0.0780 0.1840]);
hold on;
scatter(indices, b12(2,:), 'filled','MarkerFaceColor',[0.3010 0.7450 0.9330]);
scatter(indices+shift_amount, b12(3,:), 'filled','MarkerFaceColor',[0.4660 0.6740 0.1880]);
% northern shipping lane
scatter(indices-shift_amount, b12(4,:),'filled','Marker', 's','MarkerFaceColor',[0.6350 0.0780 0.1840]);
scatter(indices, b12(5,:),'filled','Marker', 's','MarkerFaceColor',[0.3010 0.7450 0.9330]);
scatter(indices+shift_amount, b12(6,:),'filled','Marker', 's','MarkerFaceColor',[0.4660 0.6740 0.1880]);
% other track
scatter(indices-shift_amount, b12(7,:),'filled','Marker', 'd','MarkerFaceColor',[0.6350 0.0780 0.1840]);
scatter(indices, b12(8,:),'filled','Marker', 'd','MarkerFaceColor',[0.3010 0.7450 0.9330]);
scatter(indices+shift_amount, b12(9,:),'filled','Marker', 'd','MarkerFaceColor',[0.4660 0.6740 0.1880]);

h1 = scatter(NaN, NaN, 100, 'filled', 'Marker', 'o', 'MarkerFaceColor', [0.6350 0.0780 0.1840]);
h2 = scatter(NaN, NaN, 100, 'filled', 'Marker', 'o', 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
h3 = scatter(NaN, NaN, 100, 'filled', 'Marker', 'o', 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
h4 = plot(NaN, NaN, 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k');
h5 = plot(NaN, NaN, 'Marker', 's', 'Color', 'k', 'MarkerFaceColor', 'k');
h6 = plot(NaN, NaN, 'Marker', 'd', 'Color', 'k', 'MarkerFaceColor', 'k');

hold off;

xlim([0.5, numel(SHIPS) + 0.5])
ylim([22 28])
ylabel('Frequency (Hz)')
title('$\beta_{1,2}$ Crossings', 'Interpreter', 'latex')
xticks(indices)
xticklabels(SHIPS)
grid on;

legend([h1, h2, h3, h4, h5, h6], 'VLA1', 'PROTEUS', 'VLA2', 'SOUTHERN LANE','NORTHERN LANE','OTHER TRACK')

%% DL1 peak SS with error bars

SS = [NaN	NaN	1821	NaN	NaN	NaN	NaN	NaN	NaN	1801	NaN	NaN	NaN	NaN	1814	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
NaN	NaN	1791	1788	1785	1814	1795	1805	1769	1791	1778	NaN	1782	1791	1795	1782	1772	1765	1765	NaN	NaN	NaN	NaN	NaN	NaN	NaN
NaN	NaN	1782	NaN	NaN	NaN	NaN	NaN	NaN	1776	NaN	NaN	NaN	NaN	1779	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1710	NaN	1729	NaN	1739	NaN
NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1729	1723	1716	1713	1713	1736	1713
NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1723	NaN	NaN	NaN	1726	NaN
1759	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
1762	1775	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1831	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
1759	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN];

SHIPS = {'SIEM ARISTOTLE'	'GLEN CANYON'	'ALS APOLLO'	'PHOENIX ADMIRAL 2'	'MAERSK KENSINGTON'	'MSC DON GIOVANNI'	'WAN HAI 301'	'PHOENIX ADMIRAL 1'	'ATLANTIC SAIL SOUTH'	'MSC CANCUN'	'HISTRIA GIADA SOUTH'	'TORRENTE'	'GSL DOROTHEA'	'DEE4 FIG'	'CARMEN'	'ADVENTURE OF THE SEAS'	'EVER FASHION'	'MAERSK KINLOSS'	'SEAMAX GREENWICH'	'AURIGA LEADER'	'MSC KATRINA'	'GSL NICOLETTA'	'HISTRIA GIADA NORTH'	'NORDBAY'	'MAERSK KLEVEN'	'ATLANTIC SAIL NORTH'};
indices = 1:numel(SHIPS);


errorbar =[1755.6	NaN	1820.1	NaN	NaN	NaN	NaN	NaN	NaN	1796.8	NaN	NaN	NaN	NaN	1806	NaN	NaN	NaN	NaN	NaN	1703.2	NaN	1721.2	NaN	1719.1	NaN
NaN	NaN	1785	1785.5	1782.1	1806.1	1794.1	1803.9	1762.8	1763.8	1778.3	NaN	1775.9	1790.3	1792.8	1741.2	1767.3	1763.7	1765.2	NaN	NaN	NaN	NaN	NaN	NaN	NaN
1758.5	1766.7	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1827.5	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1721.4	1711.3	1709.6	1706.4	1705.6	1719.2	1705
1758.5	NaN	1777.1	NaN	NaN	NaN	NaN	NaN	NaN	1764.3	NaN	NaN	NaN	NaN	1765.6	NaN	NaN	NaN	NaN	NaN	1711.5	NaN	NaN	NaN	1715.3	NaN
1764.4	NaN	1822.5	NaN	NaN	NaN	NaN	NaN	NaN	1805.6	NaN	NaN	NaN	NaN	1814.8	NaN	NaN	NaN	NaN	NaN	1710.4	NaN	1731	NaN	1742.1	NaN
NaN	NaN	1791.6	1790.7	1788.7	1826.9	1795.9	1807.3	1771.6	1801.4	1781.7	NaN	1784.1	1792.5	1795.4	1790.4	1776.3	1766.9	1765.4	NaN	NaN	NaN	NaN	NaN	NaN	NaN
1763.1	1777.5	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1830.7	NaN	NaN	NaN	NaN	NaN	NaN	NaN	1731	1724.7	1717.6	1713.4	1713.8	1741.2	1713
1760	NaN	1784.3	NaN	NaN	NaN	NaN	NaN	NaN	1794.3	NaN	NaN	NaN	NaN	1791.2	NaN	NaN	NaN	NaN	NaN	1724.3	NaN	NaN	NaN	1729.1	NaN];

figure;
shift_amount=0.1;

% southern shippign lane
scatter(indices-shift_amount, SS(1,:),'filled','MarkerFaceColor',[0.6350 0.0780 0.1840]);
hold on;
scatter(indices, SS(2,:), 'filled','MarkerFaceColor',[0.3010 0.7450 0.9330]);
scatter(indices+shift_amount, SS(3,:), 'filled','MarkerFaceColor',[0.4660 0.6740 0.1880]);
% northern shipping lane
scatter(indices-shift_amount, SS(4,:),'filled','Marker', 's','MarkerFaceColor',[0.6350 0.0780 0.1840]);
scatter(indices, SS(5,:),'filled','Marker', 's','MarkerFaceColor',[0.3010 0.7450 0.9330]);
scatter(indices+shift_amount, SS(6,:),'filled','Marker', 's','MarkerFaceColor',[0.4660 0.6740 0.1880]);
% other track
scatter(indices-shift_amount, SS(7,:),'filled','Marker', 'd','MarkerFaceColor',[0.6350 0.0780 0.1840]);
scatter(indices, SS(8,:),'filled','Marker', 'd','MarkerFaceColor',[0.3010 0.7450 0.9330]);
scatter(indices+shift_amount, SS(9,:),'filled','Marker', 'd','MarkerFaceColor',[0.4660 0.6740 0.1880]);

h1 = scatter(NaN, NaN, 100, 'filled', 'Marker', 'o', 'MarkerFaceColor', [0.6350 0.0780 0.1840]);
h2 = scatter(NaN, NaN, 100, 'filled', 'Marker', 'o', 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
h3 = scatter(NaN, NaN, 100, 'filled', 'Marker', 'o', 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
h4=plot(NaN, NaN, 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k');
h5=plot(NaN, NaN, 'Marker', 's', 'Color', 'k', 'MarkerFaceColor', 'k');
h6=plot(NaN, NaN, 'Marker', 'd', 'Color', 'k', 'MarkerFaceColor', 'k');


width=1;
for i = 1:numel(indices)
    x = indices(i);
        
    x_errorbar_VLA1 = [x, x];
    y_errorbar_VLA1 = [errorbar(1,i);errorbar(5,i)];
    
    x_errorbar_PROTEUS_S = [x, x];
    y_errorbar_PROTEUS_S = [errorbar(2,i);errorbar(6,i)];

    x_errorbar_PROTEUS_N_O = [x, x];
    y_errorbar_PROTEUS_N_O = [errorbar(3,i);errorbar(7,i)];
    
    x_errorbar_VLA2 = [x, x];
    y_errorbar_VLA2 = [errorbar(4,i);errorbar(8,i)];
    
    plot(x_errorbar_VLA1-shift_amount, y_errorbar_VLA1, 'Color', [0.6350 0.0780 0.1840],'LineWidth',width);
    plot(x_errorbar_PROTEUS_S, y_errorbar_PROTEUS_S, 'Color', [0.3010 0.7450 0.9330],'LineWidth',width);
    plot(x_errorbar_PROTEUS_N_O, y_errorbar_PROTEUS_N_O, 'Color', [0.3010 0.7450 0.9330],'LineWidth',width);
    plot(x_errorbar_VLA2+shift_amount, y_errorbar_VLA2, 'Color', [0.4660 0.6740 0.1880],'LineWidth',width);
    
    marker_length = 0.05;  % Adjust the length of the markers as needed
end

 hold off;

xlim([0.5, numel(SHIPS) + 0.5])
ylim([1700,1835])
ylabel('Sound Speed (m/s)')
title('Deep Layer 1 Peak Sound Speeds')
xticks(indices)
xticklabels(SHIPS)
grid on;

legend([h1, h2, h3, h4, h5, h6],{'VLA1','PROTEUS','VLA2', 'SOUTHERN LANE','NORTHERN LANE','OTHER TRACK'})


%% Plot SS on map of NEMP
addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\sbc-vla-data\SBCEX2022\PLOT_SHIP_TRACKS\');
addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\sbc-vla-data\SBCEX2022\PLOT_SHIP_TRACKS\bathy');  % to plot base map
load('bathy/SBCEX21bathy.mat')
%SBCEXP22_map_shallowWaterB
LONindx = find( -lon >= 70.45 & -lon <= 70.76 );
LATindx = find(  lat >= 40.41 &  lat <= 40.55 );
h1 = figure(1); clf
h1.Units = 'inches';
h1.Position = [2.5 2.8 7.33 8.0];

[~,h10] = contour(-lon(LONindx),lat(LATindx),-d(LATindx,LONindx),45:1:155,'LineColor',[0.8 0.8 0.8]);
set(gca,'xdir','reverse');

hold on;
set(h10,'LineWidth',0.5);
[C,h50] = contour(-lon,lat,-d,45:5:155,'LineColor',[0.8 0.8 0.8]);
set(h50,'LineWidth',2);
%colormap([0.8 0.8 0.8]);  % set color map to black(50%)

clabel(C,h50,[50 60 70 75 80 90 100 110 120 130 140 150],'FontSize',10,'Color',[.7 .7 .7]);

sLanes = [   40+33.5/60   70+49.0/60       40+34.0/60  70+24.0/60; ...
             40+31.5/60   70+49.0/60       40+32.0/60  70+24.0/60; ...
             40+25.5/60   70+49.0/60       40+26.0/60  70+24.0/60; ...
             40+23.5/60   70+49.0/60       40+24.0/60  70+24.0/60; ...
         ];

moorings = { '40' '28.207'  '70' '35.827' 'k' '^' 'VLA 1'
             '40' '27.547'  '70' '33.811' 'k' 's' 'PROTEUS'    %'PROTEUS ARL-UT'               
             '40' '26.507'  '70' '31.63'  'k' 'd' 'VLA 2'
            };

grid on
box on;


% plot westbound lanes
L1a = plot([sLanes(1,2) sLanes(1,4)],[sLanes(1,1) sLanes(1,3)],'color',[.5 .5 .5],'linewidth',2);
L1b = plot([sLanes(2,2) sLanes(2,4)],[sLanes(2,1) sLanes(2,3)],'color',[.5 .5 .5],'linewidth',2);
% plot eastbound lanes
L2a = plot([sLanes(3,2) sLanes(3,4)],[sLanes(3,1) sLanes(3,3)],'color',[.5 .5 .5],'linewidth',2);
L2b = plot([sLanes(4,2) sLanes(4,4)],[sLanes(4,1) sLanes(4,3)],'color',[.5 .5 .5],'linewidth',2);

axis([70.5 70.63 40.38 40.565]); grid on;
   set(gca,'dataAspectRatio',[1  cosd(40+30/60) 1 ])
   xlabel('Longitude (deg, aspect corrected for Latitude 40.5)')
   ylabel('Latitude (deg)')
   title('Deep Layer 1 Sound Speed at Discrete Coordinates')

   %  plot mooring
for k = 1:size(moorings,1)
  M = str2double(moorings(k,1:4));
  cval = char(moorings(k,5));
  sym = char(moorings(k,6));
  HM(k) = plot(M(3)+M(4)/60,M(1)+M(2)/60,'color',cval,'marker',sym,'markersize',6,'linewidth',1,'markerfacecolor',cval,'linestyle','none');
end


addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Notes\')

file = 'SBCEX_ALL_INTERESTED_SHIPS_FINAL.xlsx';

sheet = 'ALL SHIP INFO NEW';

lat_prot_range = 'O3:O28';
long_prot_range = 'P3:P28';
ss_prot_range = 'Q3:Q28';
lat_vla1_range = 'I3:I28';
long_vla1_range = 'J3:J28';
ss_vla1_range = 'K3:K28';
lat_vla2_range = 'U3:U28';
long_vla2_range = 'V3:V28';
ss_vla2_range = 'W3:W28';


lat_prot = xlsread(file,sheet,lat_prot_range);
long_prot = xlsread(file,sheet,long_prot_range);
ss_prot= xlsread(file,sheet,ss_prot_range);
lat_vla1 = xlsread(file,sheet,lat_vla1_range);
long_vla1 = xlsread(file,sheet,long_vla1_range);
ss_vla1 = xlsread(file,sheet,ss_vla1_range);
lat_vla2 = xlsread(file,sheet,lat_vla2_range);
long_vla2 = xlsread(file,sheet,long_vla2_range);
ss_vla2 = xlsread(file,sheet,ss_vla2_range);

clims=[1700,1835];
scatter(-long_prot,lat_prot,[],ss_prot,'filled')
 hold on;
scatter(-long_vla1,lat_vla1,[],ss_vla1,'filled')
scatter(-long_vla2,lat_vla2,[],ss_vla2,'filled')

colormap(jet)
caxis(clims)
colorbar;

h1=plot(NaN, NaN, 'Marker', '^', 'Color', 'k', 'MarkerFaceColor', 'k');
h2=plot(NaN, NaN, 'Marker', 's', 'Color', 'k', 'MarkerFaceColor', 'k');
h3=plot(NaN, NaN, 'Marker', 'd', 'Color', 'k', 'MarkerFaceColor', 'k');
h4=plot(NaN, NaN, 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k');

legend([h1,h2,h3,h4],'VLA1','PROTEUS','VLA2','SOUND SPEED')

%% Plot the 2022 SSP next to the 2017 SSP
wcpname  = 'SBCEX22_SSP';
wcp = load(wcpname);
ssp      = [2,3,6,11,12,14];

num_ssp  = length(ssp);
figure;
set(gca, 'YDir','reverse') 
data_17 = load('sbcexp2017_ssp.mat');
c = data_17.c_array(:,1);
c_z = data_17.c_z;
plot(c,c_z,'LineWidth',1,'Color','r')
hold on;
for i=1:num_ssp
    cur_ssp=ssp(i);
    c        = wcp.SSP{cur_ssp};         % m/s water SS array. rows: depth veriation; Columns: different SS profiles
    c_z      = wcp.C_Z{cur_ssp};
    plot(c,c_z,'LineWidth',1,'Color',[.5 .5 .5])
    set(gca, 'YDir','reverse') 
end 
h1 = plot(NaN, NaN, 'LineWidth',1,'Color','r');
h2 = plot(NaN, NaN, 'LineWidth',1, 'Color',[.5 .5 .5]);
%title('SBCEX Sound Speed Profiles')
ylabel('Depth (m)')
xlabel('Sound Speed (m/s)')
lgd=legend([h1, h2], '2017', '2022')

grid on;
xlim([1465 1510])
lgd.FontSize = 14;
ax = gca; 
ax.XAxis.FontSize = 14; 
ax.YAxis.FontSize = 14; 

%% Bounding the problem with the source and reciever DL1 SS estimation

SHIPS = {'ALS APOLLO' 'MSC DON GIOVANNI' 'ATLANTIC SAIL' 'HISTRIA GIADA' 'MAERSK KLEVEN'};
indices = 1:numel(SHIPS);

% original values
SS= [1791 1814 1769 1778 NaN; %s
    NaN NaN 1713 1713 1736;
    1814 1831 1791 1831 NaN; %south proteus
    NaN NaN 1752 1746 1775; %north proteus
    1821 NaN NaN NaN NaN; %s
    NaN NaN NaN NaN 1739;
    1847 NaN NaN NaN NaN; %south vla1
    NaN NaN NaN NaN 1765; %north vla1
    1782 NaN NaN NaN NaN;%s
    NaN NaN NaN NaN 1778;
    1778 NaN NaN NaN NaN; %south vla2
    NaN NaN NaN NaN 1762]; %north vla2


figure;
shift_amount=0.3;

scatter(indices, SS(1,:), 60, 'Marker', 'o','MarkerFaceColor','k','MarkerEdgeColor',[0.3010 0.7450 0.9330],'LineWidth',2);
hold on;
scatter(indices, SS(2,:), 60, 'Marker', 's','MarkerFaceColor','k','MarkerEdgeColor',[0.3010 0.7450 0.9330],'LineWidth',2);
scatter(indices, SS(3,:), 60, 'Marker', 'o','MarkerEdgeColor',[0.3010 0.7450 0.9330],'LineWidth',2);
scatter(indices, SS(4,:), 60, 'Marker', 's','MarkerEdgeColor',[0.3010 0.7450 0.9330],'LineWidth',2);

for i=1:5
    x=[i,i];
    y=[SS(1,i),SS(3,i)];
    line(x, y, 'Color',[0.3010 0.7450 0.9330], 'LineWidth', 1);
end

for i=1:5
    x=[i,i];
    y=[SS(2,i),SS(4,i)];
    line(x, y, 'Color',[0.3010 0.7450 0.9330], 'LineWidth', 1);
end

scatter(indices-shift_amount, SS(5,:), 60, 'Marker', 'o','MarkerFaceColor','k','MarkerEdgeColor',[0.6350 0.0780 0.1840],'LineWidth',2);
scatter(indices-shift_amount, SS(6,:), 60, 'Marker', 's','MarkerFaceColor','k','MarkerEdgeColor',[0.6350 0.0780 0.1840],'LineWidth',2);
scatter(indices-shift_amount, SS(7,:), 60, 'Marker', 'o','MarkerEdgeColor',[0.6350 0.0780 0.1840],'LineWidth',2);
scatter(indices-shift_amount, SS(8,:), 60, 'Marker', 's','MarkerEdgeColor',[0.6350 0.0780 0.1840],'LineWidth',2);
for i=1:5
    x=[i-shift_amount,i-shift_amount];
    y=[SS(5,i),SS(7,i)];
    line(x, y, 'Color',[0.6350 0.0780 0.1840], 'LineWidth', 1);
end

for i=1:5
    x=[i-shift_amount,i-shift_amount];
    y=[SS(6,i),SS(8,i)];
    line(x, y, 'Color',[0.6350 0.0780 0.1840], 'LineWidth', 1);
end
scatter(indices+shift_amount, SS(9,:), 60, 'Marker', 'o','MarkerFaceColor','k','MarkerEdgeColor',[0.4660 0.6740 0.1880],'LineWidth',2);
scatter(indices+shift_amount, SS(10,:), 60, 'Marker', 's','MarkerFaceColor','k','MarkerEdgeColor',[0.4660 0.6740 0.1880],'LineWidth',2);
scatter(indices+shift_amount, SS(11,:), 60, 'Marker', 'o','MarkerEdgeColor',[0.4660 0.6740 0.1880],'LineWidth',2);
scatter(indices+shift_amount, SS(12,:), 60, 'Marker', 's','MarkerEdgeColor',[0.4660 0.6740 0.1880],'LineWidth',2);

for i=1:5
    x=[i+shift_amount,i+shift_amount];
    y=[SS(9,i),SS(11,i)];
    line(x, y, 'Color',[0.4660 0.6740 0.1880], 'LineWidth', 1);
end

for i=1:5
    x=[i+shift_amount,i+shift_amount];
    y=[SS(10,i),SS(12,i)];
    line(x, y, 'Color',[0.4660 0.6740 0.1880], 'LineWidth', 1);
end

h1 = scatter(NaN, NaN, 100, 'Marker', 'o', 'MarkerEdgeColor', [0.6350 0.0780 0.1840],'LineWidth',2);
h2 = scatter(NaN, NaN, 100, 'Marker', 'o', 'MarkerEdgeColor', [0.3010 0.7450 0.9330],'LineWidth',2);
h3 = scatter(NaN, NaN, 100, 'Marker', 'o', 'MarkerEdgeColor', [0.4660 0.6740 0.1880],'LineWidth',2);
h4=plot(NaN, NaN, 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k');
h5=plot(NaN, NaN, 'Marker', 's', 'Color', 'k', 'MarkerFaceColor', 'k');
h6=plot(NaN, NaN, 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k');
h7=plot(NaN, NaN, 'Marker', 'o', 'Color', 'k');


xlim([0.5, numel(SHIPS) + 1])
ylim([1700,1870])
ylabel('Sound Speed (m/s)')
title('Inferred DL1 Sound Speeds')
xticks(indices)
xticklabels(SHIPS)
grid on;

lgd=legend([h1, h2, h3, h4, h5],{'VLA1','PROTEUS','VLA2', 'SOUTHERN LANE','NORTHERN LANE'});

ax = gca; 
ax.XAxis.FontSize = 14; 
ax.YAxis.FontSize = 14; 
lgd.FontSize = 12;
