%
%  plot map of shallow water site
%
%   M-file to plot SBCEXP22  shipping tracks, OP area, and
%        VLA and Temp array deployments.
%
%   reduce the lat/lon of the map to speed up the contour plot
%   generation routine.

%clear all ; % comment this out since M-file run in middle of other M-files
%close all
%
load('bathy/SBCEX21bathy.mat')



LONindx = find( -lon >= 70.45 & -lon <= 70.76 );
LATindx = find(  lat >= 40.41 &  lat <= 40.55 );

  h1 = figure(1); clf
  h1.Units = 'inches';
  h1.Position = [2.5 2.8 7.33 8.0];


  [~,h10] = contour(-lon(LONindx),lat(LATindx),-d(LATindx,LONindx),[45:1:155]);
  set(gca,'xdir','reverse');


  hold on
    set(h10,'LineWidth',0.5);
    [C,h50] = contour(-lon,lat,-d,45:5:155);
    set(h50,'LineWidth',2);
    colormap([0.8 0.8 0.8]);  % set color map to black(50%)
    clabel(C,h50,[50 60 70 75 80 90 100 110 120 130 140 150],'FontSize',10,'Color',[.7 .7 .7]);

%   OP area
OParea = [40.4792    	 70.5875	40 28.750	 70 35.250; ... % center
          40.5250     	 70.7667	40 31.500	 70 46.000; ... % 00  NW
          40.4333     	 70.7667	40 26.000	 70 46.000; ... % 01  SW
	  40.4333     	 70.4250	40 26.000	 70 25.500; ... % 02  SE
          40.5250     	 70.4250	40 31.500	 70 25.500; ... % 03  NE
         ];

% Shipping lanes:
%              western point            eastern point
sLanes = [   40+33.5/60   70+49.0/60       40+34.0/60  70+24.0/60; ...
             40+31.5/60   70+49.0/60       40+32.0/60  70+24.0/60; ...
             40+25.5/60   70+49.0/60       40+26.0/60  70+24.0/60; ...
             40+23.5/60   70+49.0/60       40+24.0/60  70+24.0/60; ...
         ];


%          color            LineStyle       LineWidth  marker markerSize
cListF = { '[     0 0.4470 0.7410]',   '-','0.5','.','6' ; ...  % Dk Blue
           '[     0 0.4470 0.7410]','none','1.0','.','6' ; ...  % Dk Blue
           '[0.8500 0.3250 0.0980]',   '-','0.5','.','6' ; ...  % Orange
           '[0.8500 0.3250 0.0980]','none','1.0','.','6' ; ...  % Orange
           '[0.9290 0.6940 0.1250]',   '-','0.5','.','6' ; ...  % golden rod
           '[0.9290 0.6940 0.1250]','none','1.0','.','6' ; ...  % golden rod
           '[0.4940 0.1840 0.5560]',   '-','0.5','.','6' ; ...  % purple
           '[0.4940 0.1840 0.5560]','none','1.0','.','6' ; ...  % purple
           '[0.4660 0.6740 0.1880]',   '-','0.5','.','6' ; ...  % green
           '[0.4660 0.6740 0.1880]','none','1.0','.','6' ; ...  % green
           '[0.3010 0.7450 0.9330]',   '-','0.5','.','6' ; ...  % lt blue
           '[0.3010 0.7450 0.9330]','none','1.0','.','6' ; ...  % lt blue
           '[0.6350 0.0780 0.1840]',   '-','0.5','.','6' ; ...  % rust
           '[0.6350 0.0780 0.1840]','none','1.0','.','6' ; ...  % rust
         };

	
moorings = { '40' '28.207'  '70' '35.827' '[.64 .08 .13]' '^' 'VLA 1'
             '40' '27.547'  '70' '33.811' '[.30 .75 .93]' 's' 'PROTEUS'    %'PROTEUS ARL-UT'               
             '40' '26.507'  '70' '31.63'  '[.47 .67 .18]' '^' 'VLA 2'
            };

%   SBCEXP 2022 Waypoints
%   Mud Patch (Shallow Site) Waypoints
waypoints = { ...
                '40' '30.288' '70' '40.541' 'WP1'; ...
                '40' '25.583' '70' '29.310' 'WP2'; ...
                '40' '30.583' '70' '29.054' 'WP3  '; ...
                '40' '26.065' '70' '36.500' 'WP4'; ...
%        Deep Site Waypoints
                '40' '00.976' '70' '39.842' 'WP7'; ...
                '40' '03.394' '70' '56.370' 'WP8'; ...
                '40' '06.738' '70' '54.631' 'WP9'; ...
                '39' '56.977' '70' '46.154' 'WP10'; ...
                '39' '57.635' '70' '58.976' 'WP11'; ...
            };
%        ********************************
%         Second Mud Patch Waypoints (between shallow and deep sites -
%         not part of AR67C)
%               '40' '16.483' '70' '34.607' 'WP5'; ...
%               '40' '08.117' '70' '21.039' 'WP6'; ...


             
%figure(1); clf
hold on
box on

%  plot operation area box
%B1 = plot(OParea([2:5 2],2),OParea([2:5 2],1),'k','linewidth',3);
hold on
set(gca,'xdir','reverse');
grid on
%axis(axisPosit)

% plot westbound lanes
L1a = plot([sLanes(1,2) sLanes(1,4)],[sLanes(1,1) sLanes(1,3)],'color',[.5 .5 .5],'linewidth',2);
L1b = plot([sLanes(2,2) sLanes(2,4)],[sLanes(2,1) sLanes(2,3)],'color',[.5 .5 .5],'linewidth',2);
% plot eastbound lanes
L2a = plot([sLanes(3,2) sLanes(3,4)],[sLanes(3,1) sLanes(3,3)],'color',[.5 .5 .5],'linewidth',2);
L2b = plot([sLanes(4,2) sLanes(4,4)],[sLanes(4,1) sLanes(4,3)],'color',[.5 .5 .5],'linewidth',2);


%  plot mooring
for k = 1:size(moorings,1)
  M = str2double(moorings(k,1:4));
  cval = char(moorings(k,5));
  sym = char(moorings(k,6));
  HM(k) = plot(M(3)+M(4)/60,M(1)+M(2)/60,'color',cval,'marker',sym,'markersize',6,'linewidth',1,'markerfacecolor',cval,'linestyle','none');
end


%hT = title('SBCEXP22 Ship Tracks');


   axis([70.4 70.85 40.33 40.605]); grid on;
   set(gca,'dataAspectRatio',[1  cosd(40+30/60) 1 ])
   xlabel('Longitude (deg, aspect corrected for Latitude 40.5)')
   ylabel('Latitude (deg)')

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %  add in 1 km arrows
  oLon = 70.7; % bottom of plot Latitude
  oLat = 40.50;  % left side of plot Longitude

 if ( 1 )
    lat = 40.5; % degrees Latitude to use for the calculation
    lat_km_min = 1.852;   % km / min  ( 1/lat_km_min = mins in km)
    lon_km_min = cosd(lat)*lat_km_min; % km / min
    lat_deg_km = 1/(lat_km_min*60); %  deg/km
    lon_deg_km = 1/(lon_km_min*60);     %  deg/km
   %
    plot([oLon-.001 oLon-0.001-(2*lon_deg_km)],[oLat+0.001 oLat+0.001],'k|-','linewidth',2')
    plot([oLon-.001 oLon-0.001-(4*lon_deg_km)],[oLat+0.001 oLat+0.001],'k|-','linewidth',2')
    plot([oLon-.001 oLon-0.001-(6*lon_deg_km)],[oLat+0.001 oLat+0.001],'k|-','linewidth',2')
    plot([oLon-.001 oLon-0.001-(8*lon_deg_km)],[oLat+0.001 oLat+0.001],'k|-','linewidth',2')
 %   plot([oLon-.001 oLon-0.001],[oLat+0.001 oLat+0.001+8*lat_deg_km],'k+-','linewidth',2')
    text(oLon-.001,oLat+0.001,'0  ','horizontalalignment','right','verticalalignment','middle')
    text(oLon-.000-8*lon_deg_km,oLat,'  8 km','horizontalalignment','left','verticalalignment','middle')
 end



hLeg = legend([L1a,HM(1),HM(2),HM(3)],'Shipping', ...
 char(moorings(1,7)), char(moorings(2,7)), ...
 char(moorings(3,7)), 'location','west');
hLeg.ItemTokenSize = [15 18];
hLeg.FontSize = 8;


%hT.Units = 'inches';
pause(1)
%hT.Position(2) = hT.Position(2) + .1;
drawnow


orient tall
print('-dpng','C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\sbc-vla-data\SBCEX2022\PLOT_SHIP_TRACKS\plots\SBCEXP22_shallowWater_map_NOwaypointsBtest.png')

