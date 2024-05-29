%
% m-file to find the depths of the MPL VLAs for the deep water site
load SBCEX21bathy.mat

imagesc(lon,lat,d)
set(gca,'ydir','normal')

%  now find depths of MPL VLAs
VLAlat = [40.048467, 39.953967 ] ; % VLA1, VLA2 Lat N
VLAlon = [-70.884167, -70.770750 ] ; % VLA1, VLA2 Lon W

VLAdepths = interp2(lon,lat,d,VLAlon,VLAlat);

fprintf(1,'High Resolution Bathymetry  VLA 1 depth: %5.1f\n',VLAdepths(1));
fprintf(1,'High Resolution Bathymetry  VLA 2 depth: %5.1f\n',VLAdepths(2));

