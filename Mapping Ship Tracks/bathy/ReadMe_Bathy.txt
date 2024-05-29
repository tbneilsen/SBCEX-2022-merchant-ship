               ReadMe_Bathy.txt

David Dall'Osto generated the SBCEXP21.mat bathymetry file.

It is parsed from the NOAA data found here:

   https://www.ngdc.noaa.gov/mgg/coastal/crm/data/arc_ascii/

and put into a *.mat file. The Resolution is 3 seconds (0.0008333 degrees).
The file name for the data used is 'se_atl_crm_v1.asc' (available as either
a zip file or gzip file).


At 40-27.5N, 70-34.0W (approximately the cross-over of the two main tracks
in the Mud Patch):

lat: 1.8507 km/min or 92.5 m / 3 sec

lon: 1.4137 km/min or 70.7 m / 3 sec


Below is a snippet of matlab code David used to extract the data from the
*.asc file.  You can glean information from the header on the size of the
map, and coordinates of the lower-left corner of the grid.



%CODE SNIPPET%
filename = 'se_atl_crm_v1.asc';
FID = fopen(filename);
textread = [];
header = 'header info';
for m= 1:6
   TLINE = fgetl(FID);
   header = [header ', ' TLINE];
end

d = rand(18001,13201); %details from header file

for m = 1:13201
    TLINE = str2num(fgetl(FID));
    d(:,m) = TLINE(1:1:end);
    if rem(m,100)==0;
        fprintf([int2str(m/100) ' percent done \n'])
    end
end
