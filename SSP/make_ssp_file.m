clear all;

%%
clear all;
% Define the directory containing your data files
Directory = 'C:\Users\alexh\Documents\UW_Acoustics\Codes\sbc-vla-data\SBCEX2022\CTD_info\CTD_files\';

% Get a list of all data files in the directory
dataFiles = dir(fullfile(Directory, '*.mat'));

% Initialize a cell array to store the data
SSP = cell(1, numel(dataFiles));
C_Z = cell(1, numel(dataFiles));

% Loop through each data file
for i = 1:numel(dataFiles)
    % Load the data from the file
    loadedData = load(fullfile(Directory, dataFiles(i).name),'ctd_Depth','ctd_ssp');

    % Assuming the data is stored in a variable called 'dataArray'
    ssp = loadedData.ctd_ssp;
    c_z = loadedData.ctd_Depth;
    [~,id]=max(c_z);
    ssp=ssp(1:id);
    c_z=c_z(1:id);
    ssp= [ssp(1),ssp,ssp(end)];
    c_z = [0,c_z,floor(c_z(end)+5)];
    % Store the data in a separate column of the cell array
    SSP{i} = ssp(:);  % Convert to a column vector
    C_Z{i} = c_z(:);
    % Optionally, you can display the size of each column
    
end

save('C:\Users\alexh\Documents\UW_Acoustics\Codes\sbc-vla-data\SBCEX2022\CTD_info\SBCEXP22_SSP.mat','C_Z','SSP')



