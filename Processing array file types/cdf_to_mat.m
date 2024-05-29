%% cdf to .mat

%   This code converts the ARL:UT (PROTEUS) 5.5-min CDF files to 5.5-min
%   .mat files. 

%   Author: Alexandra Hopps McDaniel

%% User defined variables 

% Add the data path and also path to the folder where the new files should be saved. 
cdf_filedir = 'D:\Proteus\CDF\processing\';
save_path = 'D:\Proteus\CDF\processed\';

%%
MyFolderInfo = dir(cdf_filedir);

for i=3:length(MyFolderInfo)
    file = MyFolderInfo(i).name;
    
    data = PROTEUS_ReadData_Function(cdf_filedir,file);
    fs=data.fs;
    t_unixsec=data.t_unixsec;
    hyd_data=data.hyd_data;

    save_name = append(save_path,file,'.mat');
    save(save_name,"fs","t_unixsec","hyd_data")
end





