% This script is written to convert .sio files to mat files
% by reading them with the sioread.m from Bill Hodgkiss
% then saving them to .mat file with builtin save()
% 
% This program is not very robust and just simply assumes
% you've pre organized data so that ONLY .sio files exist
% in the folder.

%%% USER INPUT %%%
% Folder to read files from
readpath = "F:\SCE17_Hodgkiss_VLA\acoustic\SCE17_VLA2_sio";
% Folder to save files to
savepath = "F:\SCE17_Hodgkiss_VLA\acoustic\SCE17_VLA2_mat";

% Define the file# to start on.
start = 1;
%%%%%%%%%%%%%%%%%%

%%% Workspace %%%
files = dir(readpath);
num_files = size(files,1);

% Load each file from the sio folder and save it as .mat file
for i = start:num_files
    i = i + 2
    f_read = readpath + '\' + files(i).name;
    
    % Change file names/loc to .mat folder and filename.
    name = strrep(files(i).name, '.sio', '.mat');
    f_save = savepath + '\' + name;
    disp(f_save)
    % Read the whole length, read all the channels
    data = sioread(f_read, 1, 0, 0);
    
    % Save to mat file
    save(f_save, 'data');
    
end


% % Randomly load in files from both the sio and npy folders and check to 
% % see if they are equivalent by subtracting them. Then they should output
% % zero.
% % Remember to change the folders that you read files from, to the npy
% % folder.
% for i = 1:50
%     z = rand(1);
%     z = z * num_files;
%     z = round(z);
%     disp(z)
%     
%     f_sio = readpath + '\' + files(z).name;
%     data_sio = sioread(f_sio, 1, 0, 0);
% 
%     name = strrep(files(z).name, '.sio', '.npy');
%     f_npy = savepath + '\' + name;
%     data_npy = readNPY(f_npy);
%     data_npy = cast(data_npy, 'double');
%     
%     output = data_sio - data_npy;
%     disp(sum(sum(output)))
% end

    