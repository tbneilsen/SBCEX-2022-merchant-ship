                                                     %%  Save 1-min files 

%   This code converts the 5.5-min .mat ARL:UT VLA/HLA (PROTEUS) data into
%   1-min .mat files. 

%   This code loads the data from a folder or directory in 10-file chunks and
%   splits the approximately 55 minutes of data into 1-min files. 
%   Each file contains data from all 52 channels as well as the time step in
%   unix time and "special" time stamp form, which is a double of the form
%   YYDDDHHMMSS.SSSSSSS, where DDD is the Julian day. 

%   In order to not miss any minutes of data, each iteration, this code
%   loads in the last 5.5-min file from the previous iteration, so that
%   each file ends at a whole minute. (Ex: PROTEUS.22127092200.mat)
%   As a result, each time this script is called, it is necessary to include 
%   the last file from the previous run in the new folder.   
%   

%   Author: Alexandra Hopps McDaniel
%   Created: 8/15/23
%   Modified: 10/4/23
%%
clear all;
close all;
%% User defined variables 

data_path = 'F:\Proteus\2022 Data\5.5-min files\CURRENT\';
save_path = 'G:\PROTEUS\2022 Data\1-min Files\';

%% Defined Variables

addpath(genpath(data_path)) 

% Creates list of files in the folder that will be processed
directory = data_path;
MyFolderInfo = dir(directory);
MyFolderInfo = MyFolderInfo(3:end); % The first two indices for 'MyFolderInfo' are special directory references and can be ignored
files_in_dir = length(MyFolderInfo); 

repeat_process = false; 
num_groups = floor(files_in_dir / 10);

if num_groups >= 1
    repeat_process = true;
end 

k=1; % will use this var later for iterating though a loop
fs = 10240; % sampling frequency 
channels = 52; 
size_5min_data = 3353600; % size of data in the 5.5-min .mat files 
size_1min_data = fs*60; % size of data to be saved in the 1-min .mat files 
file_prefix = 'PROTEUS.'; % Prefix of the saved files 

%% This case is for when there are less than 10 files in the folder
if repeat_process==false
    % Preallocate the data arrays
    combined_dataSize = size_5min_data * files_in_dir;
    hyd_Data = zeros(channels, combined_dataSize);
    unix_Data = zeros(1, combined_dataSize);
    special_Data = zeros(1, combined_dataSize);
    num_min_files = floor(combined_dataSize / size_1min_data)-1;
    starting_indexes = zeros(1,num_min_files);

    % Load and combine data from multiple files
    for i = 1:files_in_dir
        % Load your data from each file
        file = MyFolderInfo(i).name;
        data = load(file,'hyd_data', 't_unixsec');
        % Add the data to the preallocated arrays
        hyd_Data(:, (i-1) * size_5min_data + 1 : i * size_5min_data) = data.hyd_data;
        unix_Data((i-1) * size_5min_data + 1 : i * size_5min_data) = data.t_unixsec;
    end
    
    % Allocate the special_Data array with the time converted from unix
    % time to special time (YYDDDHHMMSS, with DDD being the Julian day). 
    % Also, determine the indexes where the new 1-min files will start. 
    start_at_zero = false; 
    for i=1:combined_dataSize
        
        dateTime = datetime(unix_Data(i),'ConvertFrom', 'posixtime');
                
        year = sprintf('%02d', dateTime.Year);
        julianDayNumber = num2str(floor(datenum(dateTime) - datenum(dateTime.Year, 1, 1) + 1));
        hour = sprintf('%02d', dateTime.Hour);
        min = sprintf('%02d', dateTime.Minute);
        sec = sprintf('%015.12f',dateTime.Second);
        special_format = str2double(append(year(end-1:end),julianDayNumber,hour,min,sec));
        special_Data(i) = special_format;
        check_sec = floor(str2double(sec));
        
        % Determines the index to start at for each file to begin at a
        % new minute exactly.
        if start_at_zero == false &&  check_sec == 0
            starting_indexes(k) = i;
            k=k+1;
            start_at_zero = true; 
        end
        if start_at_zero == true && check_sec == 2
            start_at_zero = false;
        end
    end 

    % Divide the data into 1-min .mat files. 
    for x=1:num_min_files
        begin = starting_indexes(x);
        finish = starting_indexes(x+1)-1;
        array_size = finish - begin + 1;
        
        t_unixsec = zeros(1, array_size);
        t_special = zeros(1, array_size);
        hyd_data = zeros(channels, array_size);

        for y=begin:finish 
            t_unixsec(y-begin+1)=unix_Data(y);
            t_special(y-begin+1)=special_Data(y);
            hyd_data(:,y-begin+1)=hyd_Data(:,y);
        end 

        file_time_stamp = num2str(round(t_special(1)));
        
        save_name = append(save_path,file_prefix,file_time_stamp,'.mat');
        save(save_name,"t_unixsec","t_special","hyd_data")
    end

end 

%% This case is for when there are more than 10 files in the folder. 
%  The code will process all files in multiples of 10. 
if repeat_process==true
    for j = 1:num_groups
        k=1; 

        if j == 1
            % Preallocate the data arrays
            num_files = 10;
            combined_dataSize = size_5min_data * num_files;
            hyd_Data = zeros(channels, combined_dataSize);
            unix_Data = zeros(1, combined_dataSize);
            special_Data = zeros(1, combined_dataSize);
            num_min_files = floor(combined_dataSize / size_1min_data);
            starting_indexes = zeros(1,num_min_files);

            % Load and combine data from multiple files
            for i = 1:10
                % Load your data from each file
                file = MyFolderInfo(i).name;
                data = load(file,'hyd_data', 't_unixsec');
                % Add the data to the preallocated arrays
                hyd_Data(:, (i-1) * size_5min_data + 1 : i * size_5min_data) = data.hyd_data;
                unix_Data((i-1) * size_5min_data + 1 : i * size_5min_data) = data.t_unixsec;
            end
    
        else 
            % Preallocate the data arrays
            num_files = 11;
            combined_dataSize = size_5min_data * num_files;
            hyd_Data = zeros(channels, combined_dataSize);
            unix_Data = zeros(1, combined_dataSize);
            special_Data = zeros(1, combined_dataSize);
            num_min_files = floor(combined_dataSize / size_1min_data);
            starting_indexes = zeros(1,num_min_files);

            % Load and combine data from multiple files
            file_start = 10*(j-1);
            file_end = 10*j;
            y=1;
            for i = file_start:file_end
                % Load your data from each file
                file = MyFolderInfo(i).name;
                data = load(file,'hyd_data', 't_unixsec');
                % Add the data to the preallocated arrays
                hyd_Data(:, (y-1) * size_5min_data + 1 : y * size_5min_data) = data.hyd_data;
                unix_Data((y-1) * size_5min_data + 1 : y * size_5min_data) = data.t_unixsec;
                y=y+1;
            end
        end 
        
        % Allocate the special_Data array with the time converted from unix
        % time to special time (YYDDDHHMMSS, with DDD being the Julian day). 
        % Also, determine the indexes where the new 1-min files will start. 
        start_at_zero = false; 
        
        for i=1:combined_dataSize
                      
            dateTime = datetime(unix_Data(i),'ConvertFrom', 'posixtime');
                    
            year = sprintf('%02d', dateTime.Year);
            julianDayNumber = num2str(floor(datenum(dateTime) - datenum(dateTime.Year, 1, 1) + 1));
            hour = sprintf('%02d', dateTime.Hour);
            min = sprintf('%02d', dateTime.Minute);
            sec = sprintf('%015.12f',dateTime.Second);
            special_format = str2double(append(year(end-1:end),julianDayNumber,hour,min,sec));
            special_Data(i) = special_format;
            check_sec = floor(str2double(sec));
            
            % Determines the index to start at for each file to begin at a
            % new minute exactly.
            if start_at_zero == false &&  check_sec == 0
                starting_indexes(k) = i;
                k=k+1;
                start_at_zero = true; 
            end
            if start_at_zero == true && check_sec == 2
                start_at_zero = false; % resets the boolean
            end
        end 
    
        % Divide the data into 1-min .mat files.      
        for x=1:(num_min_files-1)
            begin = starting_indexes(x);
            finish = starting_indexes(x+1)-1;
            array_size = finish - begin + 1;
            
            t_unixsec = zeros(1, array_size);
            t_special = zeros(1, array_size);
            hyd_data = zeros(channels, array_size);
    
            for y=begin:finish 
                t_unixsec(y-begin+1)=unix_Data(y);
                t_special(y-begin+1)=special_Data(y);
                hyd_data(:,y-begin+1)=hyd_Data(:,y);
            end 
    
            file_time_stamp = num2str(round(t_special(1)));
            
            save_name = append(save_path,file_prefix,file_time_stamp,'.mat');
            save(save_name,"t_unixsec","t_special","hyd_data")
        end
    end 

end 

