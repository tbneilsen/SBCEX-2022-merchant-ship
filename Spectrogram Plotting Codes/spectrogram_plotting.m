%% Spectrogram Plotting 

%   This code creates spectograms for the VLA data from Scripps saved in
%   1-min .npy files and the PROTEUS data from ARL:UT saved in 1-min .mat files.

%   The code is divided into sections, where the first section is the only
%   section that the user needs to edit based on their specific plotting
%   parameters (User Defined Parameters). 

%   This code reads the interested ships spreadsheets and gathers
%   information for a specified ship. Using the timestep at CPA from the
%   spreadsheet, data is loaded and then plots are made and saved as .png files
%   for each specified freqeuncy range and channel by the user. 

%   As an option for the user, a specific frequency at beta = infinity may
%   be specified and the a white dashed line at that frequency will also
%   appear on the spectrogram. 

%   For example: 
%   The user could make the plots for the ship 'CARMEN' going from 10-90 Hz
%   and 90-1100 Hz on all 16 channels of VLA1 with a single run of this script. 
% 
%   Author: Alexandra Hopps McDaniel
%   Created: 9/22/23
%   Modified: 10/4/23

%% Start up
clear all;
%close all;

%% User Defined Parameters

% Add paths to the folders that contain the time series data for PROTEUS and VLA1 and VLA2 at each location
addpath(genpath('F:\VLA1-2 (Scripps)\data_files\1-min_npy')) % Path to VLA1 and VLA2 data 
addpath(genpath('G:\PROTEUS\2022 Data\1-min Files')) % Path to PROTEUS data 
addpath(genpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Notes\SBCEX Spreadsheets')) % Path to the ship spreadsheets
addpath(genpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\sbc-vla-data\prossessing_MPL_data'))
% For processing .npy files (Scripps data), the npy-matlab library must be downloaded or cloned to your machine. 
% You can find the npy-matlab library on GitHub: https://github.com/kwikteam/npy-matlab. 
% Add path to the folder that contains the script 'readNPY.m'
addpath(genpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Notes\SBCEX Spreadsheets'));

% Input which array's data you will be plotting
which_array = 'VLA2'; % Can be 'PROTEUS', 'VLA1', OR 'VLA2'

% Include the file name for the spreadsheet with ships that approach 15km or less to the arrays. 
% The spreadsheet includes information regarding time stamp, distance, and speed for CPA.
% For 2022, spreadsheet can be: 'SHIPS_SBCEX2022_VLA1_VLA2.xlsx' or 'SHIPS_SBCEX2022_PROTEUS.xlsx'
ships_spreadsheet = 'SHIPS_SBCEX2022_VLA1_VLA2.xlsx';

% Sheet name can be: 'Proteus', 'VLA1_Location1', 'VLA1_Location2', 'VLA2_Location1', or 'VLA2_Location2',
sheet_name = 'VLA2_Location1'; 

% Index of the interested ship will be equal to (row number - 1) on the spreadsheet. 
interested_ship = 16;

% Array channels to read (any of [1:52] for PROTEUS_22 and [1:16] for VLA1 and VLA2)
channels = 2;

% Add information requarding the frequency range and frequency resolution desired for each plot. 
% Information is ordered as follows: [fmin1, fmax1, nfft1; fmin2, fmax2, nfft2; ...] 
% Note: A larger nfft provides better frequency resolution
% Suggested: [20,40,2^18; 20,100,2^18] (2^18 is a good nfft for VLA 1 & 2; 2^17 is a good nfft for PROTEUS)
spec_info = [15,80,2^17];%20,40,2^19];%[15,80,2^19;20,40,924288]; %[15,80,2^18;20,40,2^18];% 20,100,2^18]; 

% Interval of time that the spectrogram should display in minutes. 
TimeInterval = 20;

% Specify whether or not you want to save a .mat file with the spectrogram info 
save_mat = false;
save_mat_path = 'C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\matlab-bborca\Inversion\Inversion_4000_iteration_scans\Ship Specs .mat files\';

save_spec = true;
% Specify whether or not you want to make a plot with the beta line visible
show_beta = true; % Can be true or false 
beta_val = [27.8]; % Frequency where beta goes to +/- infinity 

% Specify whether or not you want the plot to have a title
has_title = false; 

% The path to the folder where the spectrograms .png files should be saved. 
% Matlab does not create folders when saving, so make sure the correct folders are already created. 
% Currently, this appends the following to the folder paths depending on the condition of show_beta as true or false. 
% If show_beta == false: 
% '...\{ship_name}\{which_array}\{fmin}_{fmax}Hz_{TimeInterval}m\'.
% Example: '...\CARMEN\PROTEUS\10_90Hz_25m\'.
% If show_beta == true:
% '...\FORMAL PLOTS\{fmin}_{fmax}Hz_{TimeInterval}m\'];
% Example: '...\FORMAL PLOTS\10_90Hz_25m\'.
save_folder_path = 'C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\matlab-bborca\Inversion\Inversion_4000_iteration_scans\actual_specs\';
%% Read ship info

warning('off', 'MATLAB:table:ModifiedVarnames'); % Turns off a warning about reading in the ship names off of the spreadsheet. 

% Reads the information on the interested ships spreadsheet.
ship_info = readtable(ships_spreadsheet, 'Sheet', sheet_name);

% Gives the information (name and timestamp, distance, and speed in knots at cpa) regarding the ship.
ship_name = cell2mat(ship_info.Vessel(interested_ship));
%ship_cpa = sprintf('%.1f', ship_info.CPA(interested_ship));
%ship_sog = sprintf('%.1f', ship_info.SOG_knots_(interested_ship));

%% Load Data

npy_file = false; % Boolean determining whether or not .npy files will be processed

% Determines the sampling frequency in Hz (Fs) and prefix and suffix for data files based on the array. 
% Also reads the timestamp from the ship spreadsheet, which reads differently on each spreadsheet. 
condition1 = strcmp(which_array, 'PROTEUS');
condition2 = strcmp(which_array, 'VLA1');
condition3 = strcmp(which_array, 'VLA2');
if condition1
    Fs = 10240;
    prefix = 'PROTEUS.';
    suffix = '00.mat';
    timeCell = ship_info.TimeStamp(interested_ship);
    time = timeCell{1};
elseif condition2
    Fs = 25000;
    prefix = 'RAVA01.';
    suffix = '00.npy';
    time = sprintf('%f', ship_info.TimeStamp(interested_ship));
    npy_file = true;
elseif condition3
    Fs = 25000;
    prefix = 'RAVA02.';
    suffix = '00.npy';
    time = sprintf('%f', ship_info.TimeStamp(interested_ship));
    npy_file = true;
else
    disp('Wrong array defined by the user.')
end

% This is the time before and after CPA time that we want to load. 
time_change = round(TimeInterval / 2);

% Separates the time stamp into the year, day, hour, min, and sec components and saves as doubles. 
year = str2double(time(1:2));
days = str2double(time(3:5));
hours = str2double(time(6:7));
mins = str2double(time(8:9));
secs = str2double(time(10:11));

% Rounds to the nearest minute. 
if secs >= 30
    secs = 0;
    mins = mins + 1;
elseif secs < 30
    secs = 0;
end

% Edits the timestamp to be start time of the spectogram (which is time - time_change).
if mins - time_change < 0
    mins = mins + 60 - time_change;
    hours = hours - 1;
elseif mins - time_change >= 0
    mins = mins - time_change;
end

if hours < 0
    days = days - 1;
    hours = 24 + hours;
end

% Finds the name of all data files required for the time interval. 
for x = 1:TimeInterval
    if mins + x >= 60
        file_mins = mod((mins + x), 60);
        file_hours = hours + 1;
    elseif mins + x < 60
        file_mins = mins + x;
        file_hours = hours;
    end 
    
    if file_hours >= 24
        file_days = days + 1;
        file_hours = mod(file_hours, 24);
    else
        file_days = days;
    end
    
    % Formats the timestamp for each data file.
    file_time = sprintf('%02d%03d%02d%02d',year,file_days,file_hours,file_mins);
    filename=[prefix,file_time,suffix] %#ok<NOPTS>

    % Reading in the data is different for .mat and .npy files. 
    % Only the specified channels are saved. 
    if npy_file == false
        matObj = matfile(filename);
        var_name = 'hyd_data';
        data = matObj.(var_name)(channels,:);
    else
        data = readNPY(filename);
        data = data';
        data = data(channels,:);
    end 

    % Takes the data array and appends it to dat each iteration. 
    if x==1
        dat=data;
    else
        dat = [dat,data];
    end
end

dat = double(dat); % The values in dat must be doubles for the spectrogram to compute correctly. 

% The data contained in the .npy files for VLA 1 & 2 was not previously calibrated. 
% The following code correctly calibrates the code before performing the fft. 
if condition2 || condition3
    calibration = 148;
    voltRange = 2.5 * 2;
    ADrange = (2^16) - 1;
    calfactor = 10^(calibration/20) * voltRange/ADrange;
    dat = dat * calfactor;
end

%% Plot Spectrograms

num_plots = size(spec_info,1); % The number of plots that will saved for each channel of the specified ship. 
numchans=length(channels); % Number of channels for each freqeuncy range. 

for x=1:num_plots

    % Specify the fmin, fmax, and nfft for each plot
    fmin = spec_info(x,1);
    fmax = spec_info(x,2);
    nfft = spec_info(x,3);

    window = hamming(nfft); % Use a hamming window on the data 
    noverlap = nfft/2; % nfft/2 is 50% overlap

    % Define frequency array to match the output of the fft
    frequ_array = (0:nfft/2)*Fs/nfft; % Full frequency axis in Hz
    frequ_index = find(frequ_array >= fmin & frequ_array <= fmax); % Finds the indexes within the full frequency axis for the specific plot 
   
    for i=1:length(channels)
        % Compute the spectrogram
        [sspec,fspe,tspec,psspec] = spectrogram(dat(i,:), window, noverlap, nfft, Fs,'PSD');
                
        plot_PSD = 10*log10(psspec(frequ_index,:)); % Compute spectral levels
        plot_freq = frequ_array(frequ_index); % The frequency array to be plotted
        
        % Settings for the colorbar
        clims(2)=max(plot_PSD(:));
        clims(1)=clims(2)-60;

        chan_num = num2str(channels(i)); % Current channel number. 
    
        spec_Title = [ship_name,', ', which_array, ', CH:', chan_num] ;
        fig = figure('name', spec_Title);
        clims(2)    = 130;
        clims(1)    = clims(2)-60;
        
        % Plots the spectrogram
        imagesc(plot_freq, tspec, plot_PSD', clims)
        % Figure parameters 
        colormap parula
        c = colorbar;
        
        xlabel('Frequency (Hz)', 'FontSize', 15)
        ylabel('Time (sec)', 'FontSize', 15)
        % Font sizes
        c.FontSize = 14;
        c.Label.String = 'LEVEL: dB re: 1 \muPa^2/Hz';
        c.Label.FontSize = 16;
        
        ax = gca; 
        ax.XAxis.FontSize = 15; 
        ax.YAxis.FontSize = 15; 
        

        png_name = [ship_name, '_', which_array, '_CH', chan_num,'BIGGER3.png'];
        
        if has_title == true
            title(spec_Title)
        end

        if show_beta == true
            % Plots the dashed line at the beta equals infinity frequency 
            ylim = get(gca, 'YLim');
            
            
            for i=1:length(beta_val)
                line([beta_val(i) beta_val(i)], ylim, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 1.5);
                % Labels the dashed line with the frequency
                [~, text_index] = min(abs(plot_freq - beta_val(i))); % Find the index in the x-axis (plot_freq) closest to the target frequency
                y_coordinate = -60;  % Adjust the vertical position as needed
                text_str = ['$f_{1,2}$ = ', sprintf('%#.1f', beta_val(i))];
                if i==2
                    %text_index=text_index+35;
                end
                text(plot_freq(text_index), y_coordinate, text_str, 'Interpreter','latex','HorizontalAlignment', 'center', 'FontSize', 18); % Add the text above the spectrogram at the specified frequency
            end
            %save_path = [save_folder_path,ship_name,'\FORMAL PLOTS\',num2str(fmin),'_',num2str(fmax),'Hz_',num2str(TimeInterval),'m\',png_name];
       
        else
            %save_path = [save_folder_path,ship_name,'\',which_array,'\',num2str(fmin),'_',num2str(fmax),'Hz_',num2str(TimeInterval),'m\',png_name];
        
        end 
        save_nam=[save_folder_path,png_name];
        if save_spec==true
            saveas(fig, save_nam)
        end 
        

    end

end

%% Save .mat file of spectrogram 
if save_mat==true
    fselect=plot_freq';
    PSDselectf=plot_PSD;
    save_mat=[save_mat_path,ship_name,'_',which_array,'.mat'];
    save(save_mat,'PSDselectf','fselect','tspec')
end 
% [15,80,2^19] This is the spec infor needed




