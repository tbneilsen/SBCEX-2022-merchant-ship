%% Start
% Add paths to the folders that contain the time series data for PROTEUS and VLA1 and VLA2 at each location
addpath(genpath('F:\VLA1-2 (Scripps)\data_files\1-min_npy')) % Path to VLA1 and VLA2 data 
addpath(genpath('G:\PROTEUS\2022 Data\1-min Files')) % Path to PROTEUS data 
addpath(genpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Notes\')) % Path to the ship spreadsheets
addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\matlab-bborca\Inversion\Inversion_4000_iteration_scans\dist_files\')
addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\matlab-bborca')
addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\sbc-vla-data\prossessing_MPL_data\')
nmodes          = 7;
num_arrays = 3;
TimeInterval = 20;
fmin = 15;
fmax = 800;
ssp_file         = 'SBCEX22_SSP';     % File name with 20 water column SSPs
% Calibration for the VLA 1&2 files
calibration = 148;
voltRange = 2.5 * 2;
ADrange = (2^16) - 1;
calfactor = 10^(calibration/20) * voltRange/ADrange;
% Spread sheet with all interested ships 
ships_spreadsheet = 'SBCEX_ALL_INTERESTED_SHIPS_FINAL.xlsx';
% Sheet name 
sheet_name = 'ALL SHIP INFO NEW'; 
all_ship_info = readtable(ships_spreadsheet, 'Sheet', sheet_name);
[num_ship, ~] = size(all_ship_info);
max_clims=zeros(num_ship,3);
CH=zeros(num_arrays,1);
Time=zeros(num_arrays,1);
CPA=zeros(num_arrays,1);
SOG=zeros(num_arrays,1);
b12=zeros(num_arrays,1);
b13=zeros(num_arrays,1);
prefix={'RAVA01.';'PROTEUS.';'RAVA02.'};
suffix={'00.npy';'00.mat';'00.npy'};
fs=[25000;10240;25000];
nfft=[2^18;2^17;2^18];
array={'VLA1';'PROTEUS';'VLA2'};
type=[0;1;0];

plotting_info = table(CH, Time, CPA, SOG, b12, b13, prefix, suffix, fs, nfft, array, type);


zmax                = 950;                          % maximum depth point for pressure field
delz                = 1;                            % depth mesh size  [m]
z                   = delz:delz:zmax;               % depth vector [m]

Basespeed           = 2350;                                         % sound speed in the basement 

svp_in.uphalf_cp    = 343.0;                                        % compressional sound speed in upper half space [m/s]
svp_in.uphalf_cs    = 0.0;                                          % shear sound speed in upper half space [m/s]
svp_in.uphalf_rho   = .00121;                                       % density in upper half space [g/cc]
svp_in.uphalf_ap    = 0.000058;                                     % compressional attenuation in upper half space [pos~db/m/kHz, neg~dB/wavelength]
svp_in.uphalf_as    = 0.0;                                          % shear attenuation in upper half space [dB/m/kHz]

svp_in.ctol         = 1.0;                                          % tolerance used in fitting SVP to eliminate layers (0 = keep all layers)


% svp_in.btm_env      = sedtemp;                                    % Commented out because this paramter will be updated during the Monte Carlo iterations 
svp_in.lowhalf_cp   = Basespeed;                                    % compressional sound speed in lower halfpsace [m/s]
svp_in.lowhalf_cs   = 0;                                            % shear sound speed in lower halfspace [m/s]
svp_in.lowhalf_rho  = 2.6;                                          % density in lower halfspace [g/cc]
svp_in.lowhalf_ap   = 0.22;                                       % compressional attenuation in lower halfspace [pos~db/m/kHz, neg~dB/wavelength]
svp_in.lowhalf_as   = 0.0;                                          % shear attenuation in lower halfspace [pos~db/m/kHz, neg~dB/wavelength]
svp_in.ntop         = 0;                                            % number of surface layers to be read in on next line
svp_in.above_sea    = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

%% orca90 .opt file -------------------------------------------------------
% Create the equivalent of ORCA.opt file
opt_in.nmode        = nmodes;                                       % number of modes
opt_in.cphmin       = 0;                                            % min phase speed (0=p-wave modes only; -1=seismic modes also);
opt_in.cphmax       = 0;                                            % max phase speed [0=use rmin; pos=speed; neg=max angle(deg)];
opt_in.rmax         = 100   ;                                       % max range of interest in km (0=use S/R geom);
opt_in.phfac        = 4;                                            % phase step parm: step by 2*pi/phfac (set to 4-8, 0=default==>4);
opt_in.dbcut        = 100;                                          % modes weaker by db_cut ignored (set to 30-60, 0=default==>50);
opt_in.Aih_l        = 135.001;                                      % gradient lower h-space: 0=default,-1=homogeneous, >0=da_bar;
opt_in.Aih_u        = 0;                                            % gradient upper h-space: 0=default,-1=homogeneous, >0=da_bar;
opt_in.nf           = 1;                                            % frequencies (nf>0 ==> List fcw's; nf<0 ==> List first,last f)

opt_in.nzm          = -length(z);                                   % depth points in modes (nzm>0 ==> List zm's; nf<0 ==> List first,last zm)
opt_in.zm           = [z(1) z(end)];                                % [z1 z2 ... ]
opt_in.zm_n         = length(opt_in.zm);                            % size of the zm array

iimf = 1;            

%%

for i=15
    ship = all_ship_info.Var1{i};
    plotting_info.CH = [all_ship_info.CH(i); all_ship_info.CH_1(i); all_ship_info.CH_2(i)];
    plotting_info.Time = [all_ship_info.TIME(i); all_ship_info.TIME_1(i); all_ship_info.TIME_2(i)];
    plotting_info.CPA = [all_ship_info.CPA(i); all_ship_info.CPA_1(i); all_ship_info.CPA_2(i)];
    plotting_info.SOG = [all_ship_info.SOG(i); all_ship_info.SOG_1(i); all_ship_info.SOG_2(i)];
    plotting_info.b12 = [all_ship_info.b12(i); all_ship_info.b12_1(i); all_ship_info.b12_2(i)];
    plotting_info.b13= [all_ship_info.b13(i); all_ship_info.b13_1(i); all_ship_info.b13_2(i)];
    ssp_num = all_ship_info.Var4(i);
    
        
    %ssp = load(ssp_file);                % Load water Sound Speed profile given by name wcpname
              
    %c        = ssp.SSP{ssp_num};         % m/s water SS array. rows: depth veriation; Columns: different SS profiles
    %c_z      = ssp.C_Z{ssp_num};             % (m) depth variation
    %rho0     = 1.05;      % water density [g/cc]
    %zbottom         = c_z(end);                                         % bottom depth
    %z_bathy         = 1:delz:zbottom;
    
    
    %cofztemp(:,1)   = interp1(c_z, c_z, z_bathy,'linear','extrap').';   % extrapolated depth points for wcssp [m]
    %cofztemp(:,2)   = interp1(c_z, c, z_bathy,'linear','extrap').';     % extrapolated soundspeed for wcssp [m/S]
    %cofz(:,1)       = linspace(cofztemp(1,1),cofztemp(end,1),250);      % interpolate new depth of water column profile
    %cofz(:,2)       = interp1(cofztemp(:,1),cofztemp(:,2),cofz(:,1));   % interpolate ssp at new depth of water column
    
    %Bwss                = cofz(end,2);                                  % Bwss=bottom water sound speed
    %svp_in.nsvp         = length(cofz);                                 % number of SVP points in the ocean to be read
    %svp_in.wssp         = cofz;                                         % SVP profile
    %svp_in.wrho         = rho0;                                         % water density [g/cc]                                      % number of sediment layers to be read in
    
    for k=3
        no_info = isnan(plotting_info.CH(k));
        if no_info == true
            continue; % skips to the next iteration
        end 

        % plot actual data
        time = num2str(plotting_info.Time(k));
        % This is the time before and after CPA time that we want to load. 
        time_change = round(TimeInterval / 2);
        
        % Separates the time stamp into the year, day, hour, min, and sec components and saves as doubles. 
        year = str2double(time(1:2));
        days = str2double(time(3:5));
        hours = str2double(time(6:7));
        mins = str2double(time(8:9));
        secs = str2double(time(10:11));
        

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
            clear min;
            % Formats the timestamp for each data file.
            file_time = sprintf('%02d%03d%02d%02d',year,file_days,file_hours,file_mins);
            
            filename=[cell2mat(plotting_info.prefix(k)),file_time,cell2mat(plotting_info.suffix(k))] %#ok<NOPTS>
            channel=plotting_info.CH(k);
            % Reading in the data is different for .mat and .npy files. 
            % Only the specified channels are saved. 
            if plotting_info.type(k) == 1
                matObj = matfile(filename);
                var_name = 'hyd_data';
                data = matObj.(var_name)(channel,:);
            else
                data = readNPY(filename);
                data = data';
                data = data(channel,:);
            end 

            % Takes the data array and appends it to dat each iteration. 
            if x==1
                dat=data;
            else
                dat = [dat,data];
            end
        end

        dat = double(dat); % The values in dat must be doubles for the spectrogram to compute correctly. 

        if plotting_info.type(k) == 0
            dat = dat * calfactor;
        end 


        nfft = plotting_info.nfft(k);
        fs = plotting_info.fs(k);
        b12 = plotting_info.b12(k);
        b13=plotting_info.b13(k);
        window = hamming(nfft); % Use a hamming window on the data 
        noverlap = nfft/2; % nfft/2 is 50% overlap

        % Define frequency array to match the output of the fft
        frequ_array = (0:nfft/2)*fs/nfft; % Full frequency axis in Hz
        frequ_index = find(frequ_array >= fmin & frequ_array <= fmax); % Finds the indexes within the full frequency axis for the specific plot 
           
        [sspec,fspe,tspec,psspec] = spectrogram(dat, window, noverlap, nfft, fs,'PSD');
        plot_PSD = 10*log10(psspec(frequ_index,:)); % Compute spectral levels
        plot_freq = frequ_array(frequ_index); % The frequency array to be plotted\
        max_psd = max(plot_PSD(:));
        if max_psd <= 100 
            clims(2) = 100;
        elseif max_psd <= 105 
            clims(2) = 105;
        elseif max_psd <= 110 
            clims(2) = 110;
        elseif max_psd <= 115 
            clims(2) = 115;
        elseif max_psd <= 120 
            clims(2) = 120;
        elseif max_psd <= 125 
            clims(2) = 125;
        else
            clims(2) = 130;
        end 
        max_clims(i,k)=clims(2);
        clims(1)=clims(2)-60;
        spec_title = [ship,'_',plotting_info.array{k},'_full_spec_CH',num2str(channel)];
        % Plots the spectrogram
        fig = figure('name', spec_title);
        imagesc(plot_freq, tspec, plot_PSD', clims)
        % Figure parameters 
        colormap parula
        % Settings for the colorbar
        c = colorbar;
        c.Label.String = 'LEVEL: dB re: 1 \muPa^2/Hz';
        xlabel('Frequency (Hz)', 'FontSize', 16)
        ylabel('Time (sec)', 'FontSize', 16)
        
        % Font sizes
        c.Label.FontSize = 16;
        c.FontSize = 12;
        ax = gca; 
        ax.XAxis.FontSize = 12; 
        ax.YAxis.FontSize = 12; 
        ylim = get(gca, 'YLim');
        %line([b12 b12], ylim, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 2);
        % Labels the dashed line with the frequency
        %[~, text_index] = min(abs(plot_freq - b12)); % Find the index in the x-axis (plot_freq) closest to the target frequency
        %y_coordinate = -60;  % Adjust the vertical position as needed
        %text_str = ['$f_{1,2}$ = ', sprintf('%#.1f', b12)];
        %text(plot_freq(text_index), y_coordinate, text_str, 'Interpreter','latex','HorizontalAlignment', 'center', 'FontSize', 18); % Add the text above the spectrogram at the specified frequency
        
        
        %line([b13 b13], ylim, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 2);
        % Labels the dashed line with the frequency
        %[~, text_index] = min(abs(plot_freq - b13));
        %text_str = ['$f_{1,4}$ = ', sprintf('%#.1f', b13)];
        %text(plot_freq(text_index), y_coordinate, text_str, 'Interpreter','latex','HorizontalAlignment', 'center', 'FontSize', 18);
        save_name=['C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\matlab-bborca\Inversion\Inversion_4000_iteration_scans\New actual specs\',spec_title,'_SPEC.png'];
        saveas(fig, save_name)
        % Plot the model 
        fselect=plot_freq';
        PSDselectf=plot_PSD;
        
        save_mat=['C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\matlab-bborca\Inversion\Inversion_4000_iteration_scans\New .mat files\',spec_title,'.mat'];
        %save(save_mat,'PSDselectf','fselect','tspec')
        %% insert cut lines
    end 
end 

%parms = info.optimal_parm;








                




