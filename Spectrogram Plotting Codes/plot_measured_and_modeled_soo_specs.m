%% 
%% Spectrogram Plotting 

%   This code can mass creates spectograms for the VLA data from Scripps saved in
%   1-min .npy files and the PROTEUS data from ARL:UT saved in 1-min .mat files.
%   It will also create the related model spectrograms using the lowest
%   cost seabed parameters from the maximum entropy inversions. 

%   This code reads the interested ships spreadsheets and gathers
%   information for a specified ship. Using the timestep at CPA from the
%   spreadsheet, data is loaded and then plots are made and saved as .png files
%   for each specified freqeuncy range and channel by the user. 

%   It reads the optimal parameter values from the distribution files saved
%   as .mat files.

%   As an option for the user, a specific frequency at beta = infinity may
%   be specified and the a white dashed line at that frequency will also
%   appear on the spectrogram. 

%   Author: Alexandra Hopps McDaniel
%   Created: 9/22/2023
%   Modified: 5/29/2024

%%
clear all; close all;
currentFilePath = mfilename('fullpath');

%% 
% Add paths to the folders that contain the time series data for PROTEUS and VLA1 and VLA2 at each location
addpath(genpath('F:\VLA1-2 (Scripps)\data_files\1-min_npy')) % Path to VLA1 and VLA2 data 
addpath(genpath('G:\PROTEUS\2022 Data\1-min Files')) % Path to PROTEUS data 

% Add path to matlab version of ORCA
addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\matlab-bborca')


%% The script accesses the correct folder paths

[currentFolderPath, ~, ~] = fileparts(currentFilePath);
git_library = fileparts(currentFolderPath);

% Add path to the ship spreadsheets
spreadsheet_path = fullfile(git_library,'Merchant Ship Spreadsheets');
addpath(spreadsheet_path) 

% Path to scripts for processing .npy files (Scripps data)
array_processing_path = fullfile(git_library,'Processing array file types');
addpath(array_processing_path)

% Path to scripts for processing .npy files (Scripps data)
dist_files_path = fullfile(git_library,'SOO dist files');
addpath(dist_files_path)

% Path to scripts for processing .npy files (Scripps data)
soo_mat_path = fullfile(git_library,'SOO .mat files');
addpath(soo_mat_path)


% Path to scripts for maximum entropy inversion codes
inversion_path = fullfile(git_library,'Inversion Codes');
addpath(inversion_path)

%% Add information for the spectrograms
nmodes              = 7;
num_arrays          = 3;
TimeInterval        = 20;
fmin                = 15;
fmax                = 80;
ssp_file            = 'SBCEX22_SSP';     % File name with 20 water column SSPs
plot_betas = true;
plot_gv = true;
% Calibration for the VLA 1&2 files
calibration = 148;
voltRange = 2.5 * 2;
ADrange = (2^16) - 1;
calfactor = 10^(calibration/20) * voltRange/ADrange;

% Spread sheet with all interested ships 
ships_spreadsheet = 'SBCEX_ALL_INTERESTED_SHIPS_FINAL.xlsx';

% Sheet name 
sheet_name      = 'ALL SHIP INFO NEW'; 
all_ship_info   = readtable(ships_spreadsheet, 'Sheet', sheet_name);
[num_ship, ~]   = size(all_ship_info);
max_clims       = zeros(num_arrays,1);
CH              = zeros(num_arrays,1);
Time            = zeros(num_arrays,1);
CPA             = zeros(num_arrays,1);
nparm           = zeros(num_arrays,1);
SOG             = zeros(num_arrays,1);
b12             = zeros(num_arrays,1);
b13             = zeros(num_arrays,1);
prefix          = {'RAVA01.';'PROTEUS.';'RAVA02.'};
suffix          = {'00.npy';'00.mat';'00.npy'};
fs              = [25000;10240;25000];
nfft            = [2^18;2^17;2^18];
array           = {'VLA1';'PROTEUS';'VLA2'};
varNames        = {'Ship', 'Array', 'B12', 'B13', 'B14'};
varType         = {'string', 'string', 'double', 'double', 'double'};
optimalInfo     = table('Size', [0 length(varNames)], 'VariableTypes', varType, 'VariableNames', varNames);
type            = [0;1;0];

plotting_info   = table(CH, Time, CPA, SOG, b12, b13, prefix, suffix, fs, nfft, array, type, max_clims);


zmax                = 950;                          % maximum depth point for pressure field
delz                = 1;                            % depth mesh size  [m]
z                   = delz:delz:zmax;               % depth vector [m]


%%

for i=1:num_ship % can be i=1:num_ship or i=#
    ship = all_ship_info.Var1{i};
    plotting_info.CH = [all_ship_info.CH(i); all_ship_info.CH_1(i); all_ship_info.CH_2(i)];
    plotting_info.Time = [all_ship_info.TIME(i); all_ship_info.TIME_1(i); all_ship_info.TIME_2(i)];
    plotting_info.CPA = [all_ship_info.CPA(i); all_ship_info.CPA_1(i); all_ship_info.CPA_2(i)];
    plotting_info.SOG = [all_ship_info.SOG(i); all_ship_info.SOG_1(i); all_ship_info.SOG_2(i)];
    plotting_info.b12 = [all_ship_info.b12(i); all_ship_info.b12_1(i); all_ship_info.b12_2(i)];
    plotting_info.b13= [all_ship_info.b13(i); all_ship_info.b13_1(i); all_ship_info.b13_2(i)];
    plotting_info.max_clims = [all_ship_info.MaxDB(i); all_ship_info.MaxDB_1(i); all_ship_info.MaxDB_2(i)];
    ssp_num = all_ship_info.Var4(i);
    
    for k=1:3
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
        nfft        = plotting_info.nfft(k);
        fs          = plotting_info.fs(k);
        b12         = plotting_info.b12(k);
        b13         = plotting_info.b13(k);
        window      = hamming(nfft); % Use a hamming window on the data 
        noverlap    = nfft/2; % nfft/2 is 50% overlap

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

        colorbar
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
        
        ylim = get(gca, 'YLim');
        
        if plot_betas==true
            line([b12 b12], ylim, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 2);
            %Labels the dashed line with the frequency
            [~, text_index] = min(abs(plot_freq - b12)); % Find the index in the x-axis (plot_freq) closest to the target frequency
            y_coordinate = -60;  % Adjust the vertical position as needed
            text_str = ['$f_{1,2}$ = ', sprintf('%#.1f', b12)];
            text(plot_freq(text_index), y_coordinate, text_str, 'Interpreter','latex','HorizontalAlignment', 'center', 'FontSize', 18); % Add the text above the spectrogram at the specified frequency
            
            %line([b13 b13], ylim, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 2);
            % Labels the dashed line with the frequency
            %[~, text_index] = min(abs(plot_freq - b13));
            %text_str = ['$f_{1,4}$ = ', sprintf('%#.1f', b13)];
            %text(plot_freq(text_index), y_coordinate, text_str, 'Interpreter','latex','HorizontalAlignment', 'center', 'FontSize', 18);
        end
        save_name=[spec_title,'.png'];
        save_name=fullfile(git_library,'Spectrograms/Measured/',save_name);
        saveas(fig, save_name)

        %% Plot the model 

        array = plotting_info.array{k};
        channel = plotting_info.CH(k);
        CPA = plotting_info.CPA(k);
        sourceKnots = plotting_info.SOG(k);
        clim2 = plotting_info.max_clims(k);

        % load in the .mat file 
        file_name           = [ship, '_', array,'.mat'];
        load(file_name, 'tspec','fselect','PSDselectf')
        file_name           = [ship, '_', array, '_dist_file_4000_iterations_01.mat'];
        load(file_name, 'info')
        params              = info.optimal_parm;
        b12                 = info.optimal_betas;
        aa                  = 3.75;                         % Distance between channels in VLA
        dist_channel1       = 7.25;                         % Distance from ocean floor to the first hydrophone 
        
        water_depth         = 75;
        zs                  = 9.2;                          % Source position (m) (note: the source comes from the propellers of the ship which depth can be estimated as the kiln depth found online)
        
        zr                  = water_depth - dist_channel1 - channel*aa;             % Reciever position (m)
        nTT                 = length(tspec);                % number times for calculating each spectrogram

        freq                = fselect;
        fmin                = fselect(1);                   % minimum frequency (Hz) 
        fmax                = fselect(length(fselect));     % maximum frequency (Hz) of analysis
        df                  = (fselect(2)-fselect(1)) ;     % frequncy resolution 
        nfreq               = length(freq);                 % number of Orca frequencies evaluated
        
        S_0                 = 237.0;                        % y-intercept of Wales-Heitmyer container ship source levels broadband
        [PSource,~]         = s00spectrum(freq,S_0);        % PSource(freq) mean complex pressure power density (uPa/Hz) @ 1 m w/ random phase
        
        zmax                = 950;                          % maximum depth point for pressure field
        delz                = 1;                            % depth mesh size  [m]
        z                   = delz:delz:zmax;               % depth vector [m]

        wcp = load(ssp_file);                               % Load water Sound Speed profile given by name wcpname
        qmultiSSP = isstruct(wcp);                          % T/F (1/0) is wcp a stuctured array?
        
        if qmultiSSP                                        % if wcp IS a stuctured array, read multiSSPs with following format
            c               = wcp.SSP{ssp_num};             % m/s water SS array. rows: depth veriation; Columns: different SS profiles
            c_z             = wcp.C_Z{ssp_num};             % (m) depth variation
            rho0            = 1.05;                         % water density [g/cc]
              
        else
            c_z             = wcp(:,1);                     % depth points for wcssp [m]
            c               = wcp(:,2);                     % sound speed points for wcssp [m/s]
            rho0            = mean(wcp(:,3));               % water density [g/cc]
        end
        
        zbottom             = c_z(end);                     % bottom depth
        z_bathy             = 1:delz:zbottom;
        a0                  = wateratten(freq);             % function: vector of frequency dependent water attenuation
        clear cofztemp cofz;
        cofztemp(:,1)       = interp1(c_z, c_z, z_bathy,'linear','extrap').';   % extrapolated depth points for wcssp [m]
        cofztemp(:,2)       = interp1(c_z, c, z_bathy,'linear','extrap').';     % extrapolated soundspeed for wcssp [m/S]
        cofz(:,1)           = linspace(cofztemp(1,1),cofztemp(end,1),250);      % interpolate new depth of water column profile
        cofz(:,2)           = interp1(cofztemp(:,1),cofztemp(:,2),cofz(:,1));   % interpolate ssp at new depth of water column profile
        
        if length(cofz) > 500
            error('Error - ORCA does not permit more than 500 wcssp points')
        end
        
        Basespeed           = 2350;                                         % sound speed in the basement 
        Bwss                = cofz(end,2);                                  % Bwss=bottom water sound speed
        
        svp_in.uphalf_cp    = 343.0;                                        % compressional sound speed in upper half space [m/s]
        svp_in.uphalf_cs    = 0.0;                                          % shear sound speed in upper half space [m/s]
        svp_in.uphalf_rho   = .00121;                                       % density in upper half space [g/cc]
        svp_in.uphalf_ap    = 0.000058;                                     % compressional attenuation in upper half space [pos~db/m/kHz, neg~dB/wavelength]
        svp_in.uphalf_as    = 0.0;                                          % shear attenuation in upper half space [dB/m/kHz]
        svp_in.nsvp         = length(cofz);                                 % number of SVP points in the ocean to be read
        svp_in.ctol         = 1.0;                                          % tolerance used in fitting SVP to eliminate layers (0 = keep all layers)
        svp_in.wssp         = cofz;                                         % SVP profile
        svp_in.wrho         = rho0;                                         % water density [g/cc]                                      % number of sediment layers to be read in
        svp_in.walphs       = a0;                                           % water attenuation [pos~db/m/kHz, neg~dB/wavelength]
        % svp_in.btm_env      = sedtemp;                                    % Commented out because this paramter will be updated during the Monte Carlo iterations 
        svp_in.lowhalf_cp   = Basespeed;                                    % compressional sound speed in lower halfpsace [m/s]
        svp_in.lowhalf_cs   = 0;                                            % shear sound speed in lower halfspace [m/s]
        svp_in.lowhalf_rho  = 2.6;                                          % density in lower halfspace [g/cc]
        svp_in.lowhalf_ap   = 0.22;                                       % compressional attenuation in lower halfspace [pos~db/m/kHz, neg~dB/wavelength]
        svp_in.lowhalf_as   = 0.0;                                          % shear attenuation in lower halfspace [pos~db/m/kHz, neg~dB/wavelength]
        svp_in.ntop         = 0;                                            % number of surface layers to be read in on next line
        svp_in.above_sea    = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
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
        opt_in.fcw          = freq;                                         % [f1 f2 ... ]
        opt_in.fcw_n        = length(opt_in.fcw);                           % size of the fcw array
        opt_in.nzm          = -length(z);                                   % depth points in modes (nzm>0 ==> List zm's; nf<0 ==> List first,last zm)
        opt_in.zm           = [z(1) z(end)];                                % [z1 z2 ... ]
        opt_in.zm_n         = length(opt_in.zm);                            % size of the zm array
        
        iimf = 1;                                               

        % Preallocate to improve speed of program 

        cphi                = zeros(length(z),nmodes);
        p2                  = ones(nTT,nfreq)*1e-12;                % initialize for speed & set pressure greens function array p2 to very small value
        p2Greens            = ones(nTT,nfreq)*1e-12;                % initialize for speed & set pressure greens function array p2 to very small value 
        ppparray            = ones(1,nTT,nfreq)*1e-12;              % initialized complex pressure array
        envinputsArray      = zeros(1,8);
        RRa                 = zeros(1,nTT);                         % RRa = 3D array of ranges (initialized here)                         
        sedtemp             = zeros(6,16);
        Vg                  = zeros(nmodes, nfreq);
        F                   = zeros(1,3);  
        V_point             = zeros(1,3);
        
        lowest_error        = 100;
        
        kn                  = zeros(nmodes,nfreq);
        
        if length(a0) == 1
            a0 = a0*ones(size(freq));
        end

        x0  = CPA * 1000;       % Source-reciever range (m)

        for iTT=1:nTT  % time index for each spectrogram 
            vsource     = sourceKnots*0.51444;              % Ship speed in m/s from Ship speed in knots
            vy          = vsource;                          % Source to Receiver Range speed [m/s] in y direction                 
            y0          = -vy*0.5*tspec(length(tspec));     % Symmetric about cpa
            YY          = y0+vy*tspec(iTT);
            RRa(iTT)    = (x0.^2+YY.^2).^0.5;               % Calculation Ranges [m] coresponding to Source to Receiver times  
        end

        Thicklayer1         = params(1);
        Layerspeed1         = params(2);
        Thicklayer2         = params(3);
        Layerspeed2         = params(4);

        for ifreq = 1:nfreq             % loop over Selected ORCA frequencies freq
            opt_in.fcw       = freq(ifreq);                 % define the frequency for this instant in ifreq loop
            svp_in.walphs    = a0(ifreq);                   % fequency dependent attenuation in water column
             
            sedtemp(1,:) = [1 9.2 1445 1446 0 0 1.62 1.62 0.04 0.04 0 0 0 0 0 0];
            sedtemp(2,:) = [1 3.0 1446 1750 0 0 1.8 1.8 0.15 0.15 0 0 0 0 0 0];
            sedtemp(3,:) = [1 7.5 1750 1750 0 0 1.83 1.83 0.15 0.15 0 0 0 0 0 0];
            sedtemp(4,:) = [1 Thicklayer1 Layerspeed1 Layerspeed1 0 0 2.0 2.0 0.15 0.15 0 0 0 0 0 0];
            sedtemp(5,:) = [1 Thicklayer2 Layerspeed2 Layerspeed2 0 0 2.2 2.2 0.15 0.15 0 0 0 0 0 0];
            sedtemp(6,:) = [1 100 2100 2100 0 0 2.2 2.2 0.15 0.15 0 0 0 0 0 0];
            
            svp_in.nlayb     = size(sedtemp,1);             % number of sediment layers to be read in
            svp_in.btm_env   = sedtemp; 
                                                      
            % Call to ORCA 
            [nmodes, kn_re, kn_im, vg, ~, phi_re, phi_im, ~] = sub_orca(svp_in, opt_in, iimf);  % Expensive call in time
            Vg(:,ifreq)             = vg;                                                       % construct vg 
            kn(1:nmodes,ifreq)      = kn_re(1:nmodes,1) - 1i*kn_im(1:nmodes,1);                 % construct kn
            cphi                    = phi_re - 1i*phi_im;                                       % construct cphi
        
            [~,idxs]                = min(abs(z-zs));                                           % find index of closest source depth on grid                        
            [~,idxr]                = min(abs(z-zr));                                           % find index of closest receiver depth on grid
            
            RR                      = transpose(RRa);
            sumTerm                 = zeros(nmodes,length(RR));
            cphi_s                  = cphi(idxs,1:nmodes);
            cphir                   = cphi(idxr,1:nmodes);
            km                      = kn(1:nmodes,ifreq);
                               
            for idx_mode = 1:nmodes
                sumTerm(idx_mode,:) = cphi_s(idx_mode).*cphir(idx_mode).*...
                sqrt(2./(pi.*(km(idx_mode).*RR))).*exp(-1i.*((km(idx_mode).*RR)-pi/4));
            end
            
            summation               = sum(sumTerm);
            p2Greens(:,ifreq)       = 1i*pi.*summation;                                         % construct the field at each frequency 
                               
        end % ifreq = 1:nfreq  


        RR                      = shiftdim(RRa(:)).';                   % ranges for short time varying time TT
        
        psrc                    = PSource(1,:);                         % complexx Source Pressure at depth # izs, for frequencies: freq
        ppp                     = (ones(nTT,1)*psrc).*p2Greens;         % pressure at depth ppp( = Source pressure * pressure greens function      
         
        TL                      = 20*log10(abs(p2Greens));
        
        spec_title = [ship,'_',plotting_info.array{k},'_model_CH',num2str(channel)];

        % Plots the spectrogram
        fig = figure('name', spec_title);
        

        TL=-TL; % makes transmission loss positive
        SL=mean(transpose(PSDselectf) + TL);
  
        TL_mea      = transpose(PSDselectf) + SL; % infer measured TL from the data and source model
        RL_mod      = SL - TL;
        clims(2)    = clim2;
        clims(1)    = clims(2)-60;
        
        imagesc(fselect, tspec, RL_mod, clims)
        colorbar
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

        ylim = get(gca, 'YLim');
        line([b12(1) b12(1)], ylim, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 2);
        % Labels the dashed line with the frequency
        [~, text_index] = min(abs(freq - b12(1))); % Find the index in the x-axis (plot_freq) closest to the target frequency
        y_coordinate = -60;  % Adjust the vertical position as needed
        text_str = ['$f_{1,2}$ = ', sprintf('%#.1f', b12(1))];
        text(freq(text_index), y_coordinate, text_str, 'Interpreter','latex','HorizontalAlignment', 'center', 'FontSize', 18); % Add the text above the spectrogram at the specified frequency
        
        save_name=[spec_title,'.png'];
        save_name=fullfile(git_library,'Spectrograms/Modeled/',save_name);
        saveas(fig, save_name)

        disp('Plotting')
        %%
        if plot_gv==true
            fr=real(kn);
            clear fprime;
            for n=1:nmodes
                fprime(n,1)=(1.0/24.0)*(1.0/df)*(-50.0*fr(n,1)+96.0*fr(n,2)-72.0*fr(n,3)+32.0*fr(n,4)-6.0*fr(n,5));
                fprime(n,2)=(1.0/24.0)*(1.0/df)*(-6.0*fr(n,2)-20.0*fr(n,3)+36.0*fr(n,4)-12.0*fr(n,5)+2.0*fr(n,6));
               
            for i=3:length(fselect)-2
            fprime(n,i)=(1.0/24.0)*(1.0/df)*(2.0*fr(n,i-2)-16.0*fr(n,i-1)+16.0*fr(n,i+1)-2.0*fr(n,i+2));
            end
             fprime(n,length(fselect)-1)=fprime(n,length(fselect)-2);
             fprime(n,length(fselect))=fprime(n,length(fselect)-2);
            end
            V=(2.0*pi ./fprime);
            MM=ones(nmodes,length(fselect));
            MM=2.0*pi*MM;
           
            clear Vinv Cinv;
            C=  MM ./fr;
            Cinv=1.0 ./C;
            Vinv=1.0 ./V;
             
            
            fr      = real(kn);
            MM      = 2.0*pi*ones(nmodes,length(fselect));      % converting from Hz to angular frequ 
            C       = MM ./fr;                                  % phase speed (c = omega/k) 
           
            betatop             = zeros(nmodes, nmodes, nfreq);          
            betabottom          = zeros(nmodes, nmodes, nfreq);
            for ii=1:nmodes
                for jj=1:nmodes
                    betatop(ii,jj,:)    = Cinv(ii,:)-Cinv(jj,:);
                    betabottom(ii,jj,:) = Vinv(ii,:)-Vinv(jj,:);
                end
            end
            
            Beta    = betatop ./betabottom;
            F = zeros(3,1);
            V_point= zeros(3,1);
            mode_crosses = [2,3,4];
            clear AA;
            % Finds the freq's where group velocity intersects and saves them as F(1), F(2), F(3), and F(4)
            for n=1:length(mode_crosses)
                mode_cross=mode_crosses(n);
                AA(:)                   = -Beta(1,mode_cross,:);
                %AA                      = AA ./transpose(fselect);
                %[pks,locs]              = findpeaks(AA);
                [Amax,imax]=max(abs(AA));
                F(n)=fselect(imax); 
                V_point(n)=V(1,imax);
            end 
    
            addInfo = {ship, array, F(1), F(2), F(3)};
            optimalInfo = [optimalInfo; addInfo];
    
            gv=figure;
            hold on;
            modes2plot=4;
            hL = legend;
            for j=1:modes2plot
                mode=V(j,:);
                plot(freq, mode, 'DisplayName', ['Mode ', num2str(j)]);     % Plots group velocity versus frequ. for each mode 
                
                set(hL,'AutoUpdate','off')
                if j<4
                    scatter(F(j),V_point(j),100,'k')
                end
                set(hL,'AutoUpdate','on')
            end
            xlim([17 30]);
            %ylim([1100 2100]);
            hold off;
            gv_title = [ship,'_',plotting_info.array{k},'_GV_CH',num2str(channel)];
            % Figure parameters
            legend('show');
            ylabel('Group velocity (m/s)')
            xlabel('Frequency (Hz)')
            
            plot_lab=['Group Velocity ', ship, ' ', array];
            title(plot_lab)
            save_name=[gv_title,'.png'];
            save_name=fullfile(git_library,'Spectrograms/GV (Modeled)/',save_name);
            saveas(gv, save_name)
            
        
        end
    end 
end 











                




