% *****************************************************************
%  *                                                              *
%  *     COPYRIGHT © 2019                                         *
%  *     Knobles Scientific and Analysis, LLC                     *                                      *
%  *     All rights reserved                                      *
%  *                                                              *
%  * Redistribution of source or binary forms is not permitted    *
%  * without the permission of KSA, LLC                           *
%  * THIS SOFTWARE IS PROVIDED BY KSA,LLC “AS IS” AND ANY EXPRESS *
%  * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    *
%  * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A      *
%  * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL KSA,LLC *
%  * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
%  * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT      *
%  * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
%  * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)     *
%  * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN    *
%  * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR *
%  * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS         *
%  * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. *
%  *                                                              *
%  ****************************************************************

%%
%   This program performs an inversion to find the thickness and sound
%   speed of the bottom three layers in the New England Mudpatch for 2017
%   and 2022. The parameters for the mudlayer, sand, transition, and
%   basement sediments are pre-set. 
%
%   This inversion technique works by finding the frequencies where beta
%   equals infinity (where the group velocity of mode one intersects modes
%   2, and 3) for each potential sediment and comparing it to the actual
%   beta equals infinity points from the 2017 or 2022 data and calculating
%   a squared error. 
%
%   The code is split in multiple sections where only the first needs to be
%   updated by the user with information regarding the particular ship,
%   VLA, SSP used, estimated upper and lower limits for each parameter
%   of interest in the inversion, and number of iterations to compute. 
%
%   Output is a modeled spectrogram using the optimal seabed parameters and
%   text display with values for the modeled beta equals infinity points
%   as well as sound speed and thickness for the optimal seabed parameters.
%   
%   This is modified code based off of 'sim_version_3po_featurebasedinversion.m' 
%   Author: David Knobles 
%   Modified by: Alexandra Hopps-McDaniel (11/22/23)

%%  start clock for execution, clear variables, & add paths
clear; %close all;  tic

addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\matlab-bborca')
addpath('C:\Users\alexh\OneDrive\Desktop\KSA Research\Codes\matlab-bborca\Inversion\Inversion_4000_iteration_scans\Ship Specs .mat files')

%%  Inversion parameters, updated by user 

% Specify the ship of interest and the given beta values from the data 
ship            = 'ALS APOLLO';           % Ship of interest
array           = 'VLA2';               % Specifiies which array is capturing the signal 
sourceKnots     = 16.3;                 % Ship speed [kn] 
CPA             = 1.6;                 % Closest point of approach [km]
betas           = [25.8,39.5];      % Frequency value of [beta12 beta13 beta14]
%beta_val = [24.6,39.9];
NMAX            = 4000;                   % Number of Monte Carlo iterations in the inversion 
nmodes          = 7;                   % Number of modes  
run_num         = '03_NEWMUD';
mud_tickness    = 7.2;
% Update the lower and upper limits for each parameter for the monte carlo iterations 
% The parameters are ordered as follows: 
% Layer 1 thickness; layer 1 sound speed; layer 2 thickness; ... layer 3 sound speed 
% Note: the closer the upper and lower bounds are to the actual seabed
% paramters, the more accurate the inversion results 
limits = [50 300; 1700 1860; 50 300; 1865 2050]; %450 550; 2000 2200];

wcpname         = 'SBCEX22_SSP';     % File name with 20 water column SSPs
ssp_num         = 3;%3;                 % Specifies which of the 20 SSP should be used

% Input array information
aa              = 3.75;             % Distance between channels in VLA
dist_channel1   = 7.25;             % Distance from ocean floor to the first hydrophone 
channel         = 1;                % Receiver channel being modeled 
water_depth     = 78;
zs              = 9.2;              % Source position (m) (note: the source comes from the propellers of the ship which depth can be estimated as the kiln depth found online)

zr              = water_depth - dist_channel1 - channel*aa;             % Reciever position (m)

% Load the .mat file that contains the ship spectrogram
% Use names saved into the .mat file for tspec fselect psdsectf (in
% decibels)
% Load in VLA2_featureInversion file 
load('ALS APOLLO_VLA2.mat','tspec','fselect','PSDselectf')
clims(2)=max(PSDselectf(:));
clims(1)=clims(2)-60;
%imagesc(fselect, tspec, PSDselectf, clims)

plot_gv = false;
%%  Time, freq, PSource, & z parameters: 

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

%% orca90 .svp file -------------------------------------------------------
% Create water column profile (maximum of 500 points allowed)

wcp = load(wcpname);                % Load water Sound Speed profile given by name wcpname
qmultiSSP = isstruct(wcp);          % T/F (1/0) is wcp a stuctured array?

if qmultiSSP                        % if wcp IS a stuctured array, read multiSSPs with following format
    c        = wcp.SSP{ssp_num};         % m/s water SS array. rows: depth veriation; Columns: different SS profiles
    c_z      = wcp.C_Z{ssp_num};             % (m) depth variation
    rho0     = 1.05;      % water density [g/cc]
      % c is currently 1st water ssp  (m/s)
else
    c_z      = wcp(:,1);            % depth points for wcssp [m]
    c        = wcp(:,2);            % sound speed points for wcssp [m/s]
    rho0     = mean(wcp(:,3));      % water density [g/cc]
end

zbottom         = c_z(end);                                         % bottom depth
z_bathy         = 1:delz:zbottom;
a0              = wateratten(freq);                                 % function: vector of frequency dependent water attenuation

cofztemp(:,1)   = interp1(c_z, c_z, z_bathy,'linear','extrap').';   % extrapolated depth points for wcssp [m]
cofztemp(:,2)   = interp1(c_z, c, z_bathy,'linear','extrap').';     % extrapolated soundspeed for wcssp [m/S]
cofz(:,1)       = linspace(cofztemp(1,1),cofztemp(end,1),250);      % interpolate new depth of water column profile
cofz(:,2)       = interp1(cofztemp(:,1),cofztemp(:,2),cofz(:,1));   % interpolate ssp at new depth of water column profile

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
opt_in.fcw          = freq;                                         % [f1 f2 ... ]
opt_in.fcw_n        = length(opt_in.fcw);                           % size of the fcw array
opt_in.nzm          = -length(z);                                   % depth points in modes (nzm>0 ==> List zm's; nf<0 ==> List first,last zm)
opt_in.zm           = [z(1) z(end)];                                % [z1 z2 ... ]
opt_in.zm_n         = length(opt_in.zm);                            % size of the zm array

iimf = 1;                                               

%% Preallocate to improve speed of program 

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
betatop             = zeros(nmodes, nmodes, nfreq);          
betabottom          = zeros(nmodes, nmodes, nfreq);
lowest_error        = 100;
num_params          = length(limits);
Parm                = zeros(1, num_params);
sq                  = zeros(1,length(betas));
dist                = zeros(NMAX, 5);
kn                  = zeros(nmodes,nfreq);
labels={'Thickness1'; 'SS1'; 'Thickness2'; 'SS2'};%'Thickness3'; 'SS3'};
info.lim = cell(1, num_params);
info.type = zeros(1, num_params);
info.label = cell(1, num_params);



for n=1:num_params
    info.lim{n}=limits(n,:);
    info.type(n)=1;
    info.label{n}=labels(n);
end 

if length(a0) == 1
    a0 = a0*ones(size(freq));
end

%% Multiloop to simulate 2d spectrogram Range array centered about each CPA range

x0  = CPA * 1000;       % Source-reciever range (m)

for iTT=1:nTT  % time index for each spectrogram 
    vsource     = sourceKnots*0.51444;              % Ship speed in m/s from Ship speed in knots
    vy          = vsource;                          % Source to Receiver Range speed [m/s] in y direction                 
    y0          = -vy*0.5*tspec(length(tspec));     % Symmetric about cpa
    YY          = y0+vy*tspec(iTT);
    RRa(iTT)    = (x0.^2+YY.^2).^0.5;               % Calculation Ranges [m] coresponding to Source to Receiver times  
end

%% Loop computing the beta equals infinity points for each monte carlo iteration
for NNN=1:NMAX

disp(['Monte Carlo iteration: ', num2str(NNN)])


for ilim = 1:num_params
    Parm(ilim) = (limits(ilim, 2)-limits(ilim, 1)).*rand(1) + limits(ilim, 1);
end 

% Defines the layer speed, mudspeed, and thicklayer with the random numbers defined above
Thicklayer1         = Parm(1);
Layerspeed1         = Parm(2);
Thicklayer2         = Parm(3);
Layerspeed2         = Parm(4);
%Thicklayer3         = Parm(5);
%Layerspeed3         = Parm(6);
    
isamp       = 0;

for ifreq = 1:nfreq             % loop over Selected ORCA frequencies freq
    opt_in.fcw       = freq(ifreq);                 % define the frequency for this instant in ifreq loop
    svp_in.walphs    = a0(ifreq);                   % fequency dependent attenuation in water column
     
    sedtemp(1,:) = [1 mud_tickness 1445 1446 0 0 1.62 1.62 0.04 0.04 0 0 0 0 0 0];
    sedtemp(2,:) = [1 3.0 1446 1750 0 0 1.8 1.8 0.15 0.15 0 0 0 0 0 0];
    sedtemp(3,:) = [1 7.5 1750 1750 0 0 1.83 1.83 0.15 0.15 0 0 0 0 0 0];
    %sedtemp(4,:) = [1 173.2 1760.3 1760.3 0 0 2.0 2.0 0.15 0.15 0 0 0 0 0 0];
    %sedtemp(5,:) = [1 185 1900 1900 0 0 2.2 2.2 0.15 0.15 0 0 0 0 0 0];
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

RR                      = shiftdim(RRa(:)).';               % ranges for short time varying time TT
isamp                   = isamp + 1; 
psrc                    = PSource(1,:);                     % complexx Source Pressure at depth # izs, for frequencies: freq
ppp                     = (ones(nTT,1)*psrc).*p2Greens;     % pressure at depth ppp( = Source pressure * pressure greens function      
ppparray(isamp,:,:)     = p2Greens(:,:);                    % array for each isamp 
TL                      = 20*log10(abs(p2Greens));
ppp_fullfield           = ppparray(:,:,1:nfreq);

%%

fr=real(kn);
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
%MM(1:25,:)=MM(1:25,:) * fselect(:);

 C=  MM ./fr;
 Cinv=1.0 ./C;
 Vinv=1.0 ./V;
 

fr      = real(kn);
MM      = 2.0*pi*ones(nmodes,length(fselect));      % converting from Hz to angular frequ 
C       = MM ./fr;                                  % phase speed (c = omega/k) 
%Cinv    = 1.0 ./C;                                  % inverse of phase speed 
%Vinv    = 1.0 ./Vg;                                 % inverse of the group speed

for ii=1:nmodes
    for jj=1:nmodes
        betatop(ii,jj,:)    = Cinv(ii,:)-Cinv(jj,:);
        betabottom(ii,jj,:) = Vinv(ii,:)-Vinv(jj,:);
    end
end

Beta    = betatop ./betabottom;

mode_crosses = [2,3];
% Finds the freq's where group velocity intersects and saves them as F(1), F(2), F(3), and F(4)
for n=1:length(mode_crosses)
    mode_cross=mode_crosses(n);
    AA(:)                   = -Beta(1,mode_cross,:);
    %AA                      = AA ./transpose(fselect);
    %[pks,locs]              = findpeaks(AA);
    [Amax,imax]=max(abs(AA));
    F(n)=fselect(imax);
    V_point(n)=Vg(1,imax);
    if n==2
        AA24=AA;
    end
    %F(n)                  = fselect(locs(length(locs)));
end 

%% Compute the error & save the optimal parameters

for i=1:length(betas)
    sq(i) = (F(i) - betas(i))^2;
end 

% Save parameters and error array for all iteration in a .mat file 
Error       = sum(sq);

% Add to dist file 
dist(NNN,:) = [Error Thicklayer1 Layerspeed1 Thicklayer2 Layerspeed2];% Thicklayer3 Layerspeed3];

% Saves the paramters and calculated TL of the iteration that gives the lowest error value
if Error < lowest_error
    optimal_TL              = TL; %save the TL for the lowest parameters 
    optimal_parameters      = [Thicklayer1 Layerspeed1 Thicklayer2 Layerspeed2];% Thicklayer3 Layerspeed3]; % update lowest error
    optimal_beta_vals       = [F(1) F(2) F(3)];
    lowest_error            = Error; % update lowest error so that at the end of the loop you have the lowest errors
end
     
end % End of monte caro iterations

%% Plotting the model 

figure

optimal_TL=-optimal_TL; % makes transmission loss positive
SL=mean(transpose(PSDselectf) + optimal_TL);

%SL          = S_0 - 10*log10(fselect.^3.594)+10*log10((1+(fselect/340).^2).^0.917); % Source Level - Wales-Heitmyer expression (2007) 
%SL          = transpose(SL);
TL_mea      = transpose(PSDselectf) + SL; % infer measured TL from the data and source model
RL_mod      = SL - optimal_TL;
clims(2)    = max(RL_mod(:));
clims(1)    = clims(2)-30;

imagesc(fselect, tspec, RL_mod, clims)
colorbar


%% Display ending information 

info.optimal_betas = optimal_beta_vals;
info.actual_betas = betas;
info.optimal_parm = optimal_parameters;
info.lowest_error = lowest_error;
info.num_iterations = NMAX;
cpa_info = ['CPA: ', num2str(x0/1000), ' km'];
knot_info = ['Speed: ', num2str(sourceKnots), ' kn'];
info.ship = {ship array cpa_info knot_info};
info.frequ = {['# Freqs: ', num2str(nfreq), ', Freq array = ', num2str(fmin),' : ',num2str(df),' : ',num2str(fmax),' Hz']};

disp('Program has completed.') 
disp(['Ship: ', ship, ', Array: ', array, ', ', cpa_info, ', ', knot_info])
disp(info.frequ{1})
disp(['# Monte Carlo Iterations: ', num2str(NMAX)])
%disp(['Actual beta values: ', num2str(betas(1)), ', ', num2str(betas(2)), ', ', num2str(betas(3)), ' Hz'])
disp(['Modeled beta values: ', num2str(optimal_beta_vals(1)), ', ', num2str(optimal_beta_vals(2)), ', ', num2str(optimal_beta_vals(3)), ' Hz'])
disp(['Optimal Parameters:', newline, 'DL1 = ', num2str(optimal_parameters(1)), ' m, ', num2str(optimal_parameters(2)), ' m/s', ...
    newline, 'DL2 = ', num2str(optimal_parameters(3)), ' m, ', num2str(optimal_parameters(4)), ' m/s']) %, newline, 'DL3 = ', ...
    %num2str(optimal_parameters(5)), ' m, ', num2str(optimal_parameters(6)), ' m/s'])
disp(['Error: ', num2str(lowest_error)])

elapsedTime = toc; 

%%
file_name = [ship, '_',array,'_dist_file_', num2str(NMAX), '_iterations_',run_num,'.mat'];
fig_name  = [ship, '_',array,'_model_spec_', num2str(NMAX), '_iterations_',run_num,'.png'];
saveas(gcf, fig_name)
save(file_name,'info','dist')
%%

if plot_gv == true
    figure
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
    hold off;
    
    % Figure parameters
    legend('show');
    ylabel('Group velocity (m/s)')
    xlabel('Frequency (Hz)')
    title('Group Velocity CARMEN VLA2')

end

%%
PLOT= false;
if PLOT==true
    hold on;
    
    colorbar
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
    %title(['Modeled Spectrogram: ', ship, ', ', array])
    
    ylim = get(gca, 'YLim');
                
                
    for i=1:length(beta_val)
        line([beta_val(i) beta_val(i)], ylim, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 1.5);
        % Labels the dashed line with the frequency
        [~, text_index] = min(abs(fselect - beta_val(i))); % Find the index in the x-axis (plot_freq) closest to the target frequency
        y_coordinate = -60;  % Adjust the vertical position as needed
        text_str = ['$f$ = ', sprintf('%#.1f', beta_val(i))];
        if i==2
            %text_index=text_index+35;
        end
        text(fselect(text_index), y_coordinate, text_str, 'Interpreter','latex','HorizontalAlignment', 'center', 'FontSize', 18); % Add the text above the spectrogram at the specified frequency
    end
end