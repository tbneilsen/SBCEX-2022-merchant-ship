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
%  ****************************************************************/
% 
% path(path,'D:\KSA\BBORCA-on-Dell') %my specific path to ORCA_MEX file: sub_orca
% cd D:\KSA\MyWorkspace-on-Dell
% 
% pscriptpath=mfilename('fullpath')  %#ok<NOPTS> %name of fullpath of currenewnt scipt
% return
% clear all; clc; close all;

%%  start clock for execution &       clear variables 
addpath(cd)
addpath('../ORCA_LIBRARY')
%load('kalamata_16')
%Data=kalamata_16;
%load('PSDselectf_04')
%PSDselectf=PSDselectf_env;
%load('SL_chirp')
%load('fselect')
%load('tspec')
tic
%clearvars ; close all; format compact,% before running, clear all variables & close all figures
%% Options for analysis, plotting, saving and acoustic source :
% nameroot='Ver12g2and5'; %root name of saved figures and data files
qts=0; % if ==1, model  time series
qpca=0;  % if ==1, peform and save PCA (Princicipal Component Aanalysis )
qimageplot=0;   %if ==1,  Plot Spectrogram Images
qsavefig=0; %if ==1, save figures to desktop, and then closes figure
qsavedat=0; % save ORCA data each time after ORCA completes
qRunOrca=1 ; % if ==1, run Orca (always should be 1)
qpressvsfreqplot=0;   % if ==1,plot individual |pressure| in dB   vs frequency (may generate many plots)
qtimeseriesPlots=0;  %if ==1, plot pressure  in microPa    vs time
qtestpressureplot=0;  %if ==1,For Testing: Plot Ship source pressure, then quit.
qRRvsTTplot=0; %if ==1, For Testing: Plot vector if Ship Ranges vs time TT, then quit
qtestPath=0             ; % if qtestPath: print name of current script & paths to dependent functions, the quit
qplotVGS=0; %if ==1,  Plot frequency dependent  VGS  sediment acoustic parameters; then quit
%%  options are set by choosing: qsource = -2,  -1, 0, 1 or 2 (search for "SourceNotes" for more details in this script)
 qsource = -1;     %% qsource = -2  generate mean ship spectrum + random enseble 
%qsource = -2;     %% qsource = -1  generate mean ship spectrum
% qsource=2; w=0.031;   %if  qsource =2 generate  SUS TNT pressure wave with  1st & 2nd bubble pulses;  w is Charge weight (kg) {for Mk64 w=0.031}
%% Important For-loop Calculation parameters
% RRv, zsa & zra arrays do Not require seperate ORCA runs:
%TT=(0:1.5171:9905.1);  % Spectrogram times [s] relative to  CPA Receiver Range. CPA range occurs
%TT=(0:0.5:8142);  % Spectrogram times [s] relative to  CPA Receiver Range. CPA range occurs
TT=(15.4998:15.4998:1844.48096);  % Spectrogram times [s] relative to  CPA Receiver Range. CPA range occurs
TT=tspec;
% sourceKnotsa= [5 7.5 10 12.5 15 17.5 20 22.5 25];% Vector of Source speed vector (kn); vsource=sourceKnots*0.51444; % Ship speed in m/s from Ship speed in knnots
aa=3.75;

sourceKnotsa= [19.9];    %KALAMATA speed from AIS data
sourceKnotsa= [2.6];    %chirp speed
sourceKnotsa= [14.8];    %viking bravery  speed 
sourceKnotsa= [2.8];    %chirp speed
sourceKnotsa= [2.8];    %chirp speed
 %sourceKnotsa= [18.9];    %nykdiana
 %sourceKnotsa= [19.9];    %KALAMATA speed from AIS data

 %sourceKnotsa= [11.5];    %matuama speed from AIS data
 sourceKnotsa= [17];    %viking bravery  spead 
  %sourceKnotsa= [19.9];    %KALAMATA speed from AIS data
  sourceKnotsa= [19.5];    %KALAMATA speed from AIS data
  %sourceKnotsa= [2.8];    %chirp spead


S_0=[237.0];   %y-intercept of Wales-Heitmyer container ship source levels broadband
%zra = [14,14+aa,14+2*aa,14+3*aa,14+4*aa,14+5*aa,14+6*aa,14+7*aa,14+8*aa,14+9*aa,14+10*aa,14+11*aa,14+12*aa,14+13*aa,14+14*aa,14+15*aa]; % vector of Receiver depths [m]  (Typically zbottom-0.5=79.5);
zra=[71.2-15.0*aa 71.2-14.0*aa 71.2-13.0*aa 71.2-12.0*aa 71.2-11.0*aa 71.2-10.0*aa 71.2-9.0*aa 71.2-8.0*aa 71.2-7.0*aa 71.2-6.0*aa 71.2-5.0*aa 71.2-4.0*aa 71.2-3.0*aa 71.2-2.0*aa 71.2-aa 71.2];

zra=71.2-5.0*aa    ;


 % Vector of CPA source-sensor  displacement ranges  in m

%zsa=[45]; %chirp
%zsa=[11];
zsa=[8.3];
zsa=[9.2];   %vikingbravery
%zsa=[10.3];
%zsa=[9.2];   %nykdiana
%zsa=[45]; %chirp

%x0a=[1186.5]; % Vector of CPA source-sensor  displacement ranges  in m

%x0a=[1027.5]; %chirp
%x0a=4215;
x0a=3500; %viking bravery  vla 1
x0a=3500; %viking bravery  vla 1
x0a=3187; %viking bravery  vla 2
x0a=3000; %viking bravery  vla 1
%x0a=2830; %viking bravery  vla 1
%x0a=3187; %viking bravery  vla 2
x0a=3087; %viking bravery  vla 2

x0a=3060; %kalamata  vla 2
%x0a=3000; %kalamata  vla 2

%x0a=[281.2]; %chirp vla2
%x0a=[8640]; %nykdiana
%x0a=[1296.0]; %chirp %chirp vla1

%x0a=3147; %viking bravery  vla 2

% NOTE: SSclassa & geoclassa require Expensive Orca calulation for each element (so keep vector size small for faster calculations
SSclassa= [2]; %%#ok<NBRAK> % %vector of water Soundspeed profile classes (read in wcpname)
NMAX=1;

        
        a=1750;
        b=1830;
        Parm_1=(b-a).*rand(NMAX,1)+a;
      
      
        a=1435;
        b=1749;
        Parm_2=(b-a).*rand(NMAX,1)+a;
        

        
        a=240;
        b=300;
       Parm_3=(b-a).*rand(NMAX,1)+a;
       
   %     a=0;
    %    b=10;
    %   Parm_4=(b-a).*rand(NMAX,1)+a;
       
      

        %%
        %%
        for NNN=1:NMAX
        %%

     
        Layerspeed1=Parm_1(NNN);
        mudspeed=Parm_2(NNN);
      Thicklayer1=Parm_3(NNN);

  
     grad=0;


              %matsuyama VLA 2 
%  Layerspeed1=1766.9;  
 % Thicklayer1=287.1733; 
  %  mudspeed=1438.8;
     Layerspeed2=1900;
     Layerspeed3=2100;
     Basespeed=2350;
     Thicklayer2=185.0;
     Thicklayer3=500.0;
     csand=1749.483;
     
 
     
     
  %corriedo 3 parm optimal VLA 2 
 % Layerspeed1=1770.996;  
 % Thicklayer1=291.34989; 
 %  mudspeed=1440.6517;
     Layerspeed2=1900;
     Layerspeed3=2100;
     Basespeed=2350;
     Thicklayer2=185.0;
     Thicklayer3=500.0;
     csand=1749.483;
     
     

     
    %viking bravery vla 2 optimal
    %Layerspeed1=1762.94244;  
    %Thicklayer1=291.2376; 
   %mudspeed=1441.94769;
     Layerspeed2=1900;
     Layerspeed3=2100;
     Basespeed=2350;
     Thicklayer2=185.0;
     Thicklayer3=500.0;
     csand=1749.483;
     
     
         %viking bravery vla 2 peak pdf
    %Layerspeed1=1771;  
     %Thicklayer1=292; 
   %mudspeed=1443;
     Layerspeed2=1900;
     Layerspeed3=2100;
     Basespeed=2350;
     Thicklayer2=185.0;
     Thicklayer3=500.0;
     csand=1749.483;
     

     
   %hafina green vla2 
  %  Layerspeed1=1766.067;  
  %   Thicklayer1=293.6046; 
  %  mudspeed=1442.7923;
     Layerspeed2=1900;
     Layerspeed3=2100;
     Basespeed=2350;
     Thicklayer2=185.0;
     Thicklayer3=500.0;
     csand=1749.483;
     
        %kalamata vla2 
 %   Layerspeed1=1788;  
  %   Thicklayer1=299.1538; 
  %  mudspeed=1446.8;
     Layerspeed2=1900;
     Layerspeed3=2100;
     Basespeed=2350;
     Thicklayer2=185.0;
     Thicklayer3=500.0;
     csand=1749.483;


geoclassa=[1];  %vector of  geo classes; possible values: 1,2,3,4, 5 defined by function sedtemp=sedimentssrb(geoclass,Bwss)(i)

%%  input paramaters that remain fixed throuhout Looping:
% wcpname = 'Squaredeal_array' ;  %%#ok<NASGU> % structured file name with 20 water column SSPs
% zmax    = 2000;              % maximum depth point for pressure field
wcpname = 'sbcexp2017_ssp' ; % structured file name with 20 water column SSPs
zmax    = 950;              % maximum depth point for pressure field
% delz= 0.1;         %#ok<NASGU> % depth mesh size  [m]
delz= 1;         % depth mesh size  [m]
% nmodes=800;  % prsetting number of modes analyzed
nmodes=35;  % prsetting number of modes analyzed
%fmin=5.00679016113; %  monimum frequency (Hz) of analysis, frequency resolution
fmin=fselect(1);
fmax=fselect(length(fselect));  %  maximum frequency (Hz) of analysis, frequency resolution df is controlled by max(x0a)
%fmax=199.9855;  %  maximum

df=(fselect(2)-fselect(1)) ;  % frequncy resolution 
%fmin=501; %  monimum frequency (Hz) of analysis, frequency resolution
%fmax=4001;  %  maximum frequency (Hz) of analysis, frequency resolution df is controlled by max(x0a)
T     =1/df ;   %time interval for ORCA analyses. Note: if T is undefined at this point, T will
%  later be redefined in terms of max(x0a), along with df, f, freq, nf,
%  nfreq after lines 
% if ~exist('T','var') %If T is not defined, then define T & redefine df, f, freq, nf, nfreq:
%     T = max(x0a(:))/1500 *1.4; %T is defined to encapsulate longest source to receiver travel time and
%     df      =1/T ;              % ORCA model frequency resolution  [Hz], delta f needed in sampling
%     f   = df:df:fmax;        % FFT frequencues evaluated [Hz]
%     freq    = fmin:df:fmax;        % ORCA frequencues evaluated [Hz]
%     fmax=max(freq);
%     nfreq=length(freq);   %Number of Orca frequencies evaluated
%     nf=length(f);   %Number of fft frequencies evaluated
% end  

%% VGS theory: frequency dependent  VGS  sediment acoustic parameters
% The input parameters are listed below. The example values used below
% are for the medium sand of SAX99. These values should be replaced
% with values for DK's sediment of interest.
dynrange=60; % dynamic range  (dB) of the plotted pressure spectrograms {= Max-Min of 20*log10(abs(ppp)) that is plotted }

%% ( code below this line should be left UNCHANGED 

% SourceNotes to generate Source spectrum: options are set by choosing: qsource = -2,  -1, 0, 1 or 2
% if qsource ==-2 % -2: generate  Ship Source spectrum from average fouier spectrum + added Gaussian Noise %  Wales's ship radiated noise paper [JASA 111 1211 (2002)]
% generate mean ship spectrum (in frequencny domain) + Gaussian Noise from Wales [JASA 111 1211 (2002)]
%     [PShipSig,PEnsemble]= shipSandN(ones(nzsa,1)*freq);
%     PSource=PShipSig +PEnsemble;  %Shipping Source Signal + complex random  normally distriuted ensemble, according to Wales
% PShipSig(freq) mean complex pressure power density (uPa/Hz) @ 1 m w/ 0 phase due to ensemble of ships 
% Wales' ship radiated noise paper [JASA 111 1211 (2002)], 
% PEnsemble(freq) is complex, normally distributed, random pressure based on Wales' standard deviation of ship radiated  levels
%
% elseif qsource ==-1
% generate mean ship spectrum with white Gaussian phase noise
% [PSource,S]= meanshipspectrum(freq);
%
% elseif qsource ==0  % generate impulse response, PSource=1 for all frequencies
% PSource = ones(size(freq));
% elseif qsource ==1 % generate Chapman_Dial_mk64SUS pressure wave with 1st bubble pulse
%
% elseif qsource ==2 % generate Chapman_Dial_mk64SUS pressure wave with 1st % 2nd bubble pulses
%
% elseif qsource ==3  % generate impulse response, PSource=1 for all frequencies
%     PSource = ones(size(freq));
% end %qsource logical
% Initialization 
% If qtestPath==1, This section of code displays name of current script & paths to dependent functions
%
%
%
if qtestPath % if qtestPath: print name of current script & paths to dependent functions, the quit
    pscriptname=mfilename; %name of current scipt
    [fList,pList] = matlab.codetools.requiredFilesAndProducts(pscriptname);
    disp([pscriptname 'Needs access to following Toolboxes:'])
    for i=1:length(pList)
        disp(pList(i))
    end
    disp(' ')
    disp(['And ',pscriptname,' Needs access to following paths\files:'])
    for i=1:length(fList)
        disp(fList(i))
    end
    return
end
% pscriptname=mfilename; %name of current script
nameroot=mfilename; %name of current script;
pscriptpath=mfilename('fullpath')  %#ok<NOPTS> %name of fullpath of current scipt
pscriptname=mfilename; %name of current script
nTT=length(TT); %number times for calculating each spectrogram  
iRR=0;
nzsa=length(zsa); %number of source depths
nzra=length(zra);  %number of receiver depths
ngeoa=length(geoclassa); %requires ORCA
nSSa=length(SSclassa);  %requires ORCA
nORCAruns=ngeoa*nSSa;  %#ok<NOPTS>
nx0=length(x0a); % nx0 = # of cpa ranges
nsh=length(sourceKnotsa);% nsh= shipping speeds
RRa=zeros(nsh,nx0,nTT); % RRa = 3D array of ranges (initialized here)
nTotalRuns=nTT*nx0*nsh*nzra*nzsa*ngeoa*nSSa;  %#ok<NOPTS>
nsamp=         nx0*nsh*nzra*nzsa*ngeoa*nSSa; %#ok<NOPTS>
nRunsPerOrca=nTT*nx0*nsh*nzsa*nzra  ; %#ok<NOPTS>
Rprime=RRa;% preallocating for speed
ix0ishp=0;
leg=cell(1,nx0*nsh); %legend
%
z       = delz:delz:zmax;        % depth vector [m]
% fmin=df; %  DO NOT CHANGE minimum frequency (Hz) of analysis
% df      =1/(max(x0a(:))/1500 *1.4) ; % frequency resolution (Hz)
if ~exist('fmin','var') %if fmin is not predefined, then 
    fmin=df; 
end % minimum frequency is frequency resolution (Hz)
%     T = max(RRa(:))/1500 *1.4; %T is defined to encapsulate longest source to receiver travel time and
%     df      =1/T ;              % ORCA model frequency resolution  [Hz], delta f needed in sampling
f   = df:df:fmax;        % FFT frequencues evaluated [Hz]
freq = fmin:df:fmax;        % ORCA frequencues evaluated [Hz]
freq=fselect;
if qplotVGS,return,end
fmax=max(freq);
nfreq=length(freq);   %Number of Orca frequencies evaluated
nf=length(f);   %Number of fft frequencies evaluated
%% multiloop to simulate 2d spectrogram Range array centered about each CPA range
for ix0=1:nx0 %# of cpa range loop
    x0=x0a(ix0);  %x0= cpa range in m; simulate Range RRa(ishipspd,ix0,iTT) & Rprime vs cpa range x0, ship speed & time (TT)
    for ishipspd=1:nsh %# of ship speeds
        sourceKnots=sourceKnotsa(ishipspd);  %Source Speed in Knots
        for iTT=1:nTT  % time index for each spectrogram (repeats for each CPA rane
            iRR=iRR+1; % Range index for each spectrogram (repeats for each CPA rane
            %             irr=(ix0*nsh)*(ishipspd)*nTT+iTT-nTT;
            %             disp([iRR irr])
            vsource=sourceKnots*0.51444; % Ship speed in m/s from Ship speed in knnots
            vy = vsource; %Source to Receiver Range speed [m/s] in y direction
                                   % y0= -vy*20.0*60.0*sourceKnots/19.47;  %source x position [m] in y direction at t=0;
  
            y0=-vy*0.5*tspec(length(tspec)); %symmetric about cpa
            %y0=-vy*tspec(length(tspec)); %do not go through CPa
     
            % GIVEN UNIFROM TIME ARRAY TT, COMPUTE RANGE ARRAY RR
            % RR=(x0+vx*TT); % Calculation Ranges [m] coresponding to Source to Receiver times
            %         XX=x0+vx*TT;
            YY=y0+vy*TT(iTT);
            RRa(ishipspd,ix0,iTT)=(x0.^2+YY.^2).^0.5; % Calculation Ranges [m] coresponding to Source to Receiver times
            %                 envinputs=[irun,x0,sourceKnots,zs,geoclass] ;
            Rprime(ishipspd,ix0,iTT)=(y0*vsource+vsource.^2*TT(iTT))./RRa(ishipspd,ix0,iTT);
        end
       
        if qRRvsTTplot  % Plot Ship range RRa(ishipspd,ix0,:) vs time TT(:)/60 and quit
            if ix0ishp==0 
                RRTTitle=[nameroot,' Ship range (km) vs time (min)']; % %#ok<UNRCH> %
                RRvsTTplot=figure('name', RRTTitle);
            end
            RR =shiftdim(RRa(ishipspd,ix0,:)).'; 
            ix0ishp=ix0ishp+1;
            leg{ix0ishp}=['Ship speed: ',num2str(sourceKnotsa(ishipspd)),' kn, CPA range: ', num2str(x0a(ix0)),' km'];
            plot(TT/60,RR/1e3 )% ,'.-')
            xlabel('TT - minutes'),ylabel('Range (km)'),grid
            hold on
            
            
        end %ifqpressvsfreqplot
    end
end
if qRRvsTTplot
    figure(RRvsTTplot) %#ok<*UNRCH>
    title(RRTTitle)
    legend(leg,'location','best')
    hold off
    return
end
%%
qShipPositionPlot=0; %if ==1, For Testing: Plot Ship Position Track Relative to Target, then quit, no longer used
if qShipPositionPlot %Test  Testing: Plot Ship Position Track Relative to Target, then quit
    RRTTitle='Ship Position Track Relative to Target'  ;
    h_position=figure('name', RRTTitle);
    subplot(1,1,1);
    plot(0,0,'R*',XX/1000,YY/1000,'x-');
    axis equal
    grid
    xlabel('X position (km)')
    ylabel('Y position (km)')
    title(  RRTTitle)
    %             disp('Program automatically stops after plotting source positions')
end
if ~exist('T','var') %If T is not defined, then define T & redefine df, f, freq, nf, nfreq:
    T = max(RRa(:))/1500 *1.4; %T is defined to encapsulate longest source to receiver travel time and
    df      =1/T   ;            % ORCA model frequency resolution  [Hz], delta f needed in sampling
    f   = df:df:fmax;        % FFT frequencues evaluated [Hz]
    freq    = fmin:df:fmax;        % ORCA frequencues evaluated [Hz]
%  freq=fselect;
    fmax=max(freq);
    nfreq=length(freq);   %Number of Orca frequencies evaluated
    nf=length(f);   %Number of fft frequencies evaluated
end  
% T=.1;
disp([num2str(nfreq),' ORCA freqs: ',num2str(fmin),' : ',num2str(df),' : ',num2str(fmax),' Hz'])


if qpca==1 % %#ok<UNRCH> % if qpca, peform and save PCA (Princicipal Component Aanalysis )
    latentarray=zeros(nTotalRuns,60); %initialize latten array for speed, (not entirely sure if this is always correct initial size)
end
wcp=load(wcpname);  % Load water Sound Speed profile given by name wcpname
qmultiSSP = isstruct(wcp); % T/F (1/0) is wcp a stuctured array?
if qmultiSSP % if wcp IS a stuctured array, read multiSSPs with following format
    c_array=wcp.c_array; % m/s water SS array. rows: depth veriation; Columns: different SS profiles,
    c_z=wcp.c_z;   % (m) depth variation
    rho0=mean(wcp.rho0);   % water density [g/cc]
    c=c_array(:,1); % c is currently 1st water ssp
    
else  % % if wcp IS NOT a stuctured array, read onlif qmultiSSP %
    SSclassa=SSclassa(1); % ensure  only  one SSp, with followind format
    wcp=load(wcpname);
    c_z     = wcp(:,1);        % depth points for wcssp [m]
    % The water depth is controlled by the last point in the bathymetry file.
    % We can change the water depth by defining z_bathy and then use linear
    % interpolation to create a new profile from the profile contained in
    %     zbottom=wcp(end,1);
    c       = wcp(:,2);        % sound speed points for wcssp [m/s]
    rho0    = mean(wcp(:,3));   % water density [g/cc]
end
zbottom=c_z(end); %bottom depth
z_bathy = 1:delz:zbottom;

% a0      = 0;                   % %constant water attenuation (previously named a0)
a0 = wateratten(freq); %function: vector of frequency dependent water attenuation
%% some preloop Bookeeping
% envinputsArray=cell(1,nTotalRuns);
ppparray= ones(nsamp,nTT,nfreq)*1e-12; % initialized  complex pressure array
envinputsArray=zeros(nsamp,8);
if length(a0)==1
    a0=a0*ones(size(freq));
end

% kn(1:nmodes,ifreq)      = kn_re(1:nmodes,1) - 1i*kn_im(1:nmodes,1);                     % construct kn
% cphi(1:length(z),1:nmodes,ifreq) = phi_re(:,1:nmodes) - 1i*phi_im(:,1:nmodes);                   % construct cphi
if qRunOrca % preallocate to speeed program
    % [tSource,Pt(izs,:),tdec,PtSUSdecim(izs,:),fSource,PSource0(izs,:),SSource0(izs,:)]=DecimatedChapmanSUS(T,zs,w);
    
    kn=zeros(nmodes,nfreq);
    cphi=zeros(length(z),nmodes);
    %     p2(izr,izs,iRR,ifreq)=ones(nzra,nzsa,nRRv,nfreq);
%     p2=ones(nzra,nzsa,nRRv,nfreq)*1e-12; %initialize for speed & set pressure greens function array p2 to very small value
    p2=ones(nzra,nzsa,nsh,nx0,nTT,nfreq)*1e-12; %initialize for speed & set pressure greens function array p2 to very small value
    p2Greens=ones(nTT,nfreq)*1e-12; %initialize for speed & set pressure greens function array p2 to very small value
%     p2pv=p2; %same initialization size for particle velocity Greensfunction array, p2pv
    
end
%% Generate  Source complex frequency spectrum: PSource
if qsource ==-2 % -2: generate  Ship Source spectrum from average fouier spectrum + added Gaussian Noise %  Wales's ship radiated noise paper [JASA 111 1211 (2002)]
    %% generate mean ship spectrum (in frequencny domain) + Gaussian Noise from Wales [JASA 111 1211 (2002)]
    [PShipSig,PEnsemble]= shipSandN(ones(nzsa,1)*freq);
    PSource=PShipSig +PEnsemble;  %Shipping Source Signal + complex random  normally distriuted ensemble, according to Wales
% PShipSig(freq) mean complex pressure power density (uPa/Hz) @ 1 m w/ 0 phase due to ensemble of ships 
% Wales' ship radiated noise paper [JASA 111 1211 (2002)], 
% PEnsemble(freq) is complex, normally distributed, random pressure based on Wales' standard deviation of ship radiated  levels
    [~,if1]    = min(abs(f-freq(1))); % index of fsource closest to freq(1)
    [~,if2]    = min(abs(f-freq(end))); % index of fsource closest to freq(end)

elseif qsource ==-1 % -1: generate  Ship Source spectrum from average fouier spectrum w/ random phase uniformly distributed: -pi to pi
    %% generate mean ship spectrum with random phase
    [PSource,SSource]= s00spectrum(freq,S_0); % function meanshipspectrum with random phase uniformly distributed between -pi  to pi
    PSource=ones(nzsa,1)*PSource;  %make shipping Source array consistent with  SUS Source array
    % PSource=real(PSource); %to remove random phase of ship source pressure, uncomment this  line
    % function [p,s]= meanshipspectrum(f)
    % PSource(freq)) mean complex pressure power density (uPa/Hz) @ 1 m w/ random phase
    % due to ensemble of ships
    % S((freq))= 20*log(|p((freq))|) Spectral Level in dB
    [~,if1]    = min(abs(f-freq(1))); % index of fsource closest to freq(1)
    [~,if2]    = min(abs(f-freq(end))); % index of fsource closest to freq(end)
elseif qsource >=0 && qsource <=2  % 0:2: generate Chapman_Dial_mk64SUS Source pressure wave with  qsource= 0 1, or 2 bubble pulses
    %% Generate SUS Source spectrum PSource at each source depth zsa
%     for izs=1:nzsa   % zs source depth loop in m
%         zs=zsa(izs);
%         [tSource,Pt(izs,:),tdec,PtSUSdecim(izs,:),fSource,PSource0(izs,:),SSource0(izs,:)]=DecimatedChapmanSUS(T,zs,w);
%     end %izs source depth loop
    for izs=1:nzsa   % zs source depth loop in m
        zs=zsa(izs);
        %[t,PtSus,tdecim,PtSUSdecim,fd,fftPtSusDec,Std]=DecimatedChapmanSUS(T,z0,w,Fsd,Nfir) {Fsd =5000,Nfir=30}
        [tSource,Pt(izs,:),tdec,~,fSource,PSource0(izs,:)]=DecimatedChapmanSUS(T,zs,w); %#ok<SAGROW>
    end %izs source depth loop
    N=length(Pt(1,:));
    [~,if1]    = min(abs(fSource-freq(1))); % index of fsource closest to freq(1)
    [~,if2]    = min(abs(fSource-freq(end))); % index of fsource closest to freq(end)
    PSource=PSource0(:,if1:if2); %confine pressure complex spectrum to freq(if1:if2) components
%     SSource=SSource0(:,if1:if2);  %confine pressure power spectrum to freq(if1:if2) components
%     ESource=ESource0(:,if1:if2); %confine pressure energy spectrum to freq(if1:if2) components
    if qpressvsfreqplot
        for izs=1:nzsa   % nzsa different source depths
            zs=zsa(izs); % zs source depth loop in m
            RRTTitle=['Estimated SUS Spectra at Source depth: ',num2str(zs),' m']; % %#ok<UNRCH> %
            hSusSourceFr=figure('name', RRTTitle);
            semilogx(fSource,20*log10(abs(PSource0)),'-+',fSource(if1:if2),20*log10(abs(PSource0(if1:if2))),'-o')
            legend('20Log10(|fft(Chapman TS)|)',['"filtered" over: ',num2str(fSource(if1)),'-',num2str(fSource(if2)),' Hz'])
            xlabel('Frequency (Hz) ' )
            grid
            ylabel('Pressure (dB re 1 \muPa^2')
            title(RRTTitle)
        end %izs source depth loop
        
    end %if qpressvsfreqplot
    if qts
        %                         tSource,Pt,fSource,PSource0
        RRTTitle=['Estimated SUS Spectra at Source depth: ',num2str(zs),' m']; % %#ok<UNRCH> %
        nullBelowifmin=min(abs(PSource0(:)))*1e-6*ones(nzsa,if1-1); % nulled  frequecy components below if1
        psrct0=ifft([PSource0,PSource0].','sym').';
        psrct=ifft([nullBelowifmin,PSource,PSource,nullBelowifmin].','sym').';
        Nt=length(psrct);
        Nt0=length(psrct0);
        t0=(1:Nt0)*T/Nt0;
        t=(1:Nt)*T/Nt;
        pppt=zeros(nTotalRuns,Nt); % initialize for speed
        if qtimeseriesPlots % plot SUS source model vs time
            hSusSourceTs=figure('name', RRTTitle); % %#ok<UNRCH>
            RRTTitle=['TS of SUS  at Source depth: ',num2str(zs),' m '];
            plot(tSource,Pt,'-+',t0,psrct0,'-x')
            legend('Chapman TS','ifft(fft(Chapman TS),sym)')
            xlabel('Time (s) ' )
            grid
            ylabel('Pressure (\muPa)')
            title(RRTTitle)
            
            RRTTitle=['TS of SUS  at Source depth: ',num2str(zs),' m '];
            hSusSourceTsF=figure('name', RRTTitle);
            plot(t,psrct,'-o')
            legend('ifft("filtered" fft(Chapman TS)),sym)')
            xlabel('Time (s) ' )
            grid
            ylabel('Pressure (\muPa)')
            title(RRTTitle)
        end %  if qtimeseriesPlots && qtimeseriesPlots1st
    end %  if qts
    %  return
    % end %if qsource >=0 && qsource <=2  % generate Chapman_Dial_mk64SUS
    
    
elseif  qsource ==3 % untested: unit qsource spectrum results in received pressure being greens function p2
    PSource = ones(size(freq));
end %qsource logical
while qtestpressureplot % For Testing: Plot Ship source pressure, then quit
    figspec=figure;
    plot(freq,20*log10(abs(PSource)))
    % semilogx(f,S),
    grid
    xlabel('Frequency (Hz)')
    ylabel('dB re 1 \muPa/Hz @ 1 m')
    title('S= Mean Acoustic Spectrum Level of an Ensemble of Ships')
    annotation(figspec,'textbox',...
        [0.208928571428571 0.807142857142858 0.617857142857143 0.102333333333338],...
        'String',['S=230-10*log10(f^3^.^5^9^4)+10*log10((1+(f/340)^2)^0^.^9^1^7)',...
        newline,'(see: JASA 111 1216 (2002).'],...
        'FitBoxToText','off');
    return
end

%% initialize orca input parameters that most of which won't change in subsequent looping (except svp_in.wssp & sedtemp)
if qRunOrca
    clearvars cofztemp cofz
    % Create the equivalent of ORCA.svp file
    
    % Create water column profile (maximum of 500 points allowed)
  
    if qmultiSSP % if wcp IS a stuctured array, read multiSSPs with following format
        c_z=wcp.c_z;   % (m) depth variation
        rho0=mean(wcp.rho0);   % water density [g/cc]
        c=c_array(:,1); % c is currently 1st water ssp  (m/s)
    end
    cofztemp(:,1) = interp1(c_z, c_z, z_bathy,'linear','extrap').'; % extrapolated depth points for wcssp [m]
    cofztemp(:,2) = interp1(c_z, c, z_bathy,'linear','extrap').';   % extrapolated soundspeed for wcssp [m/S]
    cofz(:,1) = linspace(cofztemp(1,1),cofztemp(end,1),250);    % interpolate new depth of water column profile
    cofz(:,2) = interp1(cofztemp(:,1),cofztemp(:,2),cofz(:,1)); % interpolate ssp at new depth of water column profile
    svp_in.nsvp         = length(cofz);                                         % number of SVP points in the ocean to be read
    svp_in.wssp         = cofz;                                                 % SVP profile
    svp_in.wrho         = rho0;                                                 % water density [g/cc]
    
    
    %     pofz = [rho0 a0];                                          % remainder of water column fluid properties (rho, ap)
    
    if length(cofz) > 500
        error('Error - ORCA does not permit more than 500 wcssp points')
    end
    Bwss=cofz(end,2); % Bwss=bottom water souund speed
%     sedtemp=sedimentssrb(1,Bwss);    
%     [cp,apdB,cs,asdB]= vgsdkfct(freq,porosityS,densityg,densitybw,bulkg,bulkbw,phi,depth,qplotVGS); %VGS sediment Therory
    svp_in.uphalf_cp    = 343.0;                                                % compressional sound speed in upper half space [m/s]
    svp_in.uphalf_cs    = 0.0;       
    % shear sound speed in upper half space [m/s]
    svp_in.uphalf_rho   = .00121;                                               % density in upper half space [g/cc]
    svp_in.uphalf_ap    = 0.000058;                                                  % compressional attenuation in upper half space [pos~db/m/kHz, neg~dB/wavelength]
    svp_in.uphalf_as    = 0.0;                                                  % shear attenuation in upper half space [dB/m/kHz]
    svp_in.ctol         = 1.0;                                                  % tolerance used in fitting SVP to eliminate layers (0 = keep all layers)
    svp_in.walphs       = a0;                                                   % water attenuation [pos~db/m/kHz, neg~dB/wavelength]
%     svp_in.nlayb        = size(sedtemp,1);                                      % number of sediment layers to be read in
%     svp_in.btm_env      = sedtemp;                                              % [type h cp1 cp2 cs1 cs2 rho1 rho2 ap1 ap2 as1 as2] (if type < 0, add fexpp fexps at end of line)
    svp_in.lowhalf_cp   =  Basespeed;                                                 % compressional sound speed in lower halfpsace [m/s]
    svp_in.lowhalf_cs   = 0;                                                    % shear sound speed in lower halfspace [m/s]
    svp_in.lowhalf_rho  = 2.5;                                                  % density in lower halfspace [g/cc]
    svp_in.lowhalf_ap   =0.2222;                                                  % compressional attenuation in lower halfspace [pos~db/m/kHz, neg~dB/wavelength]
    svp_in.lowhalf_as   = 0.0;                                                    % shear attenuation in lower halfspace [pos~db/m/kHz, neg~dB/wavelength]
    svp_in.ntop         = 0;                                                    % number of surface layers to be read in on next line
    svp_in.above_sea    = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    opt_in.nmode        = nmodes;                                               % number of modes
    opt_in.cphmin       = 0;                                                    % min phase speed (0=p-wave modes only; -1=seismic modes also);
    opt_in.cphmax       = 0;                                                    % max phase speed [0=use rmin; pos=speed; neg=max angle(deg)];
    opt_in.rmax         = 100   ;                                               % max range of interest in km (0=use S/R geom);
    opt_in.phfac        = 4;                                                    % phase step parm: step by 2*pi/phfac (set to 4-8, 0=default==>4);
    opt_in.dbcut        = 100;                                                  % modes weaker by db_cut ignored (set to 30-60, 0=default==>50);
    opt_in.Aih_l        = 135.001;     
    %   opt_in.Aih_l        = -1;
%
        % gradient lower h-space: 0=default,-1=homogeneous, >0=da_bar;
    opt_in.Aih_u        = 0;                                                   % gradient upper h-space: 0=default,-1=homogeneous, >0=da_bar;
    opt_in.nf           = 1;                                                    % frequencies (nf>0 ==> List fcw's; nf<0 ==> List first,last f)
    opt_in.fcw          = freq;                                                 % [f1 f2 ... ]
    opt_in.fcw_n        = length(opt_in.fcw);                                   % size of the fcw array
    opt_in.nzm          = -length(z);                                           % depth points in modes (nzm>0 ==> List zm's; nf<0 ==> List first,last zm)
    opt_in.zm           = [z(1) z(end)];                                        % [z1 z2 ... ]
    opt_in.zm_n         = length(opt_in.zm);                                    % size of the zm array
    iimf = 1;
end % if qRunOrca
fac1=(1.0/2.0)*(1/delz);     %factor in three-point depth derivative
% irun=0; %initialize orca run loop counter
irunOrka=0;
isamp=0;
%% MAIN ORCA loop: geoclass, SSclass
for igclass=1: ngeoa %  %=1:length(geoclassa)
    geoclass=geoclassa(igclass); % 
 %if geoclass ==1 
 %   grad=0;

% end 
                                    
                                   
    for iSSclass=1:nSSa %=1:length(SSclassa);
        SSclass=SSclassa(iSSclass); % ( if wcpname = 'SW54ssp' , 20) possible soundspeed classes, & structured file name with multi water column ssp
        % Create SSclass water column profile
        if qmultiSSP % if wcp IS a stuctured array, read multiSSPs with following format
            c=c_array(:,SSclass); % (m/s) c is now water ssp vector for  SSclass
            cofztemp(:,2) = interp1(c_z, c, z_bathy,'linear','extrap').';   % extrapolated soundspeed for wcssp [m/S]
            cofz(:,2) = interp1(cofztemp(:,1),cofztemp(:,2),cofz(:,1)); % interpolate ssp at new depth of water column profile
            svp_in.wssp         = cofz;                                                 % SVP profile
            if length(cofz) > 500
                error('Error - ORCA does not permit more than 500 wcssp points')
            end
        end
        Bwss=cofz(end,2); % Bwss=bottom water souund speed
        % call ORCA_MEX to find modal properties at each frequency
      
        for ifreq = 1:nfreq       % loop over Selected ORCA frequencies freq
            if qRunOrca %  make new time computation with sub_orca if qRunOrca==1 otherwise use previosly loaded values
                if ifreq==1
                    irunOrka=irunOrka+1;
                    disp(['Beginning Run # ',num2str(irunOrka),' of ',num2str(nORCAruns),' Orca Runs'])
                end
                F=freq(ifreq);   % define the frequency for this instant in ifreq loop
                opt_in.fcw       = F;   % define the frequency for this instant in ifreq loop
                svp_in.walphs = a0(ifreq); % fequency dependent attenuation in water column
                if geoclass <= 15
                 %   sedtemp=sedimentssrb(geoclass,Bwss);   % Function sedimentssrb(geoclass,Bwss)  defines geology class (sediment properties
                 %   svp_in.nlayb        = size(sedtemp,1);                                      % number of sediment layers to be read in
                %elseif geoclass ==5 %%%    mud over sand ENV2 with  VGS sound speed & attenuation and the GS shear
                 %   FackHz=1000.0/F;
                    %[cpM,apdBM1000,csM,asdB1000M]= vgsdkfct2(1000,porosityM,densitygM,densitybwM,bulkgM,bulkbwM,phiM,depthM,nM,gmpoM,gmsoM,tauM);
                    %[cpS,apdBS1000,csS,asdBS1000]= vgsdkfct2(1000,porosityS,densityg,densitybw,bulkg,bulkbw,phi,depth,n,gmpo,gmso,tau);
               %     [cpM1,apdBM1,csM1,asdBM1]= vgsdkfct2(F,porosityM1,densityM1,densitybw,bulkgM1,bulkbw,phiM1,depth,nM1,gmpoM1,gmsoM1,tauM1);
                %    [cpM2,apdBM2,csM2,asdBM2]= vgsdkfct2(F,porosityM2,densityM2,densitybw,bulkgM2,bulkbw,phiM2,depth,nM2,gmpoM2,gmsoM2,tauM2);
                %    [cpT1,apdBT1,csT1,asdBT1]= vgsdkfct2(F,porosityT1,densityT1,densitybw,bulkgT1,bulkbw,phiT1,depth,nT1,gmpoT1,gmsoT1,tauT1);
                %    [cpT2,apdBT2,csT2,asdBT2]= vgsdkfct2(F,porosityT2,densityT2,densitybw,bulkgT2,bulkbw,phiT2,depth,nT2,gmpoT2,gmsoT2,tauT2);
                %    [cpT3,apdBT3,csT3,asdBT3]= vgsdkfct2(F,porosityT3,densityT3,densitybw,bulkgT3,bulkbw,phiT3,depth,nT3,gmpoT3,gmsoT3,tauT3);
                %    [cpS,apdBS,csS,asdBS]= vgsdkfct2(F,porosityS,densityS,densitybw,bulkgS,bulkbw,phiS,depth,nS,gmpoS,gmsoS,tauS);
                   
              
 
                  %  sedtemp(1,:) = [1 L1Thickness cpM1 cpM1+grad*L1Thickness csM1 csM1 densityM1/1000. densityM1/1000. (apdBM1)*FackHz (apdBM1)*FackHz asdBM1*FackHz asdBM1*FackHz 0 0 0 0]; 
                  
          
                  %  sedtemp(2,:) = [1 L2Thickness cpM2 cpM2+grad*L2Thickness csM2 csM2 densityM2/1000. densityM2/1000. (apdBM2)*FackHz (apdBM2)*FackHz asdBM2*FackHz asdBM2*FackHz 0 0 0 0];  
    
                    
                  
                  %  sedtemp(3,:) = [1 LT1Thickness cpT1  cpT1 csT1 csT1 densityT1/1000. densityT1/1000. (apdBT1)*FackHz (apdBT1)*FackHz asdBT1*FackHz asdBT1*FackHz 0 0 0 0];
                    


                 %  sedtemp(4,:) = [1 LT2Thickness cpT2  cpT2 csT2 csT2 densityT2/1000. densityT2/1000. (apdBT2)*FackHz (apdBT2)*FackHz asdBT2*FackHz asdBT2*FackHz 0 0 0 0];

 
                %   sedtemp(5,:) = [1 LT3Thickness cpT3  cpT3 csT3 csT3 densityT3/1000. densityT3/1000. (apdBT3)*FackHz (apdBT3)*FackHz asdBT3*FackHz asdBT3*FackHz 0 0 0 0];

                  
               %    sedtemp(6,:) = [1 L4Thickness cpS  cpS csS csS densityS/1000. densityS/1000. (apdBS)*FackHz (apdBS)*FackHz asdBS*FackHz asdBS*FackHz 0 0 0 0];
                   
               %    sedtemp(7,:) = [1 8200 1830  1831 0 0 2000/1000. 2000./1000. .205 .205 0 0 0 0 0 0];
          
                 %   sedtemp(2,:) = [1 L4Thickness cpS  cpS csS csS densityS/1000. densityS/1000. (apdBS)*FackHz (apdBS)*FackHz asdBS*FackHz asdBS*FackHz 0 0 0 0];
                 
                  
                  
      %        sedtemp(1,:) = [1 9.2 1445 1530 0 0 1.612 1.612 .04 .04 0 0 0 0 0 0]; 
       %       sedtemp(2,:) = [1 3 1530 1710 0 0 1.7 1.7 .15 .15 0 0 0 0 0 0];
        %      sedtemp(3,:) = [1 7.5 1750 1750 0 0 1.83 1.83 0.15 0.15 0 0 0 0 0 0];   
         %     sedtemp(4,:) = [1 Thicklayer1 Layerspeed1 Layerspeed1 0 0 1.84 1.84 0.05 .05  0 0 0 0 0 0];
          %     sedtemp(5,:) = [1 200 1830 1830 0 0 2.0 2.0 0.05 .15  0 0 0 0 0 0];
          %      sedtemp(6,:) = [1 200 2000 2000 0 0 2.2 2.21 0.15 0.15  0 0 0 0 0 0];
                
                
            %  sedtemp(1,:) = [1 9.2 1445 1446 0 0 1.612 1.612 .04 .04 0 0 0 0 0 0]; 
           %   sedtemp(2,:) = [1 3 1446 csand 0 0 1.7 1.7 .15 .15 0 0 0 0 0 0];
            %  sedtemp(3,:) = [1 7.5 csand csand 0 0 1.83 1.83 0.15 0.15 0 0 0 0 0 0];   
            %  sedtemp(4,:) = [1 Thicklayer1 Layerspeed1 Layerspeed1 0 0 1.84 1.84 0.05 .05  0 0 0 0 0 0];
            %   sedtemp(5,:) = [1 Thicklayer2 1830 1830 0 0 2.0 2.0 0.05 .15  0 0 0 0 0 0];
             %   sedtemp(6,:) = [1 Thicklayer2 2000 2000 0 0 2.2 2.21 0.15 0.15  0 0 0 0 0 0];
             
             
                
                %model 1   mud attenuation  10/21/2021
            %    sedtemp(1,:) = [1 9.2 1445 1530 0 0 1.612 1.612 .04 .04 0 0 0 0 0 0]; 
           %   sedtemp(2,:) = [1 3 1530 csand 0 0 1.7 1.7 .15 .15 0 0 0 0 0 0];
            %  sedtemp(3,:) = [1 7.5 csand csand 0 0 1.83 1.83 0.15 0.15 0 0 0 0 0 0];   
            %  sedtemp(4,:) = [1 Thicklayer1 Layerspeed1 Layerspeed1 0 0 1.84 1.84 0.05 .05  0 0 0 0 0 0];
            %   sedtemp(5,:) = [1 Thicklayer2 Layerspeed2 Layerspeed2 0 0 2.4 2.4 0.05 .05  0 0 0 0 0 0];
            %    sedtemp(6,:) = [1 Thicklayer3 Layerspeed3 Layerspeed3 0 0 2.5 2.5 0.05 0.05  0 0 0 0 0 0];
                
                                %model 2   mud attenuation  10/21/2021
             %   sedtemp(1,:) = [1 9.2 1445 1530 0 0 1.612 1.612 .055 .055 0 0 0 0 0 0]; 
            %  sedtemp(2,:) = [1 3 1530 csand 0 0 1.7 1.7 .15 .15 0 0 0 0 0 0];
           %   sedtemp(3,:) = [1 7.5 csand csand 0 0 1.83 1.83 0.15 0.15 0 0 0 0 0 0];   
           %   sedtemp(4,:) = [1 Thicklayer1 Layerspeed1 Layerspeed1 0 0 1.84 1.84 0.05 .05  0 0 0 0 0 0];
           %    sedtemp(5,:) = [1 Thicklayer2 Layerspeed2 Layerspeed2 0 0 2.4 2.4 0.05 .05  0 0 0 0 0 0];
           %     sedtemp(6,:) = [1 Thicklayer3 Layerspeed3 Layerspeed3 0 0 2.5 2.5 0.05 0.05  0 0 0 0 0 0];
                
                
                                                %model 3   mud attenuation  10/21/2021
          %      sedtemp(1,:) = [1 9.2 1445 1446 0 0 1.612 1.612 .040 .040 0 0 0 0 0 0]; 
          %    sedtemp(2,:) = [1 3 1446 csand 0 0 1.7 1.7 .15 .15 0 0 0 0 0 0];
          %    sedtemp(3,:) = [1 7.5 csand csand 0 0 1.83 1.83 0.15 0.15 0 0 0 0 0 0];   
          %    sedtemp(4,:) = [1 Thicklayer1 Layerspeed1 Layerspeed1 0 0 1.84 1.84 0.05 .05  0 0 0 0 0 0];
          %     sedtemp(5,:) = [1 Thicklayer2 Layerspeed2 Layerspeed2 0 0 2.4 2.4 0.05 .05  0 0 0 0 0 0];
          %      sedtemp(6,:) = [1 Thicklayer3 Layerspeed3 Layerspeed3 0 0 2.5 2.5 0.05 0.05  0 0 0 0 0 0];
                
                
                                                                %model 4   mud attenuation  10/21/2021
                          c2=mudspeed+9.2*9.2; 
                          c2=mudspeed;
                sedtemp(1,:) = [1 9.2 mudspeed c2 0 0 1.612 1.612 .040 .040 0 0 0 0 0 0]; 
              sedtemp(2,:) = [1 3 c2 Layerspeed1 0 0 1.7 1.7 .15 .15 0 0 0 0 0 0];  
              sedtemp(3,:) = [1 Thicklayer1 Layerspeed1+2 Layerspeed1+2 0 0 1.84 1.84 0.05 .05  0 0 0 0 0 0];
               sedtemp(4,:) = [1 Thicklayer2 Layerspeed2 Layerspeed2 0 0 2.4 2.4 0.05 .05  0 0 0 0 0 0];
               sedtemp(5,:) = [1 Thicklayer3 Layerspeed3 Layerspeed3 0 0 2.5 2.5 0.05 0.05  0 0 0 0 0 0];
           
              mudspeed=1443.0
            ch=mudspeed
              grad=15
               beta=-0.50;
               T1=9.2;
               zzz=T1
               AA1=((1+beta)^2)*ch^2;
               AA2=2.0*ch*(1+beta)*grad*zzz
           
              cbot=sqrt(AA1+AA2) -ch*beta
            
                h=1.6705;
               zt=0.4567;
               % zt=0.72;
                cinfty=1530;
                ch=mudspeed;
                T1=0.5;
                zzz=T1+h;
                A=1+(zt/zzz)^2;
                B=A/(cinfty)^2;
               cbot=sqrt(1/B)
               
                %              AA1=((1+beta)^2)*ch^2;
               %AA2=2.0*ch*(1+beta)*grad*zzz
               %cbot=sqrt(AA1+AA2) -ch*beta;
                
              sedtemp(1,:) = [1 T1 ch cbot 0 0 1.612 1.612 .040 .040 0 0 0 0 0 0]; 
              
                T2=T1+0.5;
                zzz=T2+h;
                A=1+(zt/zzz)^2;
              B=A/(cinfty)^2;
                cbot2=sqrt(1/B)
                                         AA1=((1+beta)^2)*ch^2;
               AA2=2.0*ch*(1+beta)*grad*zzz
               %cbot2=sqrt(AA1+AA2) -ch*beta;
              
              
              sedtemp(2,:) = [1 T1 cbot cbot2 0 0 1.612 1.612 .04 .04 0 0 0 0 0 0];  
              
                T3=T2+0.5;
                zzz=T3+h;
                A=1+(zt/zzz)^2;
               B=A/(cinfty)^2;
               cbot3=sqrt(1/B)
                                        AA1=((1+beta)^2)*ch^2;
               AA2=2.0*ch*(1+beta)*grad*zzz
               %cbot3=sqrt(AA1+AA2) -ch*beta;
               
            sedtemp(3,:) = [1 T1 cbot2 cbot3 0 0 1.612 1.612 0.04 .04  0 0 0 0 0 0];
              
               T4=T3+1;
                zzz=T4+h;
               A=1+(zt/zzz)^2;
               B=A/(cinfty)^2;
               cbot4=sqrt(1/B)
                                        AA1=((1+beta)^2)*ch^2;
               AA2=2.0*ch*(1+beta)*grad*zzz
               %cbot4=sqrt(AA1+AA2) -ch*beta;
                
               sedtemp(4,:) = [1 1 cbot3 cbot4 0 0 1.612 1.612 0.04 .04  0 0 0 0 0 0];
               T5=T4+7.2;
               zzz=T5+h;
              A=1+(zt/zzz)^2;
               B=A/(cinfty)^2;
              cbot5=sqrt(1/B)
                                       AA1=((1+beta)^2)*ch^2;
               AA2=2.0*ch*(1+beta)*grad*zzz
               %cbot5=sqrt(AA1+AA2) -ch*beta;
               
              sedtemp(5,:) = [1 6.7 cbot4 cbot5 0 0 1.612 1.612 0.04 0.04  0 0 0 0 0 0];
                              sedtemp(6,:) = [1 3 cbot5 Layerspeed1 0 0 1.7 1.7 .15 .15 0 0 0 0 0 0];  
%
             sedtemp(7,:) = [1 Thicklayer1 Layerspeed1+2 Layerspeed1+2 0 0 1.84 1.84 0.14 .14  0 0 0 0 0 0];
               sedtemp(8,:) = [1 Thicklayer2 Layerspeed2 Layerspeed2 0 0 2.4 2.4 0.15 0.15  0 0 0 0 0 0];
               sedtemp(9,:) = [1 Thicklayer3 Layerspeed3 Layerspeed3 0 0 2.5 2.5 0.15 0.15  0 0 0 0 0 0];
                
                
              
              
             
       
            svp_in.nlayb = size(sedtemp,1);
                end
                svp_in.btm_env      = sedtemp;  % [type h cp1 cp2 cs1 cs2 rho1 rho2 ap1 ap2 as1 as2] (if type < 0, add fexpp fexps at end of line)
                vp_in.nlayb        = size(sedtemp,1);                                      % number of sediment layers to be read in
             
                [nmodes, kn_re, kn_im, ~, ~, phi_re, phi_im, ~] = sub_orca(svp_in, opt_in, iimf); %Expensive call in time
                %       kn      = kn_re(1:nmodes,1) - 1i*kn_im(1:nmodes,1);                     % construct kn
                %       cphi     = phi_re(:,1:nmodes) - 1i*phi_im(:,1:nmodes);                   % construct cphi
                kn(1:nmodes,ifreq)      = kn_re(1:nmodes,1) - 1i*kn_im(1:nmodes,1);                     % construct kn
                %                        cphi(1:length(z),1:nmodes,ifreq) = phi_re(:,1:nmodes) - 1i*phi_im(:,1:nmodes);                   % construct cphi
                %                         cphi(1:length(z),1:nmodes) = phi_re(:,1:nmodes) - 1i*phi_im(:,1:nmodes);                   % construct cphi
                cphi = phi_re - 1i*phi_im;                   % construct cphi
                rhoomeg=2.0*pi*freq(ifreq)*rho0*1000.0*1i;
                %                 for iRR=1:nRRv  %loop for x0 source positions along x axis  in m             
                for ix0=1:nx0 %# of cpa range loop
                    % x0=x0a(ix0);  %x0= cpa range;
                    for ishipspd=1:nsh % source speed loop
                        RR=shiftdim(RRa(ishipspd,ix0,:)).';
                        for izs=1:nzsa  % izs loop for source depth zs 
                            zs=zsa(izs);  % zs source depth loop in m
                            [~,idxs]    = min(abs(z-zs)); % find index of closest source depth on grid
                            for izr=1:nzra %receiver depth loop
                                zr=zra(izr);  %receiver depth in m
                                [~,idxr]    = min(abs(z-zr));  % find index of closest receiver depth on grid
                                %for iTT=1:nTT  % time index for each spectrogram (repeats for each CPA rane 
                                    %                             RR=RRa(ishipspd,ix0,iTT);
                                %    kr       = kn(1:nmodes,ifreq)*RR(iTT);                                           % construct kr = kn*r
                                %    H       = sqrt(2./(pi*kr.')).*exp(-1i*(kr-pi/4)).';                % large argument approximation to the Hankel function
                                    % p2 is  the Green's function that results from a unit point source
                                %    p2(izr,izs,ishipspd,ix0,iTT,ifreq) = 1i*pi*sum( cphi(idxs,1:nmodes).*cphi(idxr,1:nmodes).*H ); % construct the field at each frequency
                                    %                             p2pv(izr,izs,iRR,ifreq) = fac1*1i*pi*sum( (-1.0*cphi(idxr-1,1:nmodes)+1.0*cphi(idxr+1,1:nmodes)).*cphi(idxs,1:nmodes).*H )*(1.0/rhoomeg); % construct field  p2 at each RANGE & frequency
                               % end %  iTT=1:nTT, RR=RRa(ishipspd,ix0,iTT);
                               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               % UDEL Option 3
                               sumTerm=zeros(nmodes,length(RR));
                               cphi_s=cphi(idxs,1:nmodes);
                               cphir=cphi(idxr,1:nmodes);
                               km=kn(1:nmodes,ifreq);
                               for idx_mode=1:nmodes
                                   sumTerm(idx_mode,:)=cphi_s(idx_mode).*cphir(idx_mode).*...
                                       sqrt(2./(pi.*(km(idx_mode).*RR))).*exp(-1i.*((km(idx_mode).*RR)-pi/4));
                               end
                               summation=sum(sumTerm);
                               p2(izr,izs,ishipspd,ix0,:,ifreq)=1i*pi.*summation; % construct the field at each frequency
                               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            end  % izr loop for receiver depth  zr in m
                        end  % izs loop for source depth  zs in m
                    end %  ishipspd=1:nsh, sourceKnots=sourceKnotsa(ishipspd); %for Source speeds in knots
                end % ix0=1:nx0 %loop for x0 source positions along x axis  in m
            end  % end if qRunOrca
        end     % ifreq = 1:nfreq  % loop over frequencies
        
        disp(' sub_orca FREQ Loop has ended.')
         %%   now vary parameters that are independent of wavenumber: zr, zs, RR, x0
        dataname=[nameroot,'_',int2str(irunOrka)];
        isp=strfind(dataname,'_');
        txtname=dataname;
        txtname(isp)=' ';
        % save ORCA Modeled  data to matfile:
        if qsavedat
            disp(['Saving ORCA results into file: ' dataname '.mat'])
            clearvars str
            str.dataname=dataname;
            str.irunOrka=irunOrka;
            str.pscriptpath=pscriptpath;
            str.pscriptname=pscriptname;
            str.elapsedTime=elapsedTime;
            str.RRa=RRa;
            str.zsa=zsa;
            str.zra=zra;
            str.SSclass=SSclass;
            str.geoclass=geoclass;
            str.kr=kr;
            str.H=H;
            str.PSource=PSource;
            str.qsource=qsource;
            str.freq=freq;
            str.f=f;
            str.sedtemp=sedtemp;
            str.freq=freq;
            str.if1=if1;
            str.if2=if2;
            eval([dataname,'=str;'])
            save(dataname,dataname)
        end
        
        %         disp(['Beginning Run # ',num2str(irunOrka),' of ',num2str(nORCAruns),' Orca Runs'])
        if qRunOrca % display progress
            %                         x0 %#ok<NOPTS>
            %                         sourceKnots %#ok<NOPTS>
            geoclass   %#ok<NOPTS>
            SSclass   %#ok<NOPTS>
            disp([num2str(nfreq),' freqs: ',num2str(fmin),' : ',num2str(df),' : ',num2str(fmax),' Hz'])
            disp(['completed ORCA run # ',num2str(irunOrka),' of ',num2str(nORCAruns)])
            disp(['completed run # ',num2str(irunOrka*nRunsPerOrca),' of ',num2str(nTotalRuns)])
        end
        ienv=0;

            for izs=1:nzsa   %  source depth loop 
                zs=zsa(izs); % zs source depth  in m
                for ix0=1:nx0 % cpa Range loop
                    x0=x0a(ix0);  % simulate varying cpa Ranges x0
                    for ishipspd=1:nsh% length(sourceKnotsa) % simulate cpa Ranges in m for varying ship speed in knots
                        sourceKnots=sourceKnotsa(ishipspd);  % simulate Ranges for varying ship speeds sourceKnots
                         for izr=1:nzra %receiver depth loop 
                         zr=zra(izr);  % zr receiver depth in m
                        %for iTT=1:nTT
                        %RRa(ishipspd,ix0,iTT)=(x0.^2+YY.^2).^0.5; % Calculation Ranges [m] coresponding to Source to Receiver times
                        RR=shiftdim(RRa(ishipspd,ix0,:)).'; %RRa(ishipspd,ix0,iTT); ranges for short time varying time TT
                        isamp=isamp+1  %#ok<NOPTS>
                        ienv=ienv+1;
%                         irun=irun+ 1;% length(RR);                        
                        isamp_geoclass_SSclass_ienv_izr_izs_ix0_ishipspd=[isamp, geoclass, SSclass, ienv izr izs ix0 ishipspd];
%                         envinputsArrayEntries='isamp,geoclass,SSclass,zs(m),zr(m),sourceKnots(kn),x0/1000(km),max(RR)/1000(km)';
                        envinputsArray(isamp,:)=[isamp,geoclass,SSclass,zs,zr,sourceKnots,x0/1000,max(RR)/1000];
                        envtitle=[txtname,': zs= ',num2str(zs),', qs=',num2str(qsource),', zr= ' ...
                            ,num2str(zr),' m, geo #: ',num2str(geoclass),', ',num2str(sourceKnots),...
                            ' kn, svp #: ',num2str(SSclass),', RR=',num2str(min(RR)/1000),':',num2str(max(RR)/1000),' km'];
%                         txtenvinputs='{pscriptpath,envtitle,irunOrka,RR,TT,zs,zr,geoclass,SSclass,freq,if1,ifmin,if2}'; 
%                         envinputs={txtenvinputs,envtitle,irunOrka,RR,TT,zs,zr,geoclass,SSclass} ; % cell array containing environmental parameters
%                                           p2(izr,izs,ishipspd,ix0,iTT,ifreq)
%                         p2Greens = shiftdim(p2(izr,izs,ishipspd,ix0,:,:)).' ; %shift dimensions to a vector for p2(izr,izs,ishipspd,ix0,iTT,ifreq)
%                         p2Greens(:,1:if1-1) = min(p2Greens(:))*1e-12*ones(size(p2Greens(:,1:if1-1))); %Define Greeens to be very small for frequencies < fmin
%%  change made 9/4/2019 KSA in disucussion with UDEL group
for iTT=1:nTT
                            p2Greens(:,:) = shiftdim(p2(izr,izs,ishipspd,ix0,:,:)) ; %shift dimensions to a vector for p2(izr,izs,ishipspd,ix0,iTT,ifreq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            
          %                  p2Greens(:,1:if1-1) = min(p2Greens(:))*1e-12*ones(size(p2Greens(:,1:if1-1))); %Define Greeens to be very small for frequencies < fmin
                        end
                        psrc=PSource(izs,:); % complexx Source Pressure at depth # izs, for frequencies: freq
                        ppp= (ones(nTT,1)*psrc).*p2Greens;  %pressure at depth ppp( = Source pressure * pressure greens function
                   
                    %    ppp=ppp +10.^(65/20);
                        ppparray(isamp,:,:)=p2Greens(:,:); % array for each isamp 
%                         p2pvtemp = shiftdim(p2pv(izr,izs,ishipspd,ix0,:,:)) ;  %shift dimensions to a vector
%                         ppppv= (ones(nTT,1)*psrc).*p2pvtemp;  % %PV at depth = Source pressure * PV greens function
                       TL=20*log10(abs(p2Greens));
                   %    ppp_fullfield=PPP(:,:,4:8);
                  % ppp_fullfield=ppp(:,:,1:nfreq);
                   ppp_fullfield=ppparray(:,:,1:nfreq);

                     
                    end   %for ishipspd=1:length(sourceKnotsa) %source velocity
                end % for ix0=1:length(x0a) x0=x0a(ix0);  % simulate varying cpa Ranges x0
            end    % zs=zsa(izs); source depth loop in m
        end  %zr=zra(izr); receiver depth loop
        %
    end %  SSclass loop % 1:21
end %  geoclass loop % 1:5

iicount=0;
              for igclass=1:ngeoa
               for iSSclass=1:nSSa
                  for izs=1:nzsa
                       for ix0=1:nx0
                           for ishipspd=1:nsh
                            for izr=1:nzra
                                 iicount=iicount+1;
                                     Label(igclass,iSSclass,ix0,ishipspd,izs,izr)=iicount;
                                                         end
                                                 end
                                   end
                          end
               end 
          end
                               

%%
TL=-TL; % makes transmission loss positive

SL=mean(transpose(PSDselectf) + TL);
%SL=SL_p;
%SL=230-10*log10(fselect.^3.594)+10*log10((1+(fselect/340).^2).^0.917);
%SL=transpose(SL);
%%
%for i =1:550
 %   SL(i)=170;
%end
%%
%SL=transpose(SL);
TL_mea=transpose(PSDselectf)-SL;
RL_mod=SL+TL;
clims=[60 120];
imagesc(fselect,tspec,RL_mod,clims)
colorbar
S=230-10*log10(fselect.^3.594)+10*log10((1+(fselect/340).^2).^0.917);

%%
clear fr fprime V MM C Cinv Vinv Beta betatop betabottom AA AAA
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
 
         for ii=1:nmodes
             for jj=1:nmodes
         betatop(ii,jj,:)=Cinv(ii,:)-Cinv(jj,:);
         betabottom(ii,jj,:)=Vinv(ii,:)-Vinv(jj,:);
             end
         end

 Beta=betatop ./betabottom;
 AA(:)=-Beta(1,2,:);
AA=AA ./transpose(fselect);
[pks,locs]=findpeaks(AA);
smallPeakIndexes = pks < 0.1*max(pks);
pks(smallPeakIndexes) = [] ; %Reject Y value of peaks below this threshold
locs(smallPeakIndexes) = [] ;
F1=fselect(locs(length(locs)))

 AA(:)=-Beta(1,3,:);
AA=AA ./transpose(fselect);
[pks,locs]=findpeaks(AA);
smallPeakIndexes = pks < 0.1*max(pks);
pks(smallPeakIndexes) = [] ; %Reject Y value of peaks below this threshold
locs(smallPeakIndexes) = [] ;
F2=fselect(locs(length(locs)))


 AA(:)=-Beta(1,4,:);
AA=AA ./transpose(fselect);
[pks,locs]=findpeaks(AA);
smallPeakIndexes = pks < 0.1*max(pks);
pks(smallPeakIndexes) = [] ; %Reject Y value of peaks below this threshold
locs(smallPeakIndexes) = [] ;
F3=fselect(locs(length(locs)))

 AA(:)=-Beta(1,5,:);
AA=AA ./transpose(fselect);
[pks,locs]=findpeaks(AA);
smallPeakIndexes = pks < 0.1*max(pks);
pks(smallPeakIndexes) = [] ; %Reject Y value of peaks below this threshold
locs(smallPeakIndexes) = [] ;
F4=fselect(locs(length(locs)))




%%
clear sq1 sq2 sq3 sq4 SUM

%    Error=(transpose(PSDselectf)-(SL+TL)).^2;
 %   SQERROR=(1.0/length(tspec))*(1.0/length(fselect))*sum(sum(Error));
  %  Errormatrix(NNN)=sqrt(SQERROR);
  

  % sq1=0.0
 % sq2=(F2 - 40.5)^2; %VB vla2
 %   sq3=(F3 - 55.245)^2; %VB vla2
%sq4=0;



 %   sq1=(F1 - 24.98)^2;  %Corriedo vla2
 %   sq2=(F2 - 41.28)^2; %Corriedo vla2
  %  sq3=(F3 - 53.72)^2; %Corriedo vla2
  % sq4=0;
   
       sq1=(F1 - 25)^2;  %KALAMATA vla2
    sq2=(F2 - 40.7)^2; %KALAMATA vla2
    sq3=(F3 - 52.78)^2; %KALAMATA vla2
   sq4=0;

   
% sq1=(F1 - 24.99)^2; %matsuyama vla2
 %  sq2=(F2 - 40.72)^2; %matsuyama vla2
 %  sq3=(F3 - 53.72)^2; %matsuyama vla2
 % sq4=(F4 - 69.86)^2; %matsuyama vla2
 
 
 % sq1=(F1 - 25.3)^2; %hafnia green vla2
  % sq2=(F2 - 41.78)^2; %hafnia green vla2
  % sq3=(F3 - 53.83)^2; %hafnia green vla2
 % sq4=(F4 - 68.43)^2; %hafnia green vla2
 %sq4=0;


    

 
     SUM=sq1+sq2+sq3+sq4
     a11=Parm_1(NNN)
     b11=Parm_2(NNN)
     c11=Parm_3(NNN)
     Errormatrix(NNN)=SUM
     
        end

 dist_file=[transpose(Errormatrix) Parm_1 Parm_2 Parm_3 ];
elapsedTime = toc %#ok<NOPTS>
disp('Program has completed.') 
toc
