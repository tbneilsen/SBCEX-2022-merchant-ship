% ==============================================================
% COPYRIGHT (c) (2022)
% Applied Research Laboratories, The University of Texas at Austin (ARL:UT)
%
% This material may be reproduced by or for the U.S. Government pursuant
% to the copyright license under the clause at DFARS 252.227-7013. This
% notice must appear in all copies of this material and its derivatives.
%
% FOR OFFICIAL USE ONLY (exemption 4)
%
% This material is the intellectual property of the University of Texas at Austin
% Any unauthorized reproduction or use is prohibited by law.
%
% Distribution Statement D. Distribution authorized to the Department of Defense
% and U.S. DoD contractors only (critical information) (9 January 2015)
% Other requests shall be referred to (ONR Code 322US).
%
% THIS SOFTWARE IS PROVIDED  "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
% EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NONINFRINGEMENT.
%
% IN NO EVENT WILL ARL:UT BE LIABLE  FOR ANY DAMAGES WHATSOEVER ARISING OUT
% OF THE USE OR INABILITY TO USE THIS SOFTWARE, EVEN IF ARL:UT HAS BEEN
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
% ============================================================== 
% Jason D. Sagers, PhD Revised 9/21/2022

function data = PROTEUS_ReadData_Function(cdf_filedir,cdf_filename)
%
% Inputs 
% cdf_filedir   file directory for cdf files
% cdf_filename  file name of cdf
%
% Outputs
% data          .hyd_data   (hydrophone data in array channel order
%               .t_unixsec  (unix timestamp for data
%
%% Reads the cdf data and provides output_cdf structure
% Preallocate
cdf_output_ac   = [];
cdf_output_cnt  = [];
cdf_output_mdw1 = [];
cdf_output_mdw2 = [];

Nread = 1e6;
cdf_filename = fullfile(cdf_filedir,cdf_filename);
    
% Open CDF file
cdf_fid = fopen(cdf_filename);

% Read all messages in the file for the acoustic channels
[cdf_output_ac]     = [cdf_output_ac, cdf_fixed_reader(cdf_fid,Nread,[167,1])]; frewind(cdf_fid);
[cdf_output_cnt]    = [cdf_output_cnt, cdf_fixed_reader(cdf_fid,Nread,[167,2])];frewind(cdf_fid);
[cdf_output_mdw1]   = [cdf_output_mdw1, cdf_fixed_reader(cdf_fid,Nread,[0,3])];frewind(cdf_fid);
[cdf_output_mdw2]   = [cdf_output_mdw2, cdf_fixed_reader(cdf_fid,Nread,[0,4])];frewind(cdf_fid);

% Clean up
fclose(cdf_fid);
    
%% Reshape bus data 
num_messages        = length(cdf_output_ac);
nsamp_per_chan      = cdf_output_ac(1).fixed_hdr.base_sample_count;
total_samp          = num_messages*nsamp_per_chan;
bus_data            = zeros(total_samp,80);
t_unixsec_record    = zeros(num_messages,1);

for imnum = 1:length(cdf_output_ac)
    for ibus = 1:8
        bus_record = cdf_output_ac(imnum).fixed_aperture(ibus).data;
        bus_record = reshape(bus_record,8,length(bus_record)/8);
        
        % Scale the data from bits to volts, but don't adjust for ANY gains
        bus_record = bus_record*(cdf_output_mdw1(1).middleware.aperture(ibus).adc_range_volts/(2^cdf_output_mdw1(1).middleware.aperture(ibus).adc_bit_count));
        
        % Get the counts
        cnt_record = cdf_output_cnt(imnum).fixed_aperture(ibus).data;
        cnt_record = reshape(cnt_record,2,length(cnt_record)/2);
        
        % Stuff 8 buses + 2 buses into a PRADIOS record
        record = [bus_record; cnt_record];
        array_record(:,(ibus-1)*10+1:ibus*10) = record.';
        
        % Extract the record timestamp
        array_record_unixsec = cdf_output_ac(imnum).fixed_hdr.timestamp;
    end
        bus_data((imnum-1)*nsamp_per_chan+1:imnum*nsamp_per_chan,:) = array_record;
        t_unixsec_record(imnum) = array_record_unixsec;
end

%% Extract the gains

for ibus = 1:8;
    ESL_board_no        = cdf_output_mdw1(1).middleware.aperture(ibus).fe_board_pn;
    raw_setting         = cdf_output_mdw1(1).middleware.aperture(ibus).gain_settings;
    [gain,sensitivity]  = ESL_interpret_gain_sensitivity_setting( ESL_board_no, raw_setting );
    bus(ibus).gain_pg   = gain.pg;
    bus(ibus).gain_fg   = gain.fg;
    bus(ibus).sensitivity_index = sensitivity;
    if bus(ibus).sensitivity_index == 0
        bus(ibus).sensitivity = cdf_output_mdw2(1).middleware.mapping(1).sensitivity_0;
    else
        bus(ibus).sensitivity = cdf_output_mdw2(1).middleware.mapping(1).sensitivity_1;
    end
end

% Construct a gain and sensitivity vector for each channel
bsi = [1 11 21 31 41 51 61 71]-1; % Starting index for each bus
PG  = ones(1,80);
FG  = ones(1,80);
S   = zeros(1,80);
for ibus = 1:8
    PG(bsi(ibus)+[1:8]) = bus(ibus).gain_pg;
    FG(bsi(ibus)+[1:8]) = bus(ibus).gain_fg;
    S(bsi(ibus)+[1:8])  = bus(ibus).sensitivity;
end

%% Apply the programmable gains
% Adjust voltage for programmable gain
bus_data = bus_data.*(repmat(PG,size(bus_data,1),1));

% Extract counters into their own matrix
data.counter = bus_data(:,[9:10 19:20 29:30 39:40 49:50 59:60 69:70 79:80]);

% Extract the acoustic data into their own matrix 
bus_data = bus_data(:,[1:8 11:18 21:28 31:38 41:48 51:58 61:68 71:78]);
S = S([1:8 11:18 21:28 31:38 41:48 51:58 61:68 71:78]);

%% Apply the fixed gains from the measured transfer function 

nfft    = size(bus_data,1);
fs      = cdf_output_ac(1).fixed_hdr.base_sample_rate;
dt      = 1/fs;
df      = fs/nfft;
f       = 0:df:fs/2-df;

% Convert the stored FRF matrix into the right format for use here

load('PFE_MN1902_RevD_Stack1_2019.04.10.mat');
icount = 0;
for ibus = 1:8
    for ichan = 1:8
        icount = icount + 1;
        frf(icount,:) = PFES.frf{ibus,ichan};
    end
end

% Condition the frf to not extend below a few Hz.  Otherwise, we will
% amplify low frequency data which has already been attenuated.  The HTI90
% phones have a cutoff of 2 Hz
[~,fidx] = min(abs(PFES.f-2));
frf(:,1:fidx-1) = repmat(frf(:,fidx),1,fidx-1);

% Interpret the frf structure to the frequency vector used in the FFT
% this will account for the fixed 5x gain in the passband and the frequency response 
cal_QN = interp1(PFES.f,frf(:,:).',f,'linear',0); 
cal_QN(isnan(cal_QN) | cal_QN==0)=1/eps;

% Convert the hydrophone sensitivity from dB into lineary units
hydro_1N = 10.^(-S/20);

% Add the hydrophone sensitivity to the cal_QN
cal_QN = (bsxfun(@rdivide,cal_QN,hydro_1N));

% FFT the data so I can apply the sensitivity-dependent fixed gain
fft_data_qN = fft(bus_data,nfft,1);

% IFFT the result to produce a calibrated time series (in uPa)
data.hyd_data = ifft(fft_data_qN(1:nfft/2,:)./cal_QN,nfft,1,'symmetric').';

%% Reorder the data based on array channel order
for i = 1:52
    array_order(i)  = cdf_output_mdw2(1).middleware.mapping(i).crio_channel;
    xelem(i)        = cdf_output_mdw2(1).middleware.mapping(i).location_meters;
end

data.hyd_data   = data.hyd_data(array_order,:);
data.fs         = fs;
data.xelem      = xelem;
data.prog_gain  = PG;
data.hyd_sens   = S;

%% Make time series vector
samp_per_record = cdf_output_ac(1).fixed_hdr.base_sample_count;
t_unixsec = [];
for i = 1:length(t_unixsec_record)
    t_unixsec = [t_unixsec t_unixsec_record(i)+[0:samp_per_record-1]*dt];
end

data.t_unixsec = t_unixsec;

%% Diagnostic plots
plot_flag = 0;
if plot_flag == 1
    for ichan = 1;
        % ichan = 50;
        nfft = round(fs/100);
        [~,f,t,p] = spectrogram(data.hyd_data(ichan,:),hann(nfft),nfft/2,nfft,fs);
        figure(10)
        clf;
        imagesc(t,f,10*log10(p));colormap jet; colorbar;set(gca,'ydir','normal');
        caxis([40 150]); title(ichan)
        drawnow
        pause
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output] = cdf_fixed_reader(fid, num_messages, sensors)
%CDF_FIXED_READER Decodes one or more CDF message payload contained within
%                      a data stream.
%
% USAGE:  [output_struct] = ...
%          CDF_FIXED_READER(fid, num_messages, sensors)201-365

% fid                     The FID to a raw CDF stream recorded as a file.
%
% num_message             Number of messages to return that match filter.
%                         This is total messages that match num_messages
%                         (i.e. combined count of sensor pairs requested).
%                         Non-matching sensor pairs will be skipped and
%                         not count as an output message incrementcdf_fixed.
%                         Set to 0 if all available file messages desired.
%
% sensors                 OPTIONAL: 2D array of [source, id] pairs.
%                                   Form is [source, id; source, id; ...]
%                                   Default is [167,1] for fixed ID 1.
%
% OUTPUTS:
% output                  The parsed messages as struct.
%                         Empty if non matching content.
%
% EXAMPLE:
% fid      = fopen('cdf_test_data.cdf');
% cdf_msgs = cdf_fixed_reader(fid, 10, [167,1; 167,2; 0,0])
%


% check input arguments for options.
if nargin == 1
    num_messages = -1;
end

if nargin <= 2
    sensors = [167,1];
    data_sources = [167];
    message_ids = [1];
else
    data_sources = sensors(:,1);
    message_ids = sensors(:,2);
end

mnum = 0;
while mnum ~= num_messages && ~feof(fid)
    % increment the message counter.
    if num_messages ~= 0
        mnum = mnum + 1;
    end
    
    % read the CDF header.
    output(mnum).cdf_hdr.sync_word_start = ...
        fread(fid, 1, 'uint32', 0, 'ieee-be');
    output(mnum).cdf_hdr.time_tag_sec    = ...
        fread(fid, 1, 'uint32', 0, 'ieee-be');
    output(mnum).cdf_hdr.time_tag_nsec   = ...
        fread(fid, 1, 'uint32', 0, 'ieee-be');
    output(mnum).cdf_hdr.data_source     = ...
        fread(fid, 1, 'uint32', 0, 'ieee-be');
    output(mnum).cdf_hdr.message_id      = ...
        fread(fid, 1, 'uint32', 0, 'ieee-be');
    output(mnum).cdf_hdr.byte_count      = ...
        fread(fid, 1, 'uint32', 0, 'ieee-be');
    output(mnum).cdf_hdr.sync_word_stop  = ...
        fread(fid, 1, 'uint32', 0, 'ieee-be');
    
    % because sensors have to be treated as source-id pairs, it is
    % necessary to have logic that doesn't cross-up a unintended pairing.
    % that is to say, only message ids that that match a source should be
    % valid and not just message ids that match any sources.  this is why
    % the message_ids are sub-selected from only pairs that have the
    % matching source in the following conditionals.
    
    % read in the CDF payload.
    if ~isempty(find(output(mnum).cdf_hdr.data_source == 167)) && ...
            ~isempty(find(data_sources == output(mnum).cdf_hdr.data_source)) && ...
            ~isempty(find(message_ids(find(data_sources == 167)) == ...
            output(mnum).cdf_hdr.message_id))
        % fixed acoustic message (source 167, id <any>).
        
        % get the fixed header info.
        output(mnum).fixed_hdr.timestamp         = ...
            fread(fid, 1, 'double', 0, 'ieee-be');
        output(mnum).fixed_hdr.msg_source        = ...
            fread(fid, 1, 'uint32', 0, 'ieee-be');
        output(mnum).fixed_hdr.base_sample_rate  = ...
            fread(fid, 1, 'single', 0, 'ieee-be');
        output(mnum).fixed_hdr.base_sample_count = ...
            fread(fid, 1, 'uint32', 0, 'ieee-be');
        output(mnum).fixed_hdr.channel_data_type = ...
            fread(fid, 1, 'uint32', 0, 'ieee-be');
        output(mnum).fixed_hdr.channel_data_size = ...
            fread(fid, 1, 'uint32', 0, 'ieee-be');
        output(mnum).fixed_hdr.aperture_count    = ...
            fread(fid, 1, 'uint32', 0, 'ieee-be');
        
        % get the data type of the aperture(s).
        switch (output(mnum).fixed_hdr.channel_data_type)
            case 0
                % signed integer (2's compliment)
                switch (output(mnum).fixed_hdr.channel_data_size)
                    case 2,
                        fixed_data_type = 'int16';
                    case 4,
                        fixed_data_type = 'int32';
                    otherwise
                        disp('Data type cannot be decoded.');
                        return;
                end
                
            case 1
                % unsigned integer
                switch (output(mnum).fixed_hdr.channel_data_size)
                    case 2,
                        fixed_data_type = 'uint16';
                    case 4,
                        fixed_data_type = 'uint32';
                    otherwise
                        disp('Data type cannot be decoded.');
                        return;
                end
            case 3
                % float
                fixed_data_type = 'single';
            otherwise
                fprintf('Data type cannot be decoded.\n');
                return;
        end
        
        % get the aperture data.
        for aperture_index = [1:output(mnum).fixed_hdr.aperture_count]
            output(mnum).fixed_aperture(aperture_index).sample_rate_factor = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            output(mnum).fixed_aperture(aperture_index).channel_count      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            output(mnum).fixed_aperture(aperture_index).data = ...
                fread(fid, ...
                output(mnum).fixed_aperture(aperture_index).channel_count * ...
                output(mnum).fixed_hdr.base_sample_count, ...
                fixed_data_type, 0, 'ieee-be');
        end
    elseif ~isempty(find(output(mnum).cdf_hdr.data_source == 0)) && ...
            ~isempty(find(message_ids(find(data_sources == 0)) == ...
            output(mnum).cdf_hdr.message_id)) && ...
            ~isempty(find(output(mnum).cdf_hdr.message_id == 3))
        % middleware message (source 0, id 3).
        
        % get the PRADIOS software version.
        output(mnum).middleware.pradios_sw_ver = ...
            fread(fid, 1, 'uint16', 0, 'ieee-be');
        
        % get the number of aperture busses.
        output(mnum).middleware.number_of_aperture_bus = ...
            fread(fid, 1, 'uint32', 0, 'ieee-be');
        
        % get the aperture non-acoustic data.
        for bus_indx = 1:output(mnum).middleware.number_of_aperture_bus
            
            % parse data depending on pradios version.
            % 	    if (output(mnum).middleware.pradios_sw_ver >= 3) && ...
            %               (output(mnum).middleware.pradios_sw_ver < 1000)
            if (output(mnum).middleware.pradios_sw_ver == 3)
                
                output(mnum).middleware.aperture(bus_indx).bus_number        = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).adios_hdr_size    = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).record_size_bytes = ...
                    fread(fid, 1, 'int32' , 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).file_number       = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).record_number     = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).sample_rate_hz    = ...
                    fread(fid, 1, 'int32' , 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).fe_board_pn       = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).gain_settings     = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).adc_range_volts   = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).adc_bit_count     = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).unix_time         = ...
                    fread(fid, 1, 'int32' , 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).adcs_per_card     = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                % Enum: 0-None, 1-System, 2-External
                output(mnum).middleware.aperture(bus_indx).time_source       = ...
                    fread(fid, 1, 'uint8', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).subsecond_ns      = ...
                    fread(fid, 1, 'uint32', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).offset_from_time_ref = ...
                    fread(fid, 1, 'int32', 0, 'ieee-be');
                % true/false
                output(mnum).middleware.aperture(bus_indx).offset_from_time_ref_bool = ...
                    fread(fid, 1, 'uint8', 0, 'ieee-be');
                % true/false
                output(mnum).middleware.aperture(bus_indx).time_sync_fault_bool = ...
                    fread(fid, 1, 'uint8', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).sync_pattern      = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                
            elseif(output(mnum).middleware.pradios_sw_ver == 4)
                
                output(mnum).middleware.aperture(bus_indx).bus_number        = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).adios_hdr_size    = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).record_size_bytes = ...
                    fread(fid, 1, 'int32' , 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).file_number       = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).record_number     = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).adc_range_volts   = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).adc_bit_count     = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).unix_time         = ...
                    fread(fid, 1, 'int32' , 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).adcs_per_card     = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                % Enum: 0-None, 1-System, 2-External
                output(mnum).middleware.aperture(bus_indx).time_source       = ...
                    fread(fid, 1, 'uint8', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).subsecond_ns      = ...
                    fread(fid, 1, 'uint32', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).offset_from_time_ref = ...
                    fread(fid, 1, 'int32', 0, 'ieee-be');
                % true/false
                output(mnum).middleware.aperture(bus_indx).offset_from_time_ref_bool = ...
                    fread(fid, 1, 'uint8', 0, 'ieee-be');
                % true/false
                output(mnum).middleware.aperture(bus_indx).time_sync_fault_bool = ...
                    fread(fid, 1, 'uint8', 0, 'ieee-be');
                output(mnum).middleware.aperture(bus_indx).sync_pattern      = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
            else
                % no-op (until more versions filled-in).
            end
        end
        
        if (output(mnum).middleware.pradios_sw_ver == 4)
            
            % Read EEPROM Values
            output(mnum).middleware.eeprom.modus_vn      = ...
                fread(fid, 1, 'double', 0, 'ieee-be');
            
            % Read Battery Housing info
            for batt_house_indx = 1:fread(fid, 1, 'int32', 0, 'ieee-be')
                output(mnum).middleware.eeprom.batt_house(batt_house_indx).pn      = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.eeprom.batt_house(batt_house_indx).sn      = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.eeprom.batt_house(batt_house_indx).voltage      = ...
                    fread(fid, 1, 'double', 0, 'ieee-be');
            end
            
            % Read Watchdawg info
            % Enum: 0- 16s, 1- 4m 15s
            output(mnum).middleware.eeprom.wd_timout      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.wd_cyc_error      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            output(mnum).middleware.eeprom.wd_sys_error      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            
            % Read AFE info
            for afe_indx = 1:fread(fid, 1, 'int32', 0, 'ieee-be')
                output(mnum).middleware.eeprom.afe(afe_indx).sn      = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                output(mnum).middleware.eeprom.afe(afe_indx).pn      = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                % Enum: 0-21, each value representing a unique gain setting.
                % See PGA281 data sheet for exact values
                output(mnum).middleware.eeprom.afe(afe_indx).gain      = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
                % Enum: 0- 12V OFF, 1- 12V ON
                output(mnum).middleware.eeprom.afe(afe_indx).v12      = ...
                    fread(fid, 1, 'uint16', 0, 'ieee-be');
            end
            
            % Read Supervisor Info
            output(mnum).middleware.eeprom.supervisor.pn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.supervisor.sn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.supervisor.secs_in_a_day      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            output(mnum).middleware.eeprom.supervisor.rec_boot      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.supervisor.rec_ceanup      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.supervisor.sch_timer      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            output(mnum).middleware.eeprom.supervisor.cont_timer      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            % Read in the schedule string
            for sch_indx = 1:fread(fid, 1, 'int32', 0, 'ieee-be')
                output(mnum).middleware.eeprom.supervisor.sch(sch_indx)      = ...
                    fread(fid, 1, 'char*1', 0, 'ieee-be');
            end
            output(mnum).middleware.eeprom.supervisor.last_eeprom_flash      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            output(mnum).middleware.eeprom.supervisor.debug      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.supervisor.fail_test      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.supervisor.verbose      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            
            % Read System Info
            output(mnum).middleware.eeprom.system.modus_sn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.system.eeprom_vn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.system.num_afe_boards      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.system.num_battery_housings      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            
            % Read Battery Interface Info
            output(mnum).middleware.eeprom.battery_interface.pn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.battery_interface.sn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.battery_interface.bench_exc_pwr      = ...
                fread(fid, 1, 'double', 0, 'ieee-be');
            output(mnum).middleware.eeprom.battery_interface.idle_exc_pwr      = ...
                fread(fid, 1, 'double', 0, 'ieee-be');
            output(mnum).middleware.eeprom.battery_interface.record_exc_pwr      = ...
                fread(fid, 1, 'double', 0, 'ieee-be');
            
            % Read Ethernet Switch Info
            output(mnum).middleware.eeprom.ethernet_switch.pn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            
            % Read PTP Info
            output(mnum).middleware.eeprom.ptp.pn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            
            % Read Clock Board Info
            output(mnum).middleware.eeprom.clk_brd.pn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.clk_brd.sn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.clk_brd.last_time_sync      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            
            % Read Backup Battery Info
            output(mnum).middleware.eeprom.bu_batt.starting_voltage      = ...
                fread(fid, 1, 'double', 0, 'ieee-be');
            output(mnum).middleware.eeprom.bu_batt.nominal_voltage      = ...
                fread(fid, 1, 'double', 0, 'ieee-be');
            output(mnum).middleware.eeprom.bu_batt.insertion_date      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            output(mnum).middleware.eeprom.bu_batt.last_charge_date      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            
            % Read SIC Board Info
            output(mnum).middleware.eeprom.sic_brd.pn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.sic_brd.sn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            
            % Read ACOMMs Info
            output(mnum).middleware.eeprom.acomms.freq      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            output(mnum).middleware.eeprom.acomms.battery_start_voltage      = ...
                fread(fid, 1, 'double', 0, 'ieee-be');
            output(mnum).middleware.eeprom.acomms.battery_nominal_voltage      = ...
                fread(fid, 1, 'double', 0, 'ieee-be');
            output(mnum).middleware.eeprom.acomms.battery_insert_date      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            output(mnum).middleware.eeprom.acomms.xdcr_pn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.acomms.modem_pn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            
            % Read Recorder Info
            output(mnum).middleware.eeprom.crio.num_modules      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.crio.rad_image_num      = ...
                fread(fid, 1, 'double', 0, 'ieee-be');
            % Enum: check crio for values
            output(mnum).middleware.eeprom.crio.data_rate      = ...
                fread(fid, 1, 'uint8', 0, 'ieee-be');
            output(mnum).middleware.eeprom.crio.pradios_sw_version      = ...
                fread(fid, 1, 'double', 0, 'ieee-be');
            % Enum: 0- uSD card
            output(mnum).middleware.eeprom.crio.memory_type      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            % Enum: 0- cRIO-9045
            output(mnum).middleware.eeprom.crio.pn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.crio.autorun_timer      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            
            % Read Power Board Info
            output(mnum).middleware.eeprom.pwr_brd.pn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.pwr_brd.sn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.pwr_brd.energy      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            output(mnum).middleware.eeprom.pwr_brd.charge      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            output(mnum).middleware.eeprom.pwr_brd.time      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
            
            % Read Ext Pressure Sensor Info
            output(mnum).middleware.eeprom.ext_press_sens.sn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            % Enum: 0- GT2100
            output(mnum).middleware.eeprom.ext_press_sens.pn      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.ext_press_sens.samps_to_average      = ...
                fread(fid, 1, 'uint16', 0, 'ieee-be');
            output(mnum).middleware.eeprom.ext_press_sens.volts_to_PSI      = ...
                fread(fid, 1, 'double', 0, 'ieee-be');
            output(mnum).middleware.eeprom.ext_press_sens.zero_offset      = ...
                fread(fid, 1, 'double', 0, 'ieee-be');
            output(mnum).middleware.eeprom.ext_press_sens.sens_settle_time      = ...
                fread(fid, 1, 'uint32', 0, 'ieee-be');
        end
        
        
        
    elseif ~isempty(find(output(mnum).cdf_hdr.data_source == 0)) && ...
            ~isempty(find(message_ids(find(data_sources == 0)) == ...
            output(mnum).cdf_hdr.message_id)) && ...
            ~isempty(find(output(mnum).cdf_hdr.message_id == 4))
        % middleware message (source 0, id 4).
        
        % get the PRADIOS software version.
        output(mnum).middleware.pradios_sw_ver = ...
            fread(fid, 1, 'uint16', 0, 'ieee-be');
        
        % get the number of fields and phones.
        output(mnum).middleware.number_of_phones = ...
            fread(fid, 1, 'uint32', 0, 'ieee-be');
        output(mnum).middleware.number_of_fields = ...
            fread(fid, 1, 'uint32', 0, 'ieee-be');
        
        % get the non-acoustic data channel mapping.
        for phone_indx = 1:output(mnum).middleware.number_of_phones
            output(mnum).middleware.mapping(phone_indx).hyd_num         = ...
                fread(fid, 1, 'single', 0, 'ieee-be');
            output(mnum).middleware.mapping(phone_indx).hydrophone      = ...
                fread(fid, 1, 'single', 0, 'ieee-be');
            output(mnum).middleware.mapping(phone_indx).location_meters = ...
                fread(fid, 1, 'single', 0, 'ieee-be');
            output(mnum).middleware.mapping(phone_indx).crio_channel    = ...
                fread(fid, 1, 'single', 0, 'ieee-be');
            output(mnum).middleware.mapping(phone_indx).sensitivity_0   = ...
                fread(fid, 1, 'single', 0, 'ieee-be');
            output(mnum).middleware.mapping(phone_indx).sensitivity_1   = ...
                fread(fid, 1, 'single', 0, 'ieee-be');
            output(mnum).middleware.mapping(phone_indx).serial_number   = ...
                fread(fid, 1, 'single', 0, 'ieee-be');
            
            % not encoded, but note associated bus number is the value
            % floor(crio_channel/8) + 1 when using 8 busses and 1-based
            % indexing to represent bus number.
        end
        
    else
        if feof(fid)
            % reset the mnum back one to recycle the message index slot.
            mnum = mnum - 1;
            
            % end of file.
            break;
        end
        
        % skip over this content.
        fseek(fid, output(mnum).cdf_hdr.byte_count, 'cof');
        
        % reset the mnum back one to recycle the message index slot.
        mnum = mnum - 1;
    end
end

% in case we end up with a residual message at the end or empty, clean up.
if mnum == 0
    output = [];
else
    output = output(1:mnum);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gain,sensitivity] = ESL_interpret_gain_sensitivity_setting( ESL_board_no, raw_setting )
% [gain,sensitivity] = ESL_interpret_gain_sensitivity_setting( ESL_board_no, raw_setting )
%
% Interprets the gain and sensitivity settings set based on the ESL board
% part number
%
% Takes 2 arguments:
%
%   ESL_board_no - the board part these settings are for (logged in PRADIOS
%   header
%   raw_settings - the setting logged in the PRADIOS header to be
%   interpreted
%
% Returns 2 values:
%
%   gain - The gain in V/V (including fixed gain) that has been applied to the
%   signal
%   sensitivity  - A value that has meaning based on the board (board 1804
%   -> 1 means the 12V was sent to hydrophone, 0 means whatever was
%   supplied to the board was sent out, typically 24V)
%


switch ESL_board_no
    
    %% This section deals with the v1.0 analog front end
    case {1804, 1902}
        if(raw_setting ~= hex2dec('ff'))    %check the validity of the setting
            gain_setting = bitand(raw_setting,hex2dec('001f'));         %the gain settings are on the lower 5 bits
            sensitivity_setting = bitand(raw_setting,hex2dec('0020'));         %the sensitivity setting is on the 6th bit
            
            %This sets 1 if the 12V was sent out to the board, or 0 if it was
            %just 24/12V passthough
            if(sensitivity_setting)
                sensitivity = 1;
            else
                sensitivity = 0;
            end
            
            %The PGA uses the 5th bit to add weird settings
            if( bitand(gain_setting,hex2dec('0010')))
                switch (bitand(gain_setting,hex2dec('000f')))
                    case 0
                        programmable_gain = .172;
                    case 1
                        programmable_gain = .344;
                    case 2
                        programmable_gain = .688;
                    case 3
                        programmable_gain = 1.375;
                    case 4
                        programmable_gain = 2.75;
                    case 5
                        programmable_gain = 5.5;
                    case 6
                        programmable_gain = 11;
                    case 7
                        programmable_gain = 22;
                    case 8
                        programmable_gain = 44;
                    case 9
                        programmable_gain = 88;
                    case 10
                        programmable_gain = 176;
                    case 11
                        programmable_gain = .172;
                    case 12
                        programmable_gain = .172;
                    case 13
                        programmable_gain = .172;
                    case 14
                        programmable_gain = .172;
                    case 15
                        programmable_gain = .172;
                end
            else
                switch (bitand(gain_setting,hex2dec('000f')))
                    case 0
                        programmable_gain = .125;
                    case 1
                        programmable_gain = .25;
                    case 2
                        programmable_gain = .5;
                    case 3
                        programmable_gain = 1;
                    case 4
                        programmable_gain = 2;
                    case 5
                        programmable_gain = 4;
                    case 6
                        programmable_gain = 8;
                    case 7
                        programmable_gain = 16;
                    case 8
                        programmable_gain = 32;
                    case 9
                        programmable_gain = 64;
                    case 10
                        programmable_gain = 128;
                    case 11
                        programmable_gain = .125;
                    case 12
                        programmable_gain = .125;
                    case 13
                        programmable_gain = .125;
                    case 14
                        programmable_gain = .125;
                    case 15
                        programmable_gain = .125;
                end
            end
            
            %%Finally, there is a fixed gain of 5 on the board itself
            gain.pg = programmable_gain;
            gain.fg = 5;
        else
            gain.pg = nan;
            gain.fg = nan;
            sensitivity = nan;
        end
        
    case {2002}
        if(raw_setting ~= hex2dec('ff'))    %check the validity of the setting
            gain_setting = bitand(raw_setting,hex2dec('001f'));         %the gain settings are on the lower 5 bits
            sensitivity_setting = bitand(raw_setting,hex2dec('0020'));         %the sensitivity setting is on the 6th bit
            
            %This sets 1 if the 12V was sent out to the board, or 0 if it was
            %just 24/12V passthough
            if(sensitivity_setting)
                sensitivity = 1;
            else
                sensitivity = 0;
            end
            
            switch (gain_setting)
                case 0
                    programmable_gain = .125;
                case 1
                    programmable_gain = .172;
                case 2
                    programmable_gain = .25;
                case 3
                    programmable_gain = .344;
                case 4
                    programmable_gain = .5;
                case 5
                    programmable_gain = .688;
                case 6
                    programmable_gain = 1;
                case 7
                    programmable_gain = 1.375;
                case 8
                    programmable_gain = 2;
                case 9
                    programmable_gain = 2.75;
                case 10
                    programmable_gain = 4;
                case 11
                    programmable_gain = 5.5;
                case 12
                    programmable_gain = 8;
                case 13
                    programmable_gain = 11;
                case 14
                    programmable_gain = 16;
                case 15
                    programmable_gain = 22;
                case 16
                    programmable_gain = 32;
                case 17
                    programmable_gain = 44;
                case 18
                    programmable_gain = 64;
                case 19
                    programmable_gain = 88;
                case 20
                    programmable_gain = 128;
                case 21
                    programmable_gain = 176;
            end
            
            %%Finally, there is a fixed gain of 5 on the board itself
            gain.pg = programmable_gain;
            gain.fg = 5;
        else
            gain.pg = nan;
            gain.fg = nan;
            sensitivity = nan;
        end
        
        %% This is all other boards
    otherwise
               
end
end