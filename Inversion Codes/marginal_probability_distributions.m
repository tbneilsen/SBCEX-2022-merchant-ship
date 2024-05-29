   
% /****************************************************************
%  *                                                              *
%  *     COPYRIGHT © 2016                                         *
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
% Inputs ------------------------------------------------------------------
%
%   dist_file   - file path/name of _dist file. 
%
%   N           - number of beta=infinity points used as a feature for the
%                 maximum entropy inversion
%
%   Ndiv        - the number of bins for averaging (default is 50 bins)
%
% Outputs ------------------------------------------------------------------
%
%   dist_file   - original but with added information about the ship
%
%   dist_plot   - plot of the marginal probability distributions for the
%                 four parameters
%
% Plots probability distributions of parameters contained in a _dist file
% On the output plot, inversion parameters are displayed such that the areas 
% under each curve look the same regardless of x axis values
% On the output plot, source level parameters are displayed such that the
% areas under each curve are exactly 1, using/including the x axis values.
%
% Each of the distributions has labeled the mean, standard deviation of the
% mean, optimal (based on the lowest cost function), and the peak (based
% on the marginal probability distributions) of the four parameter values.

% This code was orginally written by Sumedh Joshi (deceased) and David
% Knobles (1956-present).  It computes marginal distributions from a PPD 
% computed by the Maximum Entropy principle as applied in ocean acoustics,  
% D. P. Knobles et al. Maximum Entropy  approach for statistical inference 
% in an ocean acoustic waveguid," Journal Acoustical Society of America 131  
% 1087-1101 (2012)

% Modified: Alexandra M. Hopps-McDaniel (5/29/2024)

% Additional MATLAB code required:
% stats.m

%%
close all; clear all; 
currentFilePath = mfilename('fullpath');

%% Specify which ship dist file
file_name = 'MSC DON GIOVANNI_PROTEUS_dist_file_4000_iterations_01.mat';
N = 2; % The number of beta=infinity points used as a feataured for the maximum entropy inversion

%% Adds the correct paths
[currentFolderPath, ~, ~] = fileparts(currentFilePath);
git_library = fileparts(currentFolderPath);

% Add path to the ship .mat files
dist_files_path = fullfile(git_library,'SOO dist files');
addpath(dist_files_path) 
load(file_name)

%%
peaks           = zeros(1,4);
Ndiv            = 50;


% Separate cost from data.
cost            = dist(:,1);
min_cost        = min(cost);

if min_cost == 0
    min_cost    = typecast(uint64(1),'double');
end

T               = 2/N*min_cost;
B               = 1/T; 
Z               = sum(exp(-cost/T));

data            = dist(:,2:end); 
[niter, nparm] = size(data);

% Determine how many inversion parameters there are (nIP) and how many source level parameters there are (nSLP)
nIP = 0;
nSLP = 0;
for i = 1:nparm
    if info.type(i) == 1
        nIP = nIP + 1;
    elseif info.type(i) == 0
        nSLP = nSLP + 1;
    end
end

label={'Thickness (m)'; 'Sound Speed (m/s)'};

% Create marginal distributions by computing M-dim integrals.
% Iterate through parameters.
for i = 1:nparm
    % Notify user which parameter marginal is being found.
    fprintf(['Processing parameter ' num2str(i) '\n']); 

    % Get the lower and upper limits for the current parameter
    % Inversion parameters have actual limits, SL parameters do not
    thisparm = data(:,i); 
    if isempty(info.lim{i}) == 1                                      % These are SL parameters
        wiLB = min(thisparm);
        wiUB = max(thisparm);
    else                                                                % These are inversion parameters
        wiLB = info.lim{i}(1);
        wiUB = info.lim{i}(2);
    end

    % Divide wi space into Ndiv subspaces
    wi = linspace(wiLB,wiUB,Ndiv); 
    binwidth = wi(2)-wi(1);
    BinW(i) = binwidth;

    % Preallocate Pwi for the probability computation. 
    Pwi = zeros(1,length(wi));

    % Compute the marginal distribution.
    for j = 1:Ndiv

        % Get all occurances of wi (the delta function). 
        ndx = find( ( thisparm >= wi(j)-binwidth/2 ) & ( thisparm < wi(j)+binwidth/2 ) ); 

        % Compute the cumulative probability (over all occurances).
        if isempty(ndx) == 1
            Pwi(j) = 0;
        else
        Pwi(j) = sum(exp(-cost(ndx)/T))/(length(ndx)*Z);
        end
    end

    % Determine the normalization factor for each parameter.
    N = trapz(wi,Pwi);
    Pwi = Pwi/N;

    %Plot this parameter. 
    if i == 1

        % Setup plotting parameters.
        figure('name',['example',': Marginal Distributions, Ndiv = ',num2str(Ndiv)],'NumberTitle','off');
        nc = 4;
        nr = ceil(nparm/nc);

        % Preallocate.
        WI = zeros(nparm,length(wi));
        PWI = zeros(nparm,length(wi));
    end

        subplot(2,2,i);
        if info.type(1) == 1
            plot(wi,Pwi*binwidth,'b-','LineWidth',1.5);                 % Mult. by binwidth makes the inversion parameter plots look like they have the same area under the curve, despite the x axis
        elseif info.type(i) == 0
            plot(wi,Pwi,'b-','LineWidth',1.5);                          % Plot SL parameters differently
        end
        if i<3
            title(label{i});
        end
        
        xlim([wi(1) wi(Ndiv)])
        [~,max_id]=max(Pwi);
        peaks(i)=wi(max_id);

    % Store data in output parameters.
    WI(i,1:length(wi)) = wi;
    PWI(i,1:length(wi)) = Pwi; 
end
%
% Extract maximum y limits for plot scaling purposes.
ylimmax1 = max(max(PWI(1:nIP,:),[],2).*BinW(1:nIP).');
ylimmax0 = max(max(PWI(nIP+1:nIP+nSLP,:),[],2));
%
% Compute statistics for each parameter.
Statistics = stats(WI,PWI,nparm);
%
% Rescale the plots to all have the same vertical scale.
% Turn on stat plotting lines if desired.
statflag = 1;
%
% Create axes handles for subplots
axesHandles = zeros(nparm, 1);
for i = 1:nparm
    axesHandles(i) = subplot(2,2, i);
    if info.type(i)== 1
        ylimmax = 1.05*ylimmax1;
    elseif info(i).type == 0
        ylimmax = 1.05*ylimmax0;
    end
    ylim([0 ylimmax]) 
    if statflag == 1
        h1=line([Statistics(i,1) Statistics(i,1)],[0 ylimmax],'Color','r','LineWidth',1.2);
        h2=line([Statistics(i,1)+Statistics(i,2) Statistics(i,1)+Statistics(i,2)],[0 ylimmax],'Color',[0.6350 0.0780 0.1840],'LineStyle',':','LineWidth',1);
        h3=line([Statistics(i,1)-Statistics(i,2) Statistics(i,1)-Statistics(i,2)],[0 ylimmax],'Color',[0.6350 0.0780 0.1840],'LineStyle',':','LineWidth',1);
        h4=line([info.optimal_parm(i) info.optimal_parm(i)],[0 ylimmax],'Color',[0.4940 0.1840 0.5560],'LineStyle','--','LineWidth',1.2);
        h5=line([peaks(i) peaks(i)], [0 ylimmax],'Color',[0.4660 0.6740 0.1880],'LineStyle','-.','LineWidth',1.2);
    end
    %set(gca,'YTick',[],'YTickLabel',[])
end


% Create a common legend to the right of the subplots
lgd = legend(axesHandles(end), [h1, h2, h4, h5], {'Mean', 'Standard Deviation', 'Optimal', 'Peak'}, 'Location', 'northeastoutside', 'Orientation', 'vertical');


drawnow
%% Re-saved the dist files with added information    
info.peak_parm      = peaks;
info.mean_and_sd    = Statistics;
info.T_val          = T;
dist_file_name      = fullfile(git_library,'SOO dist files', file_name);
plot_name           = [info.ship{1},'_',info.ship{2},'_dist_plot_T=',num2str(T),'.png'];
dist_plots_name     = fullfile(git_library,'SOO dist files',plot_name);

save(dist_file_name,'info','dist')
saveas(gcf,dist_plots_name)
    
  

 