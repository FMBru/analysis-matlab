% MATLAB/OCTAVE implementation of processing code for PICOSEC
% version date: 16 Jul 2021
addpath(genpath([fileparts(fileparts(pwd)), filesep, 'Analysis' ]));
addpath '/Users/Florian/Documents/MATLAB/Picosec/Analysis'

requireSmallAreaTrigger = false;
time_diff_sigmoid = [];
baseline = 0.3;

% clear everything
close all
% run file coordinates and number of files
minuit = 0;   % set minuit for optimizer if 1 or matlab 0
run.id = '174';
run.year = '2021';
run.name = ['BEAM ' run.year ' RUN ' run.id];
%run.path=['M:\CommonProject\PICOSEC\Testbeams\July2021\Oscilloscopes\GDD\Run' run.id '\']; % use mounted Cernbox folder
run.path=['/Users/Florian/Desktop/OpticalPicosec/Data/Testbeam/Run108/Data/']; % use local folder with scope data
%%run.nfiles = find_fileNo(run.path);
%run.nfiles = 21; % override number of file sets to process (see how many trc files are in the folder)
run.lecroy_name = '-Trace-'; %['Run' run.id];

% options for processing micromegas channel (lot magic numbers from
% processing functionhave to be added here)
opts_MM.t_dead = 1;      % blanking time from start of the recording samples
opts_MM.t_prerms = 20  % blanking time before global maximum samples
%opts_MM.Ts = 1/20e9;     % sampling speed
opts_MM.Rin=50;          % input impedance
opts_MM.invert = 1;      % is inverted (1-inverted, 0-non-inverted)
opts_MM.type=0;          % detector type (0-DUT, 1-REF) for future use
%opts_MM.ch_name = ['C2' run.lecroy_name]; % LeCroy file name format
opts_MM.en_plot = 1;     % enable debugging plots

% options for processing MCP channel (lot magic numbers from processing
% functionhave to be added here)
opts_MCP.t_dead = 1;
opts_MCP.t_prerms = 20;
%opts_MCP.Ts = 1/20e9;
opts_MCP.Rin=50;
opts_MCP.invert = 1;
opts_MCP.type=1;
%opts_MCP.ch_name = ['C1' run.lecroy_name];
opts_MCP.en_plot = 0;

numberEvents = size(matchedEvents,1);
for pos = 1:numberEvents
    if mod(pos,1000)==0
        pos
    end
    
    %% DO PROCCEISNG
    tic;         % start time measurement
    k=pos;         % valid data counter
    
    event = matchedEvents(pos);
    
    if (requireSmallAreaTrigger && event.matchedScint == 1) || requireSmallAreaTrigger==false
        %event is matched with small area scintillator
        
        % read files using third party function
        ch_mm = event.signal.waveform;
        ch_mcp = event.mcp.waveform;
        
        Ts = ch_mcp(2,1) - ch_mcp(1,1);
        opts_MM.Ts = Ts;
        opts_MCP.Ts = Ts;
        
        % go trough all of the events in the file
        
        % generate virtual time vector (important to take care for trigger
        % offset in LeCroy scope)
        t_vec_mm=ch_mm(:,1);
        t_vec_mcp=ch_mcp(:,1);
        
        % subtract the earliest time
        etime=min([t_vec_mm; t_vec_mcp]);
        
        % process MM Picosec first to see if signal is valid
        
        if(minuit==1)
            MM_temp = process_signal_minuit(ch_mm(:,1)-etime,ch_mm(:,2),opts_MM);
        else
            MM_temp = process_signal_sampic(ch_mm(:,1)-etime,ch_mm(:,2),opts_MM,baseline);
            %pause(1);
            %close all
        end
        
        % if signal is valid process MCP
        if(MM_temp.fail==0)
            if(minuit==1)
                MCP_temp = process_signal_minuit(t_vec_mcp-etime,ch_mcp.y(:,m),opts_MCP);
            else
                MCP_temp = process_signal_sampic(ch_mcp(:,1)-etime,ch_mcp(:,2),opts_MCP,baseline);
            end
            
            % if MCP is valid store data to structure array
            if(MCP_temp.fail==0)
                % process tracker ID channel
                % MM_temp.event_id = process_tr_bitstream(t_vec_mcp, ch_tr.y(:,m), opts_TR);
                %MCP_temp.event_id =MM_temp.event_id;
                % store valid data into array of structures
                MM_data(k)= MM_temp;
                MCP_data(k)= MCP_temp;
                time_diff(k) = MM_data(k).cfd.time-MCP_data(k).cfd.time;
                time_diff_sigmoid(k) = MM_data(k).sigmoid.timepoint-MCP_data(k).sigmoid.timepoint;
                MCP_maxy(k) = MCP_data(k).sig.max.y;
                MM_maxy(k) = MM_data(k).sig.max.y;
                
                k=k+1;
            end
        end
        %toc
        
    end
    
end

time_diff_sigmoid_selected = [];

for i = 1:length(time_diff_sigmoid)
    %if time_diff_sigmoid(i)>-1000 && time_diff_sigmoid(i)<0
    if time_diff_sigmoid(i)>-115 && time_diff_sigmoid(i)<-100
        time_diff_sigmoid_selected = [time_diff_sigmoid_selected;time_diff_sigmoid(i)];
    end
end

%% save data to mat file
%save(['Run' run.id '.mat'], 'run', 'MM_data', 'MCP_data', 'time_diff', 'time_diff_sigmoid', 'MCP_maxy', 'MM_maxy');

%% do some initial resolution measurement
% make simple cut with respect to median time and minimum amplitude magic
% numbers :-)

time_avg_raw = median(time_diff_sigmoid_selected);  % assume mu from the median vaule
time_min = time_avg_raw-0.3;     % predicted resolution 100ps cut 3 sigma
time_max = time_avg_raw+0.3;     % left and 3 sigma right from median


% find events within the cut and cut with respect to amplirude ugly :-)
idx_cut = time_diff_sigmoid_selected>time_min & time_diff_sigmoid_selected<time_max;
%idx_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & MCP_maxy>0.1*max(MCP_maxy) & MM_maxy>0.1*max(MM_maxy);
time_diff_cut = time_diff_sigmoid_selected(idx_cut);
%time_diff_cut = time_diff_sigmoid_selected;

% plot time difference histogram
figure
histfit(time_diff_cut,50);
pd=fitdist(time_diff_cut,'normal')
xlabel('Time difference, ns');
ylabel('events');
grid
title_str = sprintf('%s: \\mu = %4.4f ns \\sigma = %4.4f ns',run.name, pd.mu, pd.sigma);
title(title_str)

% plot dual gauss histogram
figure
rng(abs(pd.mu));
h=histogram(time_diff_cut,50);
gmoptions = statset('MaxIter',10000,'TolFun',1e-9);
gm_time = fitgmdist(time_diff_cut,2,'Options',gmoptions);
g2.p=gm_time.ComponentProportion(1);
g2.mu = gm_time.mu;
g2.sigma = [gm_time.Sigma(1,1,1) gm_time.Sigma(1,1,2)];
g2.sigma_all = sqrt(g2.p*g2.sigma(1)+(1-g2.p)*g2.sigma(2)+g2.p*(1-g2.p)*(g2.mu(1)-g2.mu(2))^2);
g2.mu_all = g2.p*g2.mu(1)+(1-g2.p)*g2.mu(2);
g2.xbins=linspace(h.BinLimits(1),h.BinLimits(2),1000);
g2.pdf1 = normpdf(g2.xbins,g2.mu(1),sqrt(g2.sigma(1)));
g2.pdf1 = g2.pdf1./sum(g2.pdf1)*g2.p; % normalize first gauss
g2.pdf2 = normpdf(g2.xbins,g2.mu(2),sqrt(g2.sigma(2)));
g2.pdf2 = g2.pdf2./sum(g2.pdf2)*(1-g2.p); % normalize second gauss
g2.k_mult = sum(h.Values);
g2.binscaling = length(g2.xbins)/h.NumBins;

hold on
plot(g2.xbins, g2.pdf1*g2.k_mult*g2.binscaling);
plot(g2.xbins, g2.pdf2*g2.k_mult*g2.binscaling);
plot(g2.xbins, g2.pdf2*g2.k_mult*g2.binscaling+g2.pdf1*g2.k_mult*g2.binscaling,'LineWidth',2);
xlabel('Time difference, ns');
ylabel('events');
legend('RAW hist', 'Gauss 1', 'Gauss 2', 'Gauss combined');
grid
title_str = sprintf('2Gauss %s: \\mu = %4.4f ns \\sigma = %4.4f ns',run.name, g2.mu_all, g2.sigma_all);
title(title_str)

% plot peak amplitude histogram cut included
figure
histogram(MCP_maxy(idx_cut))
hold on
histogram(MM_maxy(idx_cut))
xlabel('Signal amplitude, V')
ylabel('events');
grid on
legend('MCP','MM');
title_str = sprintf('%s - electron peak distribution',run.name);
title(title_str)
%{
% plot electron charge
for i=1:length(MM_data)
    MM_e_charge(i)=MM_data(i).sig.charge.e_peak;
end
figure;
plot(MM_maxy, MM_e_charge,'.');
grid on
xlabel('Signal amplitude, V')
ylabel('Electron charge, pC');
title_str = sprintf('%s - Electron charge plot',run.name);
title(title_str)

%% plot endpoint for integration
for i=1:length(MM_data)
    MM_e_peak_end(i) = 1e9*Ts*(MM_data(i).sig.e_peak_end.idx-MM_data(i).sig.max.idx);
end
figure;
plot(MM_maxy, MM_e_peak_end,'.');
grid on
xlabel('Signal amplitude, V')
ylabel('Time, ns');
title_str = sprintf('%s - Electron charge end time',run.name);
title(title_str)


figure
histogram(MM_e_peak_end);
grid on
ylabel('events')
xlabel('Time, ns');
title_str = sprintf('%s - Electron charge end time',run.name);
title(title_str)
%}
toc