%% Sample Script to Read Kiel Smarting Data
%
% 1. Load the preprocessed data
% 2. Apply additional filters
% 3. Compute the ERP across trials
% 4. Apply a baseline correction
% 5. Compute the grand average across participants
% 6. Extract ERP peaks for statistics outside of matlab
% 7. Alternative: Statistics in Matlab

%% Set the basics
% Go to Current Folder
clear all
close all
ft_defaults; % Set the defaults of the FieldTrip Toolbox
behavtab = readtable('C:\Users\...');

% Where are the data?
inpath = ('C:\Users\...');
% What are the data called?
indat = dir([inpath,'*_preproc.mat']);
%Standardize names across behavtab and indat
for n = 1 : length(indat)
    indat(n).name = erase(indat(n).name,".xdf_preproc.mat");
end

%% Match indat name to behavtab name in order to include behavioral data
for v = 1:size(behavtab,1)
    for e = 1:size(indat,1)
        ismatch(e) = strcmpi(behavtab{v,1},indat(e).name);
    end
    if ~any(ismatch)
        removerow(v) = 1;
    end
end
behavtab(find(removerow),:) = [];
behavtab = sortrows(behavtab,'vpcode');

%% Loop participants
for v = 1:length(indat)
    %% 1. Load the preprocessed data
    load([inpath,indat(v).name,'.xdf_preproc.mat']);
    if exist('data_ci','var')
        data_cif = data_ci
        clear data_ci
    end
    
    % Fix time vector
    for t = 1:length(data_cif.trial)
        data_cif.time{t} = -1.5:.002:1.5;
    end

    %% 2. Apply additional Low-Pass-Filter for ERPs
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 15;
    cfg.lpfilttype = 'firws';
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 1;
    cfg.hpfilttype = 'firws';

    data_cif = ft_preprocessing(cfg,data_cif);

    %% Split pre and post
    samplediff = diff(data_cif.sampleinfo(:,1));    %Ableitung der zeitlichen Verteilung der Samplepunkte
    [pks, locs] = max(samplediff);                  %Speichern des höchsten Wertes und seiner Position
                                                    %Nachfolgend verwenden
                                                    %der höchsten Steigung
                                                    %als Cutoff
    cfg = [];
    cfg.trials = 1:locs;

    data_pre = ft_selectdata(cfg,data_cif);

    cfg = [];
    cfg.trials = locs+1:length(data_cif.trial);

    data_post = ft_selectdata(cfg,data_cif);
    %% Trial Selection (1 = kongruent, 2 = inkongruent, 3 = korrekt, 4 = inkorrekt)
    %all
    all_congtrials = find(data_cif.trialinfo == 13);    %korrekte Trials in der kongruenten Bedingung
    all_incongtrials = find(data_cif.trialinfo == 23);  %korrekte Trials in der inkongruenten Bedingung

    cfg =[];
    cfg.trials = all_congtrials;

    data_all_cong = ft_selectdata(cfg,data_cif);

    cfg.trials = all_incongtrials;

    data_all_incong = ft_selectdata(cfg,data_cif);
    
    %pretest
    pre_congtrials = find(data_pre.trialinfo == 13);    %korrekte Trials in der kongruenten Bedingung
    pre_incongtrials = find(data_pre.trialinfo == 23);  %korrekte Trials in der inkongruenten Bedingung

    cfg =[];
    cfg.trials = pre_congtrials;

    data_pre_cong = ft_selectdata(cfg,data_pre);

    cfg.trials = pre_incongtrials;

    data_pre_incong = ft_selectdata(cfg,data_pre);

    %posttest
    post_congtrials = find(data_post.trialinfo == 13);    %korrekte Trials in der kongruenten Bedingung
    post_incongtrials = find(data_post.trialinfo == 23);  %korrekte Trials in der inkongruenten Bedingung

    cfg =[];
    cfg.trials = post_congtrials;

    data_post_cong = ft_selectdata(cfg,data_post);

    cfg.trials = post_incongtrials;

    data_post_incong = ft_selectdata(cfg,data_post);

    %% 3. Average across trials: ERP
    %all
    ERP_all_cong{v} = ft_timelockanalysis([],data_all_cong);
    ERP_all_incong{v} = ft_timelockanalysis([],data_all_incong);
    %pretest
    ERP_pre_cong{v} = ft_timelockanalysis([],data_pre_cong);
    ERP_pre_incong{v} = ft_timelockanalysis([],data_pre_incong);
    %posttest
    ERP_post_cong{v} = ft_timelockanalysis([],data_post_cong);
    ERP_post_incong{v} = ft_timelockanalysis([],data_post_incong);
    %% 4. Correct for offset: Baseline-Correction
    cfg = [];
    cfg.baseline = [-.2 -.05]; % Keep a bit away from stimulus onset
    %all
    ERP_all_cong_bl{v} = ft_timelockbaseline(cfg,ERP_all_cong{v});
    ERP_all_incong_bl{v} = ft_timelockbaseline(cfg,ERP_all_incong{v});
    %pretest
    ERP_pre_cong_bl{v} = ft_timelockbaseline(cfg,ERP_pre_cong{v});
    ERP_pre_incong_bl{v} = ft_timelockbaseline(cfg,ERP_pre_incong{v});
    %posttest
    ERP_post_cong_bl{v} = ft_timelockbaseline(cfg,ERP_post_cong{v});
    ERP_post_incong_bl{v} = ft_timelockbaseline(cfg,ERP_post_incong{v});
end % Loop across participants is done
    %ft_math differenz zwischen conditions
    %für jede Versuchsperson mindEx/MAAS stringmatching

%% 5. Compute the Grand Average across participants
%all
GA_all_cong = ft_timelockgrandaverage([],ERP_all_cong_bl{:});
GA_all_incong = ft_timelockgrandaverage([],ERP_all_incong_bl{:});
%pretest
GA_pre_cong = ft_timelockgrandaverage([],ERP_pre_cong_bl{:});
GA_pre_incong = ft_timelockgrandaverage([],ERP_pre_incong_bl{:});
%posttest
GA_post_cong = ft_timelockgrandaverage([],ERP_post_cong_bl{:});
GA_post_incong = ft_timelockgrandaverage([],ERP_post_incong_bl{:});

    % 5.2. You can also look at all channels at the same time
    %play with arguments to get a nice image for publication
    cfg = [];
    %cfg.channel = ft_channelselection('F*',GA_all_cong);
    cfg.xlim = [-.5 .5]; % Set the interval to display
    cfg.layout = 'EEG1020.lay'; % Set the channel layout
    cfg.showlabels    = 'yes';
    %all
    ft_multiplotER(cfg,GA_all_cong,GA_all_incong);
    %pretest
    ft_multiplotER(cfg,GA_pre_cong,GA_pre_incong);
    %posttest
    ft_multiplotER(cfg,GA_post_cong,GA_post_incong);
    
    %% 6. Extract Peaks for R

% FT Math für Berechnung der inhibitorischen Kontrolle

%all
for v = 1:length(ERP_all_cong_bl) %loop across participants

    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';

    ERP_all_inhcon_bl{v} = ft_math(cfg,ERP_all_incong_bl{v},ERP_all_cong_bl{v});
end

%pretest
for v = 1:length(ERP_pre_cong_bl) %loop across participants

    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';

    ERP_pre_inhcon_bl{v} = ft_math(cfg,ERP_pre_incong_bl{v},ERP_pre_cong_bl{v});
end

%posttest
for v = 1:length(ERP_post_cong_bl)

    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';

    ERP_post_inhcon_bl{v} = ft_math(cfg,ERP_post_incong_bl{v},ERP_post_cong_bl{v});
end
    
%N200 Peaks
%pre
pre_starti_inhcon = nearest(ERP_pre_inhcon_bl{1}.time,.18); % Start point
pre_endi_inhcon = nearest(ERP_pre_inhcon_bl{1}.time,.325); % End point

%pre incong_peak
pre_starti_incong = nearest(ERP_pre_incong_bl{1}.time,.18); % Start point
pre_endi_incong = nearest(ERP_pre_incong_bl{1}.time,.325); % End point

%post
post_starti_inhcon = nearest(ERP_post_inhcon_bl{1}.time,.18); % Start point
post_endi_inhcon = nearest(ERP_post_inhcon_bl{1}.time,.325); % End point

%post incong_peak
post_starti_incong = nearest(ERP_post_incong_bl{1}.time,.18); % Start point
post_endi_incong = nearest(ERP_post_incong_bl{1}.time,.325); % End point

% Fz electrode for N2 Peaks
for v = 1:length(ERP_pre_incong_bl) % Loop VPn
    cfg = [];
    cfg.channel = {'Fp1'};
    
    ROI_pre_inhcon_bl{v} = ft_selectdata(cfg,ERP_pre_inhcon_bl{v});
    ROI_pre_incong_bl{v} = ft_selectdata(cfg,ERP_pre_incong_bl{v});
    ROI_post_inhcon_bl{v} = ft_selectdata(cfg,ERP_post_inhcon_bl{v});
    ROI_post_incong_bl{v} = ft_selectdata(cfg,ERP_post_incong_bl{v});
end

    % Pick individual Peaks for N2
    peaks = zeros(length(ROI_pre_inhcon_bl),1,4);
    lat = zeros(length(ROI_pre_inhcon_bl),1,4);
    for v = 1:length(ROI_pre_inhcon_bl) % Loop VPn
        for c = 1:size(ROI_pre_inhcon_bl{v}.avg,1) % Loop Channels
            [peaks(v,c,1) lat(v,c,1)] = max(ROI_pre_inhcon_bl{v}.avg(c,pre_starti_inhcon:pre_endi_inhcon)); %N2 pre inhcon
            [peaks(v,c,2) lat(v,c,2)] = max(ROI_pre_incong_bl{v}.avg(c,pre_starti_incong:pre_endi_incong)); %N2 pre incong
            [peaks(v,c,3) lat(v,c,3)] = max(ROI_post_inhcon_bl{v}.avg(c,post_starti_inhcon:post_endi_inhcon)); % N2 post inhcon
            [peaks(v,c,4) lat(v,c,4)] = max(ROI_post_incong_bl{v}.avg(c,post_starti_incong:post_endi_incong)); %N2 pre
        end % for c
    end % for v

    % 6.2. Option 2: Average over interval
    peaks = zeros(length(ROI_pre_inhcon_bl),3,2);
    for v = 1:length(ROI_pre_inhcon_bl) % Loop VPn
        for c = 1:size(ROI_pre_inhcon_bl{v}.avg,1) % Loop Channels
            peaks(v,c,1) = mean(ROI_pre_inhcon_bl{v}.avg(c,starti:endi));
            peaks(v,c,2) = mean(ROI_post_inhcon_bl{v}.avg(c,starti:endi));
        end % for c
    end % for v

    % 6.3. Option 3: Pick Peaks based on GA and average around them
    GA_ROI_all = ft_timelockgrandaverage([],ROI_all_bl{:});
    for c = 1:size(GA_ROI_all.avg,1) % Loop Channels
        [tmppeaks(c) tmplat(c)] = max(abs(GA_ROI_all.avg(c,starti:endi)));
    end

    peaks = zeros(length(ROI_pre_inhcon_bl),3,2);
    for v = 1:length(ROI_pre_inhcon_bl) % Loop VPn
        for c = 1:size(ROI_pre_inhcon_bl{v}.avg,1) % Loop Channels
            peaks(v,c,1) = mean(ROI_pre_inhcon_bl{v}.avg(c,starti+tmplat(c)-10 : starti+tmplat(c)+10));
            peaks(v,c,2) = mean(ROI_post_inhcon_bl{v}.avg(c,starti+tmplat(c)-10 : starti+tmplat(c)+10));
        end % for c
    end % for v

    % 6.4. Save Values to use in R
    vp = {indat.name};
    %1 = N2 pre inhcon; 2 = N2 pre incong; 3 = N2 post inhcon; 4 = N2 post
    %incong
    save('N2_Fp1_peaks_lats.mat','peaks','lat','vp');
%% 7. Stats in Matlab
% In Matlab, you can either work with the extracted peaks, or you can work
% with the entire dataset. The latter is called "Mass Univariate Approach"

    % 7.2. Option 2: Mass Univariate Approach
    % Correct for temporal and spatial clusters (c.f. Maris & Oostenveld,
    % 2007)
    
        % 7.2.1. Define Neighbours
        cfg = []; 
        cfg.method = 'distance'; % how should the neighbors be selected?
        cfg.neighbourdist = 50; % I have no Idea what range this has, just make sure, that you get meaningful neighbors
        cfg.elec = 'standard_1020.elc'; % Here we need the 3d-positions!
        
        neigh = ft_prepare_neighbours(cfg); % between 5 and 10 neighbors is a good value, always good to check!
        
        % 7.2.2. Compute Stats
        cfg = [];
        cfg.parameter = 'avg'; % On which parameter?
        cfg.method = 'montecarlo'; % non-parametric stats based on montecarlo simulation
        cfg.numrandomization = 2000; % Number of steps in the simulation
        cfg.correctm = 'cluster'; % cluster-correction
        cfg.neighbours = neigh; % Define neighbors
        cfg.statistic = 'depsamplesT'; % Dependent-samples t-test
        cfg.correcttail = 'alpha'; % Correct the alpha-level for 2-tailed test
        %cfg.uvar = 1; % How many "units" = VPn
        cfg.ivar = 2; % Who belongs to which group
        % Set the design to correspondent to uvar and ivar
        % First row: Two vectors from 1 to the number of VPn
        % Second row: Two vector, one of 1s and one of 2s
        cfg.design = [1:length(ERP_all_incong_bl),1:length(ERP_all_cong_bl);...
                        ones(1,length(ERP_all_incong_bl)),ones(1,length(ERP_all_cong_bl)).*2];
        % Restrict the data of interest
        cfg.channel = 'all'; % all channels (you can also set the ROI here)
        cfg.avgoverchan = 'no'; % Should we average actross a ROI?
        cfg.latency = [0 1]; % The post-stimulus interval
        cfg.avgovertime = 'no'; % Should we average across time

        % And compute the stats
        stats = ft_timelockstatistics(cfg,ERP_all_incong_bl{:},ERP_all_cong_bl{:});

        % 7.2.3. Plot Unrestricted Stats
        cfg = [];
        cfg.layout = 'EEG1020.lay';
        cfg.parameter = 'stat';
        cfg.maskparameter = 'mask';

        ft_multiplotER(cfg,stats);

        % 7.2.4 Find the Cluster
        [chan lat] = find(stats.posclusterslabelmat == 1);
        chan = unique(chan)
        lat = unique(lat);

        starti = stats.time(lat(1))
        endi = stats.time(lat(end))
        
        figure; hold;
        plot(GA_pre_cong.time,squeeze(mean(GA_pre_cong.avg(chan,:))),'b');
        plot(GA_pre_incong.time,squeeze(mean(GA_pre_incong.avg(chan,:))),'r');

%% Correlation inhibitorische Kontrolle pre ~ mindEx
% FT Math für die Differenz
for v = 1:length(ERP_pre_cong_bl)

    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';

    ERP_pre_inhcon_bl{v} = ft_math(cfg,ERP_pre_incong_bl{v},ERP_pre_cong_bl{v});
end

% 7.2.2. Compute Stats
        cfg = [];
        cfg.parameter = 'avg'; % On which parameter?
        cfg.method = 'montecarlo'; % non-parametric stats based on montecarlo simulation
        cfg.numrandomization = 2000; % Number of steps in the simulation
        cfg.correctm = 'cluster'; % cluster-correction
        cfg.neighbours = neigh; % Define neighbors
        cfg.statistic = 'vt_statfun_corbrainbeh'; % Dependent-samples t-test
        cfg.correcttail = 'alpha'; % Correct the alpha-level for 2-tailed test
        %cfg.uvar = 1; % How many "units" = VPn
        cfg.ivar = 1; % Who belongs to which group
        % Set the design to correspondent to uvar and ivar
        % First row: Two vectors from 1 to the number of VPn
        % Second row: Two vector, one of 1s and one of 2s
        cfg.design = [1:length(ERP_pre_inhcon_bl);table2array(behavtab(:,5))'];

        % Restrict the data of interest
        cfg.channel = 'Fz'; % all channels (you can also set the ROI here)
        cfg.avgoverchan = 'no'; % Should we average actross a ROI?
        cfg.latency = [0 1]; % The post-stimulus interval
        cfg.avgovertime = 'no'; % Should we average across time

        % And compute the stats
        stats_cor = ft_timelockstatistics(cfg,ERP_pre_inhcon_bl{:});

        stats_cor.mask2 = stats_cor.posclusterslabelmat == 1;
        % 7.2.3. Plot Unrestricted Stats
        cfg = [];
        cfg.layout = 'EEG1020.lay';
        cfg.parameter = 'stat';
        cfg.maskparameter = 'mask2';

        ft_multiplotER(cfg,stats_cor);
        
        % 7.2.4 Find the Cluster
        [chan lat] = find(stats_cor.posclusterslabelmat == 1);
        chan = unique(chan)
        lat = unique(lat);

        starti = stats_cor.time(lat(1))
        endi = stats_cor.time(lat(end))
        
        figure; hold;
        plot(GA_pre_cong.time,squeeze(mean(GA_pre_cong.avg(chan,:))),'b');
        
        % Scatterplot Mittelwert jeder Person in dem Zeitraum des
        % signifikanten Clusters x behavtab
        t1 = nearest(ERP_pre_inhcon_bl{1}.time, starti)
        t2 = nearest(ERP_pre_inhcon_bl{1}.time, endi)

        % Hier einen for-loop
        mean(ERP_pre_inhcon_bl{1}.avg(6,t1:t2))
