%% TwoMonkeys_stats_mod_TvsLS
%% To be able to run this file, access to the data with same trial numbers LS & T is needed (see line 50,51) 
% this mfile is a modification of (Elie's original Mfile) to assess Coh, PPC, Granger and
% store the Power during the Movement period (500ms) for T & LS. 
% It also corrects the voltage in Thor, so to have everything in
% microVolts. 
% Modified to store both Granger directions (sept/2023)
% Data & Figures stored in:'/Users/yuriria/Documents/MATLAB/Thor_Analysis/LFP_analysis/Coh/Barney_FromEli/ElieSCode/figs/TwoMonkeys/2M_Movement'

% how to know if these were the same session ?
% and how to equalise trials? start manually

% 210704-110137 has 169 trials (lipsmack) 1
% 210704-111302 has 246 (lipsmack) 2
% 210704-113250 has 114 (threat)
% 210704-115809 has 61 (threat)

% 210706-111817 has 252 (lipsmack)
% 210706-121537 has 215 (lipsmack)
% 210706-112757 has 93 (threat)
% 210706-123456 has 101 (threat)

% 210713-112225 has 147 (lipsmack)
% 210713-114719 has 69 (threat)
% 210713-125345 has 77 (threat)

% maybe best strategy is to match manually then rename such that they can
% easily be read by the stats script

% match lower number of trials within a day with lower number (and same for higher)
% rename them with number at the beginning of the name
% so 1_ in lipsmack matches with 1_ in threat, etc
clear all
close all

% this Mfile was modified  
% to plot for the Coh, PPC, & Granger among the following
% pairs of areas: M1, M3, PMV, and S1. 
% 2023,
% Figures Stored: '/Users/yuriria/Documents/MATLAB/Thor_Analysis/LFP_analysis/Coh/Barney_FromEli/ElieSCode/figs/thor')

addpath('/Users/yuriria/Documents/MATLAB/fieldtrip/fieldtrip-20170101'); %addpath /home/common/matlab/fieldtrip/
addpath('/Users/yuriria/Documents/MATLAB/fieldtrip/fieldtrip-20170101/external/spm12')
addpath('/Users/yuriria/Documents/MATLAB/Thor_Analysis/func')
ft_defaults;

%% load files, equalise trials, and select chans of interest
clear
% define folders
lipsmack_folder = '/Volumes/Seagate/Barney_CollabwElie/TwoMonkeys/lipsmack/'; %lipsmack_folder = '/project/3035003.01/Barney_ToShareWEli/PreProc_sessions/lipsmack/';
threat_folder  = '/Volumes/Seagate/Barney_CollabwElie/TwoMonkeys/threats/'; % threat_folder = '/project/3035003.01/Barney_ToShareWEli/PreProc_sessions/threat/';

% days used LS & threats
% Thor: 171005, 171010, 171027 ( 3 days, 4 sessions)  
% Barney: 210704, 210706, 210713, 210708, 210709, 210712, 210714, 210716 ( 8 days/ 10 sessions)  
% StoreDir='/Users/yuriria/Documents/MATLAB/Thor_Analysis/LFP_analysis/Coh/Barney_FromEli/ElieSCode/figs/TwoMonkeys/2M_15Sess'

% get list of files in these folder
% careful: files should match between the folders
lipsmack_files = dir(fullfile(lipsmack_folder, '*mat'))
threat_files = dir(fullfile(threat_folder, '*mat'))
assert(numel(lipsmack_files)==numel(threat_files)) % make sure equal number of files

% loop on files
fi=1;
s_trls=[];
s_pow1=[];
s_pow2=[];
s_fname=[];


for j=1:numel(lipsmack_files)

    % load Chew file
    load([lipsmack_folder lipsmack_files(j).name])
    disp([lipsmack_files(j).name])
    % rename file to avoid overwriting (they are all called data_mmov)
    data1=data_mmov; clear data_mmov
    labeld1='ls';

    % load threat file
    load([threat_folder threat_files(j).name])
    disp([threat_files(j).name])
    % rename
    data2=data_mmov; clear data_mmov
    labeld2='t';
 
    % if Thor * voltage (although for the coupling doesnt matter) 
    if isempty(strfind(lipsmack_files(j).name,'T'))==0 % thor files
        [data1]=Get_ThorVoltageInMicroVolts(data1); % *1e6
        [data2]=Get_ThorVoltageInMicroVolts(data2);
        disp(['-----', lipsmack_files(j).name, 'Voltage to microVolts---'])
    elseif isempty(strfind(lipsmack_files(j).name,'T'))==1 % barney files
    end

    % equalise number of trials
    mintrls = min(numel(data1.trial), numel(data2.trial));
    trls1 = randperm(numel(data1.trial), mintrls);
    trls2 = randperm(numel(data2.trial), mintrls);
    s_trls = [s_trls; mintrls];

    s_fname = [s_fname; {lipsmack_files(j).name} {threat_files(j).name}];

    cfg = [];
    cfg.trials = trls1;
    data1 = ft_selectdata(cfg, data1);
    cfg.trials = trls2;
    data2 = ft_selectdata(cfg, data2);

    cf=lipsmack_files(j).name(11:end);
    tf=threat_files(j).name(11:end);

    thisDayC=cf(1:end-17)
    thisDayT=tf(1:end-17)

    % select ROIs
    % as an example I will try M1m and PMVm
    cfg = [];
    if isempty(strfind(lipsmack_files(j).name,'T'))==0
         cfg.channel = 'PMV*';
    elseif isempty(strfind(lipsmack_files(j).name,'T'))==1 % barney files
         cfg.channel = 'PMVm*';
    end

    area1=cfg.channel(1:end-1);
    cfg.avgoverchan = 'yes';
    data1_ROI1 = ft_selectdata(cfg, data1); % LS
    data2_ROI1 = ft_selectdata(cfg, data2); % threat
    data1_ROI1.label = {'ROI1'};
    data2_ROI1.label = {'ROI1'};

    % select the 2nd ROI from both mov data
    cfg = [];
    cfg.channel='S1*'

    area2=cfg.channel(1:end-1);
    cfg.avgoverchan = 'yes';
    data1_ROI2 = ft_selectdata(cfg, data1); % LS
    data2_ROI2 = ft_selectdata(cfg, data2); % threat
    data1_ROI2.label = {'ROI2'};
    data2_ROI2.label = {'ROI2'};

    % combine the 2 ROIs in 1 structure for each mov type
    data1ROIs = ft_appenddata([], data1_ROI1, data1_ROI2); % LS both PMV & M1
    data2ROIs = ft_appenddata([], data2_ROI1, data2_ROI2); % T both PMV & M1

    cfg = [];
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.output = 'fourier'; % 'pow'-> keep all trials
    %cfg.keeptrials= 'yes' --> keep power per trial
    cfg.pad = 2;
    freq1_4conn = ft_freqanalysis(cfg, data1ROIs) % Chew (both areas)
    freq2_4conn = ft_freqanalysis(cfg, data2ROIs) % Threat (both areas)
    
    % granger
    % important to do granger on entire available freq spectrum,
    % can select freqs in next step
    cfg = [];
    cfg.method = 'granger';
    gran1 = ft_connectivityanalysis(cfg, freq1_4conn) % C both areas
    gran2 = ft_connectivityanalysis(cfg, freq2_4conn) % T both areas

    % reduce freqs now to make next steps faster
    cfg = [];
    cfg.frequency = [4 50];

    gran1 = ft_selectdata(cfg, gran1)
    gran2 = ft_selectdata(cfg, gran2)

    freq1_4conn = ft_selectdata(cfg, freq1_4conn)
    freq2_4conn = ft_selectdata(cfg, freq2_4conn)

    % power
    cfg = [];
    cfg.keeptrials = 'yes';

    pow1 = ft_freqdescriptives(cfg, freq1_4conn); % abs
    pow1.powspctrm = log10(pow1.powspctrm); % log
    pow1 = ft_freqdescriptives([], pow1); % mean

    pow2 = ft_freqdescriptives(cfg, freq2_4conn); % abs
    pow2.powspctrm = log10(pow2.powspctrm); % log
    pow2 = ft_freqdescriptives([], pow2); % mean
   

    figure(1); clf; hold on;
    plot(pow1.freq,pow1.powspctrm(1,:),'-g'); plot(pow1.freq, pow1.powspctrm(2,:),'--g')
    plot(pow2.freq,pow2.powspctrm(1,:),'-r'); plot(pow2.freq, pow2.powspctrm(2,:),'--r')
    legend(area1,area2, area1, area2); xlabel('Freq (Hz)'); ylabel('Power (log10)')

    s_pow1{fi} = pow1; % Modify above to get power stored per trial 
    s_pow2{fi} = pow2;

    % coherence
    cfg = [];
    cfg.method = 'coh';
    coh1 = ft_connectivityanalysis(cfg, freq1_4conn)
    coh2 = ft_connectivityanalysis(cfg, freq2_4conn)

    % ppc
    cfg = [];
    cfg.method = 'ppc';
    ppc1 = ft_connectivityanalysis(cfg, freq1_4conn);
    ppc2 = ft_connectivityanalysis(cfg, freq2_4conn);

    % store the chew connectivity measures in the existing pow struct
    conn1{fi} = pow1; % so we already have power there as the first 'channel'
    % name the 'channels'
    conn1{fi}.label = [{'pow'}; {'coh'}; {'ppc'}; {'gran 1->2'}; {'gran 2->1'}];
    % chan 2 is coh
    conn1{fi}.powspctrm(2,:) = squeeze(coh1.cohspctrm(1,2,:));
    % chan 3 is ppc
    conn1{fi}.powspctrm(3,:) = squeeze(ppc1.ppcspctrm(1,2,:));
    % chan 4 is gran 1-2
    conn1{fi}.powspctrm(4,:) = squeeze(gran1.grangerspctrm(1,2,:));
    % chan 5 is gran 2-1
    conn1{fi}.powspctrm(5,:) = squeeze(gran1.grangerspctrm(2,1,:));

    % same for the threat data
    conn2{fi} = pow2;
    % name the 'channels'
    conn2{fi}.label = [{'pow'}; {'coh'}; {'ppc'}; {'gran 1->2'}; {'gran 2->1'}];
    % chan 2 is coh
    conn2{fi}.powspctrm(2,:) = squeeze(coh2.cohspctrm(1,2,:));
    % chan 3 is ppc
    conn2{fi}.powspctrm(3,:) = squeeze(ppc2.ppcspctrm(1,2,:));
    % chan 4 is gran 1-2
    conn2{fi}.powspctrm(4,:) = squeeze(gran2.grangerspctrm(1,2,:));
    % chan 5 is gran 2-1
    conn2{fi}.powspctrm(5,:) = squeeze(gran2.grangerspctrm(2,1,:));

    fi=fi+1;

end % still need to index the outcomes conn1{fi} and conn2

% plot Coherence
s_coh1=[]; s_coh2=[];
s_ppc1=[]; s_ppc2=[];
s_g1_12=[]; s_g2_12=[];
s_g1_21=[]; s_g2_21=[];
s_f1=[]; s_f2=[];

for i=1:length(conn1)
    s_coh1(i,:)=conn1{i}.powspctrm(2,:);
    s_coh2(i,:)=conn2{i}.powspctrm(2,:);
    s_ppc1(i,:)=conn1{i}.powspctrm(3,:);
    s_ppc2(i,:)=conn2{i}.powspctrm(3,:);
    s_g1_12(i,:)=conn1{i}.powspctrm(4,:);
    s_g2_12(i,:)=conn2{i}.powspctrm(4,:);
    s_g1_21(i,:)=conn1{i}.powspctrm(5,:);
    s_g2_21(i,:)=conn2{i}.powspctrm(5,:);
    s_f1(i,:)=conn1{i}.freq;
    s_f2(i,:)=conn2{i}.freq;
end

s_f1=unique(s_f1)'; s_f2=unique(s_f2)';

disp(['Ntotal Trials=', num2str(sum(s_trls))])

figure(2); clf; hold on;
plot(s_f1,mean(s_coh1,1),'-g','LineWidth',2);
plot(s_f2,mean(s_coh2,1),'-r','LineWidth',2);
legend(['LS ', area1,'-', area2], ['T ', area1,'-', area2]);
xlabel('Frequency (Hz)','Fontsize', 18); ylabel('Coherence', 'Fontsize', 18)
title([ area1,'-', area2, '  Movement  (n=',num2str(length(lipsmack_files)), ' sess 2M, ntr=', num2str(sum(s_trls)), ')'],'FontSize',17); axis tight

figure(3); clf; hold on;
plot(s_f1,mean(s_ppc1,1),'-g','LineWidth',2);
plot(s_f2,mean(s_ppc2,1),'-r','LineWidth',2);
xlabel('Frequency (Hz)', 'Fontsize', 18); ylabel('PPC', 'Fontsize', 18)
title([area1,'-', area2, ' Movement (n=', num2str(length(lipsmack_files)),' sess 2M, ntr=', num2str(sum(s_trls)), ')'], 'FontSize',17); axis tight

figure(4); clf; hold on;
plot(s_f1,mean(s_g1_12,1),'-g','LineWidth',2);
plot(s_f2,mean(s_g2_12,1),'-r','LineWidth',2);
plot(s_f1,mean(s_g1_21,1),'--g','LineWidth',2.3);
plot(s_f2,mean(s_g2_21,1),'--r','LineWidth',2.3);
%legend( 'LS PMV-> M1', 'T PMV-> M1','LS M1--> PMV', 'T M1--> PMV');
lgh4=legend(['LS ', area1, '-->', area2], ['T ', area1, '-->', area2], ['LS ', area2 '-->' area1], ['T ', area2 '-->' area1])
xlabel('Frequency (Hz)', 'Fontsize', 18); ylabel('Granger', 'Fontsize', 18)
title([area1,'-', area2,' Movement  (n=', num2str(length(lipsmack_files)), ' sess 2M, ntr=', num2str(sum(s_trls)), ')'],'FontSize',17); axis tight


%% stat

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';

cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;

%cfg.clusterstatistic = 'wcm'; % 'weight Cluster mass' crashes this depends
%on ft_statistics_montecarlo
%cfg.clusterstatistic = 'maxsum'; % runs although slightly diff results
cfg.clusterstatistic = 'maxsize'

cfg.frequency = [4 50]; % original [4 34]

cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.alpha = 0.05;      % alpha level of the permutation test
cfg.numrandomization = 1000;

if cfg.numrandomization==1000
    disp('-------------  Number of Permutations = 1000 ---------------- OK')
else
    disp('TESTING SESSION nPermutations<1000')
    pause
end

% for t test
ll = numel(conn1);
design = [];
design(1,1:ll) = 1:ll;
design(1,(ll+1:2*ll)) = 1:ll;
design(2,1:ll) = 1;
design(2,(ll+1:2*ll))= 2;

cfg.design = design;             % design matrix
cfg.ivar  = 2;
cfg.uvar = 1;
cfg.neighbours = [];

for thisCh=[2:5]

    % power (cfg.channel = 1)
    % coherence (cfg.channel =2)
    % ppc (cfg.channel = 3)
    % granger PMV - > M1  % ( cfg.channel = 4 )
    % granger M1 -> PMV ( cfg.channel = 5 )
    cfg.channel = thisCh; % change the channel according to which measure you want to test

    stat = ft_freqstatistics (cfg, conn1{:}, conn2{:});

    %% check stats
    lw = 15
    if cfg.channel==2 % Coherence
        var1=mean(s_coh1,1); var2=mean(s_coh2,1);
        n=2;
        storename = ['Coh_FX_', area1,'-',area2];
        nl=2;
    elseif cfg.channel==3 % PPC
        var1=mean(s_ppc1,1); var2=mean(s_ppc2,1);
        n=3;
        storename = ['PPC_FX_', area1,'-',area2];
        nl=2;
    elseif cfg.channel==4 % Granger area1--> area 2 
        var1=mean(s_g1_12,1); var2=mean(s_g2_12,1);
        n=4;
        nl=4;
        storename = ['Granger_FX_', area1,'-',area2];
    elseif cfg.channel==5 % Granger area2 --> area 1 
        var1=mean(s_g1_21,1); var2=mean(s_g2_21,1);
        n=4;
        lw = 20;
        storename = ['Granger_FX_', area2,'-',area1];
        nl=4;
    end
    

   % Assess Significant Negative Clusters
    if isfield(stat,'negclusters')==1 && isempty(stat.negclusters)~=1
        figure(n); gcf; hold on;
        for i=1:size(stat.negclusters,2)
            if stat.negclusters(i).prob<=0.05 % changed from 0.5
                [y,idx]=ismember(stat.freq(stat.negclusterslabelmat==i), stat.freq);
                plot(stat.freq(idx),var1(idx), 'Color',[0.1, 0.9, 0.2, 0.2],'LineWidth', lw);
                plot(stat.freq(idx),var2(idx), 'Color',[1, 0, 0, 0.2],'LineWidth',lw)
            end
        end
    end

    % Assess Significant Positive Clusters
    if  isfield(stat,'posclusters')==1 && isempty(stat.posclusters)~=1 
        figure(n); gcf; hold on;
        for i=1:size(stat.posclusters,2)
            if stat.posclusters(i).prob<=0.05 % changed from 0.5
                [y,idx]=ismember(stat.freq(stat.posclusterslabelmat==i), stat.freq);  
                plot(stat.freq(idx),var1(idx), 'Color',[0.1, 0.9, 0.2, 0.2],'LineWidth', lw);
                plot(stat.freq(idx),var2(idx), 'Color',[1, 0, 0, 0.2],'LineWidth',lw);
            end
        end
    end


    if cfg.channel<4
        figure(n); gcf; hold on;
        lgh=legend(['LS ' area1,'-', area2], ['T ', area1,'-', area2]);
        lgh.String=lgh.String(1:2);
        lgh.FontSize=14
    elseif cfg.channel==5
        figure(n); gcf; hold on;
        lgh4.String=lgh4.String(1:4);
        lgh.FontSize=14
    end

    hold off

    PrintOn=1;
    xlim([4 50])

    if PrintOn==1 && cfg.channel<4 
        figure(n);
        ax = gca;
        ax.FontSize = 20; 
        cd('/Users/yuriria/Documents/MATLAB/Thor_Analysis/LFP_analysis/Coh/Barney_FromEli/ElieSCode/figs/TwoMonkeys/2M_Movement')
        % set(gcf, 'Position', get(0, 'Screensize'))
        % print(gcf,'-dpng',storename);
        % savefig(storename)
        % save(storename, 'var1', 'var2', 'stat','area1', 'area2','labeld1', 'labeld2','s_trls')
        % 
        % if cfg.channel==2
        %     save(['Power_FX_', area1,'-',area2], 's_pow1', 's_pow2','area1', 'area2', ...
        %         'labeld1', 'labeld2', 's_trls', 's_fname')
        % end

    elseif PrintOn==1 && cfg.channel==4 % modify so you store both directions! 
          %save(storename, 'var1', 'var2', 'stat', 'area1', 'area2','labeld1', 'labeld2','s_trls','s_fname')

    elseif PrintOn==1 && cfg.channel==5 % modify so you store both directions! 
        figure(n);
        ax = gca;
        ax.FontSize = 20; 
        ylim([0 0.7])
        cd('/Users/yuriria/Documents/MATLAB/Thor_Analysis/LFP_analysis/Coh/Barney_FromEli/ElieSCode/figs/TwoMonkeys/2M_Movement')
        %set(gcf, 'Position', get(0, 'Screensize'))
        %print(gcf,'-dpng',storename);
        %savefig(storename)
        %save(storename, 'var1', 'var2', 'stat', 'area1', 'area2','labeld1', 'labeld2')

    end
       
end









