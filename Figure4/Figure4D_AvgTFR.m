%% Figure4_AvgTFR
% This mfile creates the average Time Frequency Representation with its corresponsing Stats.
% It loads the GranAvg Power Normalized with premov epoch (n=2 monkeys) &  the stats with are calculated based on the raw power
% Plots the mean TFR with positive & negative significant statistical clusters  p<=0.05
% low and high frequencies.
% valiable1 to modify thisF=1(<50 Hz), thisF=2 (>50 Hz)
% variable2 to modify thisB=1(threat),2(lipsmacking),3 (chew)
% local file in Mac: FigTwoMonkeys_GrandAvg_w_Stats_noAvgTF
% It also creates Supplementary 4E-H

% updated April 30th 2025

clear all
close all

addpath('/Users/yuriria/Documents/MATLAB/fieldtrip/fieldtrip-20170101')

% ------ Change this for the Data & Stats Dir ----- % 
%TargetDir='/Users/yuriria/Dropbox/FreiwaldLabMine/code/A_GitHub/Figure_4/Data_Figure4';
%StatsDir='/Users/yuriria/Dropbox/FreiwaldLabMine/code/A_GitHub/Figure_4/Data_Figure4';

TargetDir='/Volumes/Seagate/TwoMonkeys_TB/TFR/GrandAvg_nStats';
StatsDir ='/Volumes/Seagate/TwoMonkeys_TB/TFR/TFR_RawPower_AvgSess/stats_noAvgTF'


areas= {'S1', 'M1m','PMV','M3'};
beh = {'Threat', 'LipSmacking', 'Chew'}
fst={'hann','dpss'}

%% Variables to Modify depending on the desired output
thisF=2; % use 1 Freq:5-40Hz, use 2 Freq: 40-100Hz
thisB=1; % use 1 (T),  2(LS), 3 (C)
PrintOn =0;

for thisA=1:4

    fname = ['GrandAvg',areas{thisA},'_Days*_', beh{thisB},'_*_',fst{thisF},'.mat']
    disp(fname)
    sname = ['stat_2M_500ms_',beh{thisB},'_',areas{thisA},'_',fst{thisF},'.mat']
    disp(sname)

    cd(TargetDir)
    f=dir([fname]);
    load(f(1).name);
    grandavg;


    % Stats were compared  cfg.latency = [-0.8 0]; to  cfg.latency = [0 0.8];
    load([StatsDir,'/',sname])
    stat;


    %--------------- plot this Area, this Behavior ------------

    cfg=[];
    cfg.colorbar     = 'yes';
    if strcmp(fst,'hann')
        cfg.ylim         =  []%[5 40]; % hanning
    end
    cfg.xlim         = [];% [-0.6 0.6]
    cfg.zlim         = 'maxabs';

    figure(thisA);
    ff=ft_singleplotTFR(cfg, grandavg);
    xlabel('Time (s)','Fontsize', 18); ylabel('Frequency (Hz)','FontSize',18);

    idxT0 =find(grandavg.time==stat.ActTime(1)); % indx start time
    idxT1 =find(grandavg.time==stat.ActTime(end)) % indx end time


    if strcmp(areas(thisA),'M3') || strcmp(areas(thisA),'PMV') && thisF==1
        idxT1 =31; % indx end time
    end

    if strcmp(areas(thisA),'M3') || strcmp(areas(thisA),'PMV') && thisF==2
        idxT0=21
        idxT1 =31; % indx end time
    end


    yf = flipud(stat.freq)'; % this are the freq


    for thisC=1:size(stat.posclusters,2) % assess each of the clusters
        disp(['p=', num2str(stat.posclusters(thisC).prob)])

        if isempty(stat.posclusterslabelmat)==0  && isempty(stat.posclusters(thisC))~=1 && (stat.posclusters(thisC).prob<=0.05)==1
            disp(['-----',num2str(thisC), '#ThisPosCluster, pvalue =', num2str(stat.posclusters(thisC).prob)])
            poscluster = squeeze(stat.posclusterslabelmat);

            if size(stat.posclusters,2)>1
                if thisC==1
                    poscluster(find(poscluster~=thisC))=0;
                elseif thisC>1
                    poscluster(find(poscluster==thisC))=10; % make the current cluster 10
                    poscluster(find(poscluster<10))=0; % set the rest to zero
                    poscluster(find(poscluster==10))=1; % make current cluster 1
                end
            elseif size(stat.posclusters,2)==1
                poscluster(find(poscluster>0))=1;
            end

            % Matrix for Positive Clusters
            mP=repmat(yf,1,size(poscluster,2)); % remap
            mP = mP.*poscluster;
            gcf; hold on;
            s = plot(grandavg.time(idxT0:idxT1), mP,'*r', 'Markersize', 4.5);
            pause(1)
            clear poscluster
        end

    end



    for thisC=1:size(stat.negclusters,2)
        disp(['p=', num2str(stat.negclusters(thisC).prob)])
        if isempty(stat.negclusterslabelmat)==0 && isempty(stat.negclusters)~=1 && (stat.negclusters(thisC).prob<=0.05)==1
            disp(['-----', num2str(thisC), '#NegCluster, pvalue =', num2str(stat.negclusters(thisC).prob)])
            negcluster = squeeze(stat.negclusterslabelmat);
            if size(stat.negclusters,2)>1
                if thisC==1
                    negcluster(find(negcluster~=thisC))=0;
                elseif thisC>1
                    negcluster(find(negcluster==thisC))=10;% make the current cluster 10
                    negcluster(find(negcluster<10))=0; % set the rest to zero
                    negcluster(find(negcluster==10))=1; % make current cluster 1
                end
            elseif size(stat.negclusters,2)==1
                negcluster(find(negcluster>0))=1;
            end
            % Matrix for Negative Clusters
            mN = repmat(yf,1,size(negcluster,2));
            mN =mN.*negcluster;
            gcf; hold on;
            s = plot(grandavg.time(idxT0:idxT1), mN,'*b', 'Markersize', 4.5);
            clear negcluster
        end
    end

    hold off
    xlim([-0.501 0.501])
    if thisF==1
        ylim([5 40])

    elseif thisF==2
        ylim([40 121])
    end

    title([areas{thisA}, '   ' ,beh{thisB}],'FontSize',21)
    set(gca,'TickDir','out','Fontsize', 18);

    ftoP = ['GrandAvg_NoAvgTF_',areas{thisA},'_Days_', beh{thisB},'_',fst{thisF}];

    if PrintOn==1

        set(gcf, 'Position', get(0, 'Screensize'));
        set(gcf,'renderer','painters') % this is to vectorize & open in Illustrator
        print(gcf, '-depsc', ftoP);
        savefig(ftoP)
    end
    pause(3)

    close

    clear grandavg stat ix
end




