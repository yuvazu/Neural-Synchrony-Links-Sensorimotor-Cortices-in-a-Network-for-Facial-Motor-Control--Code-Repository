%% Get_dPT_AUC_InCleanwM3Rec_meLFPS_Fix6msTime_Rebuttal2
% this file includes read, calculates Pk-Thr, and AUC
% for S1, M1m, PMVm,M3 as a stim area, and recordings in S1, M1m, PMVm, M3
% Peaks and Through ( black lines)
% AUC on the range of the peak and through (red lines)
% The output of this mfile is for each Stim area, the pks/trs within the
% same onset time (6ms) and variable offset. 
% - rebuttal round 2 (single time onset 6 ms (107) 
% - yvz oct 17th, 2025

clear all
close all
PrintOn=0; % set to 1 for saving files

store_ICSa = [];
store_A = [];
store_AUC = [];
store_dPT = [];
store_LatT =[];
store_LatP =[];
store_tAUC = [];
store_tPT = [];
store_t = [];

TargetDir='/Users/yuriria/Documents/MATLAB/Fx_Con_Analysis/LFP_related/Cleaned_eLFP/test_241127';
cd(TargetDir)
thisA='S1'; % this should be edited according to the area you are looking at (S1, M1m, PMVm, M3) 
load(['Barney_ICS_', thisA, '_eLFP_wM3Record.mat'])


% Colors of the Arrays in Barney !
c1=[0.7176    0.2745    1.0000];
c2= [ 0.8510    0.3255    0.0980];
c3= [ 0.6353    0.0784    0.1843];
c4 = [1 0 1];
c5 = [0.9294    0.6941    0.1255] % PMV
c6 = [0.9294    0.6941    0.1255];
c7= [0.4667    0.6745    0.1882];% M3
c = [c1; c2; c3; c4; c5; c6; c7];


% Default timing of Signal
t=[-100:200]; %time axis
% index 101 is time zero


if strcmp(thisA,'S1')==1 
    all_z = [{s_z1},{s_z4},{s_z6}];
    labelAreas;
elseif strcmp(thisA,'M1m')==1 
    all_z = [{s_z1}, {s_z4},{s_z6}];
    labelAreas;  
elseif strcmp(thisA, 'PMVm')==1 
     all_z = [{s_z1}, {s_z2}, {s_z6}];
     labelAreas; 
elseif strcmp(thisA,'M3')==1 
    all_z = [{s_z1} {s_z2}, {s_z5}];
    labelAreas;
end

figure(1); clf;

for a=1:length(all_z)  %all Z is 3 
    toffset = 107; % 6 ms from StimON (ok for M3, M1m
    tend = 154; % 53 ms away from Stim ON

    avg = mean(all_z{a},1);
    thisAl=labelAreas(a);

    if strcmp(thisAl{:},'S1')==1
        thisC=c1;
        ip=1;
    elseif strcmp(thisAl{:},'M1m')==1
        thisC=c2;
        ip=2;
    elseif strcmp(thisAl{:},'PMVm')==1
        thisC=c5;
        ip=3;
    elseif strcmp(thisAl{:},'M3')==1
        thisC=c7;
        ip=4;
    end

    if strcmp(thisA,'S1')==1 && a==2; % ok % ICMS in S1  & PMV (receiver) 
         tend = 190; % (89 )
    end

    if strcmp(thisA,'S1')==1 && a==3; % ok % ICMS in S1  & M3 (receiver) 
         tend = 154;
    end

    if strcmp(thisA,'M1m')==1 && a==2; % ok
        tend = 200; % 55 ms
    end

    if strcmp(thisA,'M1m')==1 && a==3; % ok
        tend = 151; % 50 ms
    end
    
    if strcmp(thisA,'PMVm')==1 && a==1; % 
        tend = 181; %80
    end

    if strcmp(thisA,'PMVm')==1 && a==3; % 
        tend = 161; % 60 ms
    end

    trange= toffset:tend;  
   
    subplot(2,2,ip); hold on;
    xline(t(trange(1)),'LineWidth',2.2);
    xline(t(trange(end)), 'LineWidth',2.2);
    thisT = t(trange); % time within the trange (5-150 ms) in one area 2 
    thisAvg = avg(trange); %  Signal in Time Range where peaks/Through will be found
  

    %  find the Peak to Through after the Artefact
    [pks,pklocs]=findpeaks(thisAvg); % this will take only the peaks after 2ms after StimOffset
    [troughs,trlocs] = findpeaks(-thisAvg);


    %% calculate Max & Minium Peaks but not necessary the first ones
    pktr = [pklocs(:), pks(:); trlocs(:) -troughs(:)];
    pktrs = sortrows(pktr);
    PeakMax=max(pktrs(:,2)); % assess the Max Value = max Pk
    ThrMin=min(pktrs(:,2)); % assess the Minimum Value = max Thr
    ixP=find(thisAvg==PeakMax);
    ixT=find(thisAvg==ThrMin,1);
    dPT = abs(PeakMax - ThrMin);

    %%  AUC in the timerange: 3 to 50 ms after Stimulation Onset
    if ixP < ixT; 
        smaller_number=ixP; larger_number=ixT; 
    else 
        smaller_number=ixT; larger_number=ixP; 
    end
    toffset1=smaller_number; % 3ms 
    tend1= larger_number; % 50 ms 

    trange1= trange(smaller_number):trange(larger_number); % 3-53 ms
    subplot(2,2,ip); hold on;
    xline(t(trange1(1)),'LineWidth',0.75,'color',[1 0 0]);
    xline(t(trange1(end)), 'LineWidth',0.75, 'color',[1 0 0]);
    thisT1 = t(trange1); % time within the trange (5-150 ms)
    thisAvg1 = avg(trange1); %  Signal in Time Range where peaks/Through will be found

    gt0 = thisAvg1>0; % thisAvg is the voltage trace within the

    if length(thisT1(gt0))>1 && length(thisT1(~gt0))>1
        posArea = trapz(thisT1(gt0), thisAvg1(gt0));
        negArea = trapz(thisT1(~gt0), thisAvg1(~gt0));
        auc = posArea + abs(negArea);
    elseif length(thisT1(gt0))<2 && length(thisT1(~gt0))>1
        negArea = trapz(thisT1(~gt0), thisAvg1(~gt0));
        auc = abs(negArea);
    elseif length(thisT1(gt0))>1 && length(thisT1(~gt0))<2
        posArea = trapz(thisT1(gt0), thisAvg1(gt0));
        auc = posArea;
    end
  

    % signal, timeRange, Peak & Through
    plot(t,avg,'-','Color',thisC,'linewidth', 3)
    plot(thisT(ixP),thisAvg(ixP),'*r','MarkerSize', 15, 'LineWidth',2);
    plot(thisT(ixT),thisAvg(ixT),'*b','MarkerSize', 15, 'LineWidth',2);
    plot(thisT,thisAvg,'--','Color',[0.9412    0.9412    0.9412],'LineWidth',0.75)
    yline(0); xline(0); legend(thisAl{:},'','','')
    title(['ICS ', thisA, ' eLFP ', thisAl{:},'  d=', num2str(dPT,2), '   tP=', num2str(thisT(ixP)), '   tT=', num2str(thisT(ixT))])
    xlabel('time (ms)'); 
    ylabel('z-score')
  
    
    store_ICSa = [store_ICSa; thisA];
    store_A = [store_A; thisAl];
    store_dPT = [store_dPT; dPT];
    store_AUC = [store_AUC; auc];
    store_tAUC = [store_tAUC; [t(trange1(1)) t(trange1(end))]]; 
   
    store_tPT = [store_tPT; [t(trange(1)) t(trange(end))]]; % time range where signal was found
    store_LatT = [store_LatT; thisT(ixT)];
    store_LatP = [store_LatP; thisT(ixP)];
 
    xlim([-50 150])
    ylim([-20 20])
  
end

keyboard

Storedir ='/Users/yuriria/Documents/MATLAB/Fx_Con_Analysis/LFP_related/Cleaned_eLFP/test_241127/Clean_dPT_AUC_Fixed6mstoVarOffset/SameOffset_afterStim';
cd(Storedir)

if PrintOn==1
    set(gcf, 'Position', get(0, 'Screensize'));
    print(gcf,'-depsc',[thisA,'Clean_dPT_AUC_sameToffset'])

    %related to these Original figure names
    print(gcf,'-dpng',[thisA,'Clean_dPT_AUC_sameToffset'])
    savefig([thisA,'Clean_dPT_AUC_sameToffset'])
    save([thisA,'Clean_dPT_AUC_Times_sameToffset'], 'store_*')
end