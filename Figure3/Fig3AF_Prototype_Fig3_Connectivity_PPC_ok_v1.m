function Fig3AF_Prototype_Fig3_Connectivity_PPC_ok_v1
%% Fig3AF_Prototype_Fig3_Connectivity_PPC_ok_v1
% this mfile creates Figure 3A-C for PPC Threats vs LipSmacks 
% It loads the mean PPC across all sessions for Threats &
% Lipsmacks, and plots them.  A graph is created showing the mean
% PPC for both behaviors for medial-lateral and lateral-lateral motor
% areas, note S1 PPC has its own mfile. 

% It also creates the corresponding Granger combinations for each pairs of
% areas. 
% All the stats for PPC & Granger combinations were run using the following
% parameters:
% alpha=0.05, nperm=1000. var(1): LS, var(2):T 
% All PPC & Granger calculations were done using: TwoMonkeys_stats_mod_TvsLS (based on Rassi E code)
% (n = 15 sessions , 2 monkeys), 
% local file in Mac: Prototype_Fig3_Connectivity_PPC_ok_v1
% updated April 30th, 2025


clear all
close all

% --change this path to where the data should be read from ---
%TargetDir='/Users/yuriria/Documents/MATLAB/Thor_Analysis/LFP_analysis/Coh/Barney_FromEli/ElieSCode/figs/TwoMonkeys/2M_Movement';
% local dir
TargetDir='/Users/yuriria/Dropbox/FreiwaldLabMine/code/A_GitHub/Figure_3/Data_PPC_Fig3'; % DB dir
cd(TargetDir)

strM='PPC';
lw=10;
pvalue=0.05; % for plotting purposes 



% -------- M1-PMV PPC/Coherence -----------
strAA={'M1m-PMVm',  'M1m-M3', 'M3-PMVm'};

for thisP=1:3
    % Load PPC for this pair of areas 
    load([strM,'_FX_', strAA{thisP},'.mat']);
    disp([strM,'_FX_', strAA{thisP}])
    var1; % LS
    var2; %  T

    [pks, locs] = findpeaks(var1); [ii,ia]=sort(var1(locs),'descend');
    FreqMaxPk = stat.freq(locs(ia(1:5))); PPCMaxPk = var1(locs(ia(1:5)));
    disp('LS')
    disp([FreqMaxPk' PPCMaxPk'])

    [pks, locs] = findpeaks(var2); [ii,ia]=sort(var2(locs),'descend');
    FreqMaxPk = stat.freq(locs(ia(1:5))); PPCMaxPk = var2(locs(ia(1:5)));
    disp('T')
    disp([FreqMaxPk' PPCMaxPk'])
 
    hold off

    stat.freq;
    disp([strM,'_FX_', strAA{thisP}])
    disp(['Stats PPC alpha =', num2str(stat.cfg.alpha)])
    disp(['Stats PPC nperm =', num2str(stat.cfg.numrandomization)])

    figure(1); subplot(2,3,thisP); hold on;
    plot(stat.freq,var1,'-','linewidth', 2,'Color', [0.9882    0.9255    0.0196]) % yellow
    plot(stat.freq,var2,'-','linewidth', 2,'Color', [0.8510    0.3255    0.0980]) % red

    % Assess Significant Negative Clusters T>LS
    if isfield(stat,'negclusters')==1 && isempty(stat.negclusters)~=1
        for i=1:size(stat.negclusters,2)
            alpha=stat.negclusters(i).prob;
            if alpha<=pvalue % if this cluster is significant
                alpha;
                stat.freq(stat.negclusterslabelmat==i);
                disp('T>L')
                [y,idx1]=ismember(stat.freq(stat.negclusterslabelmat==i), stat.freq);
                plot(stat.freq(stat.negclusterslabelmat==i),var1(idx1), 'Color',[0.9882    0.9255    0.0196, 0.9],'LineWidth', lw); % yellow
                plot(stat.freq(stat.negclusterslabelmat==i),var2(idx1), 'Color',[0.8510    0.3255    0.0980, 0.9],'LineWidth',lw); % red
                %pause
                clear y
                clear idx1
            end
        end
    end

    clear alpha

    xlabel('Freq (Hz)','fontsize',14); ylabel(strM, 'fontsize',14, 'fontname', 'Arial')
    title(strAA{thisP},'FontSize', 12,'fontname', 'Arial')
    ax = gca(gcf)
    ax.XAxis.FontSize = 13;
    ax.YAxis.FontSize = 13;
    ax.FontName='Arial'
    xlim([5 50])

    clear var*
    clear idx*
    clear stat
    clear ii
    clear ia

    if thisP==1
        ylim ([0 0.35])
    elseif thisP==2
        ylim ([0 0.15])
    elseif thisP==3
        ylim ([0 0.15])
    end

 
end

clear a
clear alpha
clear ind*

% --------- Granger --------

strAA={'M1m-PMVm', 'PMVm-M1m', 'M1m-M3', 'M3-M1m', 'M3-PMVm', 'PMVm-M3'};
strG ='Granger';
lw=8;

figure(1);

for pair=1:6

    load([strG,'_FX_',strAA{pair},'.mat'])
    var1; % mean coherence LS (?)
    var2; % mean coherence T
    stat.freq;


    disp([strG,'_FX_',strAA{pair},'.mat'])
    disp(['Stats Granger alpha =', num2str(stat.cfg.alpha)])
    disp(['Stats Granger nperm =', num2str(stat.cfg.numrandomization)])

    %% Get 2 first max of Granger 
    [pks, locs] = findpeaks(var1); [ii,ia]=sort(var1(locs),'descend');
    FreqMaxPk = stat.freq(locs(ia(1:5))); GMaxPk = var1(locs(ia(1:5)));
    disp('Granger MaxPkk LS')
    disp([FreqMaxPk' GMaxPk'])

    [pks, locs] = findpeaks(var2); [ii,ia]=sort(var2(locs),'descend');
    FreqMaxPk = stat.freq(locs(ia(1:5))); PPCMaxPk = var2(locs(ia(1:5)));
    disp('Granger MaxPkk T')
    disp([FreqMaxPk' GMaxPk'])

   
    if pair==1 || pair==2
        subplot(2,3,4); hold on;
    elseif pair==3 || pair==4
        subplot(2,3,5); hold on;
    elseif pair==5 || pair ==6
        subplot(2,3,6); hold on;
    end

    temp=rem(pair,2);
    disp([strAA(pair)]);
    
    if pair==1 || pair==4 || pair==5  % for even numbers area1--> area2
        plot(stat.freq,var1,'-','linewidth', 2, 'Color',[0.9882    0.9255    0.0196]); % yellow
        plot(stat.freq,var2,'-','linewidth', 2, 'Color',[0.8510    0.3255    0.0980]); % red
    elseif pair==2 || pair==3 || pair==6 % for pair numbers area2--> area1
        plot(stat.freq,var1,'--','linewidth', 2, 'Color',[0.9882    0.9255    0.0196])
        plot(stat.freq,var2,'--r','linewidth', 2, 'Color',[0.8510    0.3255    0.0980])
        xlabel('Freq (Hz)','fontsize', 14)
        ylabel('Granger','fontsize', 14)
    end

    if isfield(stat,'negclusters')==1 && isempty(stat.negclusters)~=1
        for i=1:size(stat.negclusters,2)
            alpha=stat.negclusters(i).prob;
            if alpha<=pvalue
                alpha
                stat.freq(stat.negclusterslabelmat==i)
                disp(['T>L i=', num2str(i)]) 
                [y,idx1]=ismember(stat.freq(stat.negclusterslabelmat==i), stat.freq);
                plot(stat.freq(stat.negclusterslabelmat==i),var1(idx1), 'Color',[0.9882    0.9255    0.0196, 0.9],'LineWidth', lw); % yellow
                plot(stat.freq(stat.negclusterslabelmat==i),var2(idx1), 'Color',[0.8510    0.3255    0.0980, 0.9],'LineWidth',lw); % red
                pause(1)
                clear y
                clear idx1
                clear alpha
            end
        end
    end
    
    xlim([5 50])
    ax = gca(gcf)
    ax.XAxis.FontSize = 13;
    ax.YAxis.FontSize = 13;
    ax.FontName='Arial'
    
end


x=[4 5 6];
thisP = [2 4 6];

for i=1:3

    subplot(2,3,x(i))

    if x(i)==4 || x(i)==6

        ylim([0 0.4])
        plot([25 28], [0.37 0.37],'-','LineWidth', 1.5, 'Color',[0.9882    0.9255    0.0196]) % change for yellow 
        text(30, 0.37,['LS-',strAA{thisP(i)-1}]);

        plot([25 28], [0.35 0.35],'-r', 'LineWidth', 1.5, 'Color',[0.8510    0.3255    0.0980]) % change for smoked red 
        text(30, 0.35,['T-',strAA{thisP(i)-1}]);

        plot([25 28], [0.30 0.30],'--','LineWidth', 1.5, 'Color',[0.9882    0.9255    0.0196]) % change for yellow 
        text(30, 0.30,['LS-',strAA{thisP(i)}])

        plot([25 28], [0.27 0.27],'--r','LineWidth', 1.5, 'Color',[0.8510    0.3255    0.0980]) % change for smoked red 
        text(30, 0.27,['T-',strAA{thisP(i)}])
        ylim([0 0.4])

    elseif x(i)==5
        ylim([0 0.7])
        plot([25 28], [0.65 0.65],'--g','LineWidth', 1.5, 'Color',[0.9882    0.9255    0.0196]) % yellow 
        text(30, 0.65,['LS-',strAA{thisP(i)-1}]);

        plot([25 28], [0.60 0.60],'--r', 'LineWidth', 1.5, 'Color',[0.8510    0.3255    0.0980]) % smoked red
        text(30, 0.60,['T-',strAA{thisP(i)-1}]);

        plot([25 28], [0.55 0.55],'-g','LineWidth', 1.5, 'Color',[0.9882    0.9255    0.0196]) % yellow 
        text(30, 0.55,['LS-',strAA{thisP(i)}])

        plot([25 28], [0.50 0.50],'-r','LineWidth', 1.5, 'Color',[0.8510    0.3255    0.0980]) % smoked red 
        text(30, 0.50,['T-',strAA{thisP(i)}])
        ylim([0 0.7])

        
    end 

end

PrintOn=0;

% set(gcf,'Position', get (0,'Screensize')); % do not do Full Screen
% so that the size matches in Illustrator. In Ilustrator 30% opasity for
% Red, yellow opasity to 50%% 

if PrintOn==1
    print(gcf,'-depsc', 'Fig3_PPC-Granger_v1f')
    print('Fig3_PPC-Granger_v1f','-dsvg')
    savefig(['Fig3_PPC-Granger_v1f', strM,'-Granger.fig'])
    close
end





