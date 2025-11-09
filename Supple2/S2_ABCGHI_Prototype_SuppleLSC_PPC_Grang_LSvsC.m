function S2_ABCGHI_Prototype_SuppleLSC_PPC_Grang_LSvsC
% this function loads the mean Coh/PPC across all sessions for LS &
% Chews and plots them (Ntrials=2Monkeys, with auto chew detections).
% A graph (2X3) is created showing the mean
% PPC/Coh for both behaviors for medial-lateral and lateral-lateral motor
% areas, S1 is excluded from this graph -it has its own mfile!
% Stats corrected pvalue<=0.05
% alpha=0.05, nperm=1000. var(1): Threats, var(2):Chew
% the results were done using TwoMonkeys_stats_mod_Chew_vs_T (11 sessions , 2 monkeys)
% Run and stored % July 17th 2025
clear all
close all

% modify TargetDir
TargetDir='/Users/yuriria/Documents/MATLAB/Thor_Analysis/LFP_analysis/Coh/Barney_FromEli/ElieSCode/figs/TwoMonkeys/2M_Mov_ChewLS/'
cd(TargetDir)

strM='PPC';
lw=12;
pvalue=0.05; % for plotting purposes



% -------- M1-PMV PPC   -----------
strAA={'M1m-M3', 'PMV-M3', 'M1m-PMV'};

for thisP=1:3
    % Load PPC for this pair of areas
    load([strM,'_FXC_', strAA{thisP},'.mat']);
    disp([strM,'_FXC_', strAA{thisP}])
    var1; % LS
    var2; %  Chew

    [pks, locs] = findpeaks(var1); [ii,ia]=sort(var1(locs),'descend');
    FreqMaxPk = stat.freq(locs(ia(1:5))); PPCMaxPk = var1(locs(ia(1:5)));
    disp('LS')
    disp([FreqMaxPk' PPCMaxPk'])

    [pks, locs] = findpeaks(var2); [ii,ia]=sort(var2(locs),'descend');
    FreqMaxPk = stat.freq(locs(ia(1:5))); PPCMaxPk = var2(locs(ia(1:5)));
    disp('C')
    disp([FreqMaxPk' PPCMaxPk'])

    hold off

    stat.freq;
    disp([strM,'_LSC_', strAA{thisP}])
    disp(['Stats PPC alpha =', num2str(stat.cfg.alpha)])
    disp(['Stats PPC nperm =', num2str(stat.cfg.numrandomization)])

    figure(1); subplot(2,3,thisP); hold on;
    plot(stat.freq,var1,'-','linewidth', 2,'Color', [0.9882    0.9255    0.0196]) % % yellow
    plot(stat.freq,var2,'-','linewidth', 2,'Color', [0    0.4471    0.7412]) % blue   (Chew)
   
    % % Assess Significant Negative Clusters C>T
    if isfield(stat,'negclusters')==1 && isempty(stat.negclusters)~=1
        for i=1:size(stat.negclusters,2)
            alpha=stat.negclusters(i).prob;
            if alpha<=pvalue % if this cluster is significant
                alpha;
                stat.freq(stat.negclusterslabelmat==i);
                disp('LS>C')
                [y,idx1]=ismember(stat.freq(stat.negclusterslabelmat==i), stat.freq);
                plot(stat.freq(stat.negclusterslabelmat==i),var1(idx1), 'Color',[0.9882    0.9255    0.0196, 0.6],'LineWidth', lw); % red
                plot(stat.freq(stat.negclusterslabelmat==i),var2(idx1), 'Color',[0    0.4471    0.7412, 0.6],'LineWidth',lw); % blue
                %pause
                clear y
                clear idx1
            end
        end
    end
    clear alpha

    % Assess Significant Negative Clusters T>C
    if isfield(stat,'posclusters')==1 && isempty(stat.posclusters)~=1
        for i=1:size(stat.posclusters,2)
            alpha=stat.posclusters(i).prob;
            if alpha<=pvalue % if this cluster is significant
                alpha;
                stat.freq(stat.posclusterslabelmat==i);
                disp('C>LS')
                [y,idx1]=ismember(stat.freq(stat.posclusterslabelmat==i), stat.freq);
                plot(stat.freq(stat.posclusterslabelmat==i),var1(idx1), 'Color',[0.9882    0.9255    0.0196, 0.6],'LineWidth', lw); % red
                plot(stat.freq(stat.posclusterslabelmat==i),var2(idx1), 'Color',[0    0.4471    0.7412, 0.6],'LineWidth',lw); %blue
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

    disp('stat.freq p<0.5')
    stat.freq(stat.prob<=0.05)

    clear var*
    clear idx*
    clear stat
    clear ii
    clear ia
    % 
     if thisP==1; % M1-PMVm
         ylim ([-0.01 0.30])
     elseif thisP==2  % M1-M3
         ylim ([-0.01 0.30])
     elseif thisP==3 % M3-PMVm
         ylim ([-0.01 0.4])
     end

end



clear a
clear alpha
clear ind*

% --------- Granger --------

strAA={'M1m-M3', 'M3-M1m','M3-PMV', 'PMV-M3','M1m-PMV', 'PMV-M1m'}


strG ='Granger';
lw=12;

figure(1);

for pair=1:6

    load([strG,'_FXC_',strAA{pair},'.mat'])
    var1; % mean PPC T
    var2; % mean PPC C
    stat.freq;


    disp([strG,'_FXC_',strAA{pair},'.mat'])
    disp(['Stats Granger alpha =', num2str(stat.cfg.alpha)])
    disp(['Stats Granger nperm =', num2str(stat.cfg.numrandomization)])

    %% Get 2 first max of Granger 
    [pks, locs] = findpeaks(var1); [ii,ia]=sort(var1(locs),'descend');
    FreqMaxPk = stat.freq(locs(ia(1:length(pks)))); GMaxPk = var1(locs(ia(1:length(pks))));
    disp('Granger MaxPkk LS')
    disp([FreqMaxPk' GMaxPk'])

    [pks, locs] = findpeaks(var2); [ii,ia]=sort(var2(locs),'descend');
    FreqMaxPk = stat.freq(locs(ia(1:length(pks)))); GMaxPk = var2(locs(ia(1:length(pks))));
    disp('Granger MaxPkk Chews')
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
    
    if pair==2 || pair==3 || pair==5  % for even numbers area1--> area2
        plot(stat.freq,var1,'-','linewidth', 2, 'Color',[0.9882    0.9255    0.0196]); % yellow
        plot(stat.freq,var2,'-','linewidth', 2, 'Color',[0    0.4471    0.7412]); % blue 
    elseif pair==1 || pair==4 || pair==6 % for pair numbers area2--> area1
        plot(stat.freq,var1,'--','linewidth', 2, 'Color',[0.9882    0.9255    0.0196]); % yellow
        plot(stat.freq,var2,'--r','linewidth', 2, 'Color',[0    0.4471    0.7412]); %blue 
        xlabel('Freq (Hz)','fontsize', 14)
        ylabel('Granger','fontsize', 14)
    end

    if isfield(stat,'negclusters')==1 && isempty(stat.negclusters)~=1
        for i=1:size(stat.negclusters,2)
            alpha=stat.negclusters(i).prob;
            if alpha<=pvalue
                alpha
                stat.freq(stat.negclusterslabelmat==i)
                disp(['C>T i=', num2str(i)])
                [y,idx1]=ismember(stat.freq(stat.negclusterslabelmat==i), stat.freq);
                plot(stat.freq(stat.negclusterslabelmat==i),var1(idx1), 'Color',[0.9882    0.9255    0.0196, 0.6],'LineWidth', lw); % yellow
                plot(stat.freq(stat.negclusterslabelmat==i),var2(idx1), 'Color',[0         0.4471    0.7412, 0.6],'LineWidth', lw); % blue 
                clear y
                clear idx1
                clear alpha
            end
        end
    end

     if isfield(stat,'posclusters')==1 && isempty(stat.posclusters)~=1
        for i=1:size(stat.posclusters,2)
            alpha=stat.posclusters(i).prob;
            if alpha<=pvalue
                alpha
                stat.freq(stat.posclusterslabelmat==i)
                disp(['T>C i=', num2str(i)]) 
                [y,idx1]=ismember(stat.freq(stat.posclusterslabelmat==i), stat.freq);
                plot(stat.freq(stat.posclusterslabelmat==i),var1(idx1), 'Color',[0.9882    0.9255    0.0196, 0.6],'LineWidth', lw); % yellow
                plot(stat.freq(stat.posclusterslabelmat==i),var2(idx1), 'Color',[0         0.4471    0.7412, 0.6],'LineWidth',lw); % blue
                %pause
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
    ax.FontName='Arial';

    
end

x=[4 5 6];
thisP = [2 4 6];



for i=1:3

    subplot(2,3,x(i))

    if x(i)==4 

        plot([25 28], [0.45 0.45],'--','LineWidth', 1.5, 'Color',[0.9882    0.9255    0.0196]) % blue
        text(30, 0.45,['LS-',strAA{thisP(i)-1}]);

        plot([25 28], [0.40 0.40],'--', 'LineWidth', 1.5, 'Color',[0         0.4471    0.7412]) % change for blue
        text(30, 0.40,['C-',strAA{thisP(i)-1}]);

        plot([25 28], [0.35 0.35],'-','LineWidth', 1.5, 'Color',[0.9882    0.9255    0.0196]) % change for red
        text(30, 0.35,['LS-',strAA{thisP(i)}])

        plot([25 28], [0.30 0.30],'-','LineWidth', 1.5, 'Color',[0         0.4471    0.7412]) % change for blue
        text(30, 0.30,['C-',strAA{thisP(i)}])
        ylim([0 0.6])

    elseif x(i)==5 || x(i)==6
       
        plot([25 28], [0.45 0.45],'-','LineWidth', 1.5, 'Color',[0.9882    0.9255    0.0196]) % yellow 
        text(30, 0.45,['LS-',strAA{thisP(i)-1}]);

        plot([25 28], [0.40 0.40],'-', 'LineWidth', 1.5, 'Color',[0         0.4471    0.7412]) % smoked red
        text(30, 0.40,['C-',strAA{thisP(i)-1}]);

        plot([25 28], [0.35 0.35],'--','LineWidth', 1.5, 'Color',[0.9882    0.9255    0.0196]) % yellow 
        text(30, 0.35,['LS-',strAA{thisP(i)}])

        plot([25 28], [0.30 0.30],'--','LineWidth', 1.5, 'Color',[0         0.4471    0.7412]) % smoked red 
        text(30, 0.30,['C-',strAA{thisP(i)}])
        ylim([0 0.60])

        
    end 

end

PrintOn=0;

% do not do Full Screen
% so that the size matches in illustrator. In Ilustrator 30% oposity for
% Red, and then copy the yellow and opasity to 50%

% 
if PrintOn==1
    print(gcf,'-depsc', 'SuppleCLS_PPC-Granger_r1')
    print('SuppleCLS_PPC-Granger_r1','-dsvg')
    savefig('SuppleCLS_PPC-Granger_r1.fig')
    close
end





