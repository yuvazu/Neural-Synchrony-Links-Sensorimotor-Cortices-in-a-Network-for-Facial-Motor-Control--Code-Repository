% MakeFigures_rebuttal2
% This Mfile plots the eLFP for all receiver areas given a sender area (S1,
% M1, PMV, M3). 
% And plots the connectivity matrix (6 ms delta |pk-tr|) 
% for all combinations of sender(s1,m1,pmv,m3) vs receiver (s1, m1, pmv, m3)
% used for rebuttal 
% nov 2025
% yvz
clear all
close all

Storedir ='/Users/yuriria/Documents/MATLAB/Fx_Con_Analysis/LFP_related/Cleaned_eLFP/test_241127/Clean_dPT_AUC_Fixed6mstoVarOffset/SameOffset_afterStim';
cd(Storedir) 
files = dir('*Clean_dPT_AUC_sameToffset.fig');


 for i= 1:length(files)
     fname = files(i)
     openfig(fname.name);
     gcf;
    % Get all axes handles in the current figure
    ax_handles = findall(gcf, 'Type', 'axes');

    % Set properties for all axes at once
    set(ax_handles, 'XLim', [-50 105]);
    set(ax_handles, 'YLim', [-20 20]);

    if strncmp (fname(1).name,'M1m',3)==1
        subplot(2,2,3)
        set(gca,'Ylim',[-5 5])
    end
    
    % Or set different properties
    set(ax_handles, 'FontSize', 14,'FontName', 'Arial','TickLength', [0.02 0.02]);
    %set(gcf, 'Position', get(0, 'Screensize'));
    print(gcf, '-dsvg', fname.name(1:end-4));
    pause(3)
    close
 end
  
 
  

%  % Confusion matrix based on fixed onset 6 ms
Storedir ='/Users/yuriria/Documents/MATLAB/Fx_Con_Analysis/LFP_related/Cleaned_eLFP/test_241127/Clean_dPT_AUC_Fixed6mstoVarOffset/SameOffset_afterStim';
cd(Storedir) 
openfig('ConfusionMatrix_dPK-Tr_6FixedtoVar.fig')
% Set data aspect ratio
pbaspect([1 1 1]);  % Plot box aspect ratio

%set(gcf, 'Position', get(0, 'Screensize'));
print(gcf, '-depsc', 'ConfusionMatrix_dPK-Tr_Fixed6mstoVarOffset');
print(gcf, '-dtiff', 'ConfusionMatrix_dPK-Tr_Fixed6mstoVarOffset');