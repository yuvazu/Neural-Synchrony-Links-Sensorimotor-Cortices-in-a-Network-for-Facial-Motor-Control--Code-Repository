%% Plot_ConfusionMatrix_wFixedOnset6ms
% this mfiles plots the confusion matrix depicting the connectivity btw
% sender and receiver areas measured by the delta(pl-thr) using
% for the figure article we used SameOffset_afterStim (6ms) 
% based on the original script (Plot_ConfusionMatrix)
% for the onsets are stored (feb 2025) 
% Generating Figure 2 (rebuttal2)
% oct/2025

TargetDir='/Users/yuriria/Documents/MATLAB/Fx_Con_Analysis/LFP_related/Cleaned_eLFP/test_241127/Clean_dPT_AUC_Fixed6mstoVarOffset/SameOffset_afterStim';
cd(TargetDir)

clear all;
close all;

toPrint=0;

m=nan(4,4);
areas = {'S1', 'M1m', 'PMVm', 'M3'};

load('S1Clean_dPT_AUC_Times_sameToffset.mat');
store_dPT
m(2:end,1)= store_dPT;

load('M1mClean_dPT_AUC_Times_sameToffset.mat')
store_dPT
m(1,2)= store_dPT(1);
m(3:4,2)= store_dPT(2:end);

load('PMVmClean_dPT_AUC_Times_sameToffset.mat')
store_dPT
m(1:2,3)= store_dPT(1:2);
m(4,3)= store_dPT(end);

load('M3Clean_dPT_AUC_Times_sameToffset.mat')
store_dPT
m(1:3,4)= store_dPT;


figure(100);
imagesc(m)
cmap = parula(256);
hold on

[row,col] = find(isnan(m));
for i = 1:length(row)
    rectangle('Position',[col(i)-0.5,row(i)-0.5,1,1],'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
end

colorbar
minVal = 0;
maxVal = max(m(:));
clim([minVal maxVal])

xticks(1:length(areas)+1); set(gca, 'XtickLabels', areas)
yticks(1:length(areas)+1); set(gca, 'YtickLabels', areas)
a=colorbar; a.Label.String=['\delta', '|Pk-Tr|', '  (au)']; a.FontSize=14;
xlabel('Sender of ICS','FontSize',14);
ylabel('Receiver Area','FontSize',14);
set(gca,'FontSize',14)

if toPrint==1
    title (['fixed Onset 6ms -var', '\delta', '|Pk-Tr|', '  (au)'])
    savefig('ConfusionMatrix_dPK-Tr_6FixedtoVar')
    save('ConfusionMatrix_dPK-Tr_6FixedtoVar', 'm', 'areas')

    set(gcf, 'Position', get(0, 'Screensize'));
    print(gcf,'-depsc','ConfusionMatrix_dPK-Tr_6FixedtoVar')
end