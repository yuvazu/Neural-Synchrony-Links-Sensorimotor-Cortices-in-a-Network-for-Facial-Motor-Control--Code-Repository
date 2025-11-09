% This mfile was used to  find the peak coupling frequency among the coupled areas
% (Table 1, 2 monkeys, all behaviors)
%  You will need to edit the paths to the dir data.
% for each figure PPC Data in the following dir
% Figure3/Data_PPC_Fig3
% Supple2/Data_2M_Mov_ChewT
% Supple2/Data_2M_Mov_ChewLS
% updated 2025


close all
clear all

Cmb = {'TLS', 'CT','CLS'}; % Behavior Combination
thisC =1; % edit to get the Combination of Behavior you want to measure the PPC freq peaks
thispair = 1; % edit to get the Area Pair you want to measure the PPC freq peaks


if strcmp(Cmb(thisC),'TLS')==1
    %% Threat vs LS
    cd('Figure3/Data_PPC_Fig3')
    AreaPair={'M1m-PMVm', 'M1m-M3', 'M3-PMVm'}; % original set of areas in LS/T
    %AreaPair = {'M1m-S1', 'M3-S1', 'PMVm-S1'}
    thisA=AreaPair{thispair};
    a=dir(['PPC_FX_',thisA,'*.mat']);
elseif strcmp(Cmb(thisC),'CT')==1
    %% Chew vs Threat
    cd('Supple2/Data_Supple2_PPC_CT ')
    AreaPair={'M1m-PMVm', 'M1m-M3', 'M3-PMVm'}; % original set of areas in LS/T%AreaPair = {'M1m-S1', 'M3-S1', 'PMVm-S1'}
    %AreaPair = {'M1m-S1', 'M3-S1', 'PMVm-S1'}
    thisA=AreaPair{thispair};
    a=dir(['PPC_FX_',thisA,'*.mat']);
elseif strcmp(Cmb(thisC),'CLS')==1
    %% Chew vs LS
    cd('Supple2/Data_Supple2_PPC_CLS')%
    AreaPair={'M1m-PMV', 'M1m-M3', 'PMV-M3'}; % or
    %AreaPair = {'M1m-S1', 'M3-S1', 'PMV-S1'}
     thisA=AreaPair{thispair};
    a=dir(['PPC_FXC_',thisA,'*.mat']);
end



load(a.name)
disp(['filename = ', a.name])
% 
disp(thisA)

% label d1 (T/LS, ) 
[pks, locs] = findpeaks(var1); [ii,ia]=sort(var1(locs),'descend');
FreqMaxPk = stat.freq(locs(ia(1:5))); PPCMaxPk = var1(locs(ia(1:5)));
disp(labeld1)
disp([FreqMaxPk' PPCMaxPk'])

% label d2 (Chew) 
[pks, locs] = findpeaks(var2); [ii,ia]=sort(var2(locs),'descend');
FreqMaxPk = stat.freq(locs(ia(1:5))); PPCMaxPk = var2(locs(ia(1:5)));
disp(labeld2)
disp([FreqMaxPk' PPCMaxPk'])
figure; hold on;
plot(stat.freq,var1,'-m', 'LineWidth',4);
plot(stat.freq,var2,'-k', 'LineWidth',4);
axis tight; legend(labeld1,labeld2)
title(['PPC_FX_',thisA])


