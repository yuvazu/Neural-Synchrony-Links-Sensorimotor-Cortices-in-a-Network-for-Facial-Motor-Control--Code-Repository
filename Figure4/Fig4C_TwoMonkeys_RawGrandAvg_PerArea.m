%% Fig4C_TwoMonkeys_RawGrandAvg_PerArea
% this mfile plots the Avg Raw power for all 2 (or 3)  behaviors within an Area
% with sem of all the sessions. The idea is to see in which areas the
% mean oscillatory activity really differs. 
% this to see if within an area we actually have differences in the
% spectrum depending on the behavior. 
% the power is calculated 1 second before and after the onset of the
% movement. 
% Also plots Supplementary 3A-D (w Chew) and Supplementary 4A-D for high
% freq+Chew
% For plotting Figure S3 low freq (5-40 Hz), uncomment line 22, and comment line 23
% For plotting Figure S4 high freq ( 40-120 Hz), uncomment line 23, and comment line 22. 


% Updated April 30th, 2025

clear all
close all
restoredefaultpath
addpath('/Users/yuriria/Documents/MATLAB/fieldtrip/fieldtrip-20170101')
addpath('/Users/yuriria/Documents/MATLAB/Thor_Analysis/func')
ft_defaults

% Set the frequency range
%freqrange ='lowfreq'; % 5-40 Hz
freqrange = 'highfreq'; % 40-120 Hz 
PrintOn=0;

if strcmp(freqrange,'lowfreq')
    fROI='all';
    %fROI=[5 10]
    fst='hann';
else strcmp(freqrange,'highfreq')
    fROI = 'all';
    fst='dpss';
end


% --- Set this to the Directory where the Data Files are located(raw power 2M collapsed) ---- %
TargetDir='/Volumes/Seagate/TwoMonkeys_TB/TFR/TFR_RawPower_AvgSess'; % local directory
%TargetDir'/Users/yuriria/Dropbox/FreiwaldLabMine/code/A_GitHub/Figure_4/Data_Figure4'; % DB directory



AllAreas = {'S1','M1m','PMV','M3'}
beh = {'Threat', 'LipSmacking'} % 
%beh = {'Threat', 'LipSmacking',Chew} %  uncomment this line to plot FigS3A-D

for a=1:length(AllAreas)
    thisArea=AllAreas{a};
    cd(TargetDir)
    % loading all Raw Power for Threat (Nsessions X 1 X time X Freq)
    disp(['-----AvgPowerPerSess_', beh{1},'_',thisArea,'---'])
    load([TargetDir,'/AvgPowerPerSess_2M_', beh{1},'_',thisArea,'_',fst], 'freq2M')
    freqT=freq2M;

    pT = squeeze(log10(freqT.powspctrm));
    nT = size(pT,1); % nSessions
    mpT = (squeeze(nanmean(squeeze(nanmean(pT,1)),2))); % mean across sessions & time
    spT = std(squeeze(nanmean(pT,3)),1)/sqrt(nT); % sem across sessions

    % loading all Raw Power for LipSmacking
    disp(['-----AvgPowerPerSess_', beh{2},'_',thisArea,'---'])
    load([TargetDir,'/AvgPowerPerSess_2M_', beh{2},'_',thisArea,'_',fst], 'freq2M')
    freqLS=freq2M;
    pLS = squeeze(log10(freqLS.powspctrm));
    nLS = size(pLS,1); % nSessions
    mpLS=squeeze(nanmean(squeeze(nanmean(pLS,1)),2));
    spLS = std(squeeze(nanmean(pLS,3)),1)/sqrt(nLS); % sem across sessions

    if length(beh)>2

        disp(['-----AvgPowerPerSess_', beh{3},'_',thisArea,'---'])
        load([TargetDir,'/AvgPowerPerSess_2M_', beh{3},'_',thisArea,'_',fst], 'freq2M')
        freqC=freq2M;
        pC = squeeze(log10(freqC.powspctrm));
        nC = size(pC,1)
        mpC = squeeze(nanmean(squeeze(nanmean(pC,1)),2));
        spC = std(squeeze(nanmean(pC,3)),1)/sqrt(nC); % sem across sessions
    end

    figure; clf; hold on % Raw Power
    plot(freqT.freq, mpT, 'r','LineWidth', 3.5); % Threat
    plot(freqT.freq, mpT'+spT, '--r','LineWidth', 1.0); % Threat
    plot(freqT.freq, mpT'-spT, '--r','LineWidth', 1.0); % Threat

    plot(freqLS.freq, mpLS, '-','LineWidth', 3.5,'Color', [0.9020    0.9020    0.2078]) % LS
    plot(freqLS.freq, mpLS'+spLS, '--','LineWidth', 1.05,'Color', [0.9020    0.9020    0.2078]) % LS
    plot(freqLS.freq, mpLS'-spLS, '--','LineWidth', 1.05, 'Color', [0.9020    0.9020    0.2078]) % LS

    if length(beh)>2

    plot(freqLS.freq, mpC, '-','LineWidth', 3.5, 'Color', [0    0.4471    0.7412]) % C
    plot(freqLS.freq, mpC'+spC, '--','LineWidth', 1.05,'Color', [0    0.4471    0.7412]) % C
    plot(freqLS.freq, mpC'-spC, '--','LineWidth', 1.05, 'Color', [0    0.4471    0.7412]) % C
    end
    axis tight

    ax=axis;



    if strcmp(fst,'hann')
        xlim([4 41])
        dy=0.04
        mUp=ax(4)-dy;
        ic = dy;
        ic2 = ic*2
        plot([30 33], [mUp mUp],'-r','LineWidth',2.0);
        text( 35, mUp, 'threat','FontSize',14, 'FontName','Arial')
        plot([30 33], [mUp-ic mUp-ic],'-','LineWidth',2.0, 'Color', [0.9020    0.9020    0.2078]);
        text( 35, mUp-ic, 'lipsmacking','FontSize',14, 'FontName','Arial')
        if length(beh)>2
            plot([30 33], [mUp-ic2 mUp-ic2],'-b','LineWidth',2.0, 'Color', [0    0.4471    0.7412]);
            text( 35, mUp-ic2, 'chew',  'FontSize',14, 'FontName','Arial');
        end
    elseif strcmp(fst,'dpss')

        mUp=ax(4)-.06;
        ic = .09;
        ic2 = ic*2
        plot([100 103], [mUp mUp],'-r','LineWidth',2.0);
        text( 105, mUp, 'threat','FontSize',16, 'FontName','Arial')
        plot([100 103], [mUp-ic mUp-ic],'-','LineWidth',2.0, 'Color',[0.9020    0.9020    0.2078]);
        text( 105, mUp-ic, 'lipsmacking','FontSize',16, 'FontName','Arial')
        if length(beh)>2
            plot([100 103], [mUp-ic2 mUp-ic2],'-b','LineWidth',2.0, [0    0.4471    0.7412]);
            text( 105, mUp-ic2, 'chew','FontSize',16, 'FontName','Arial');
        end
    end



    xlabel('Frequency (Hz)', 'FontSize',18, 'FontName','Arial')
    ylabel('Power (db)', 'FontSize',18, 'FontName','Arial')
    title (thisArea, 'FontSize',24,'FontName','Arial');

    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Arial','fontsize',18)

    b = get(gca,'YTickLabel');
    set(gca,'YTickLabel',b,'FontName','Arial','fontsize',18)

    if PrintOn==1
        cd('/Volumes/Seagate/TwoMonkeys_TB/TFR/TFR_RawPower_AvgSess/2M_RawPower_vs_Freq/AllBeh')
        savefig(['AvgRawPower_AllBeh_',thisArea, '_',fst])
        print(gcf,'-dpng',['AvgRawPower_AllBeh_',thisArea, '_',fst])
    end


end



       

     






















        


