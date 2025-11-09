function Fig4AB_Granger_Connectivity_FourNodes
%% Fig4AB_Granger_Connectivity_FourNodes
% this function creates Figure 4A or 4B Max Granger in Node Graphs (D-F and G-I)
% It reads & plots:
% 1) Mean Granger values for freq range (4-50 Hz) per behavior
% fig1 (LS Lat--> Medial), fig2 LS Medial-->Lat, fig 3 (LS Lat--Lat) 
% fig 1001 (T Lat--> Medial), fig 1002 (T Medial --> Lat),  fig 1003 (T Lat--Lat)

% Medial (M3) to Lateral Areas ( M1, PMV, S1)
% Lateral (M1, PMV, S1) to Medial Area (M3)
% Lateral-Lateral (M1--> PMV, M1--> S1, S1--> PMV) 
% Lateral-Lateral (PMV-->M1, PMV--> S1, S1-->M1) 


% 2)Node graphs are created for each of this combinations
 % fig10-40 (raw node graphs)
 % fig 11 Node Graph: Medial (M3) to Lateral Areas ( M1, PMV, S1)
 % fig 21 Node Graph: Lateral (S1, PMV, M1) to Medial (M3)
 % fig 31 Node Graph: Lateral (M1, S1) to Lateral (PMV)
 % fig 41 Node Graph: Lateral (PMV,S1) to Lateral (M1) 

% for Illustrator print('-dsvg')
% local file in Mac: Granger_Connectivity_FourNodes. 
%Target=/Users/yuriria/Dropbox/FreiwaldLabMine/code/A_GitHub/Figure_3/Data_PPC_Fig3
% updated April 30th, 2025

clear all
close all

%% Modify Path to targets the Data! 
%TargetDir='/Users/yuriria/Documents/MATLAB/Thor_Analysis/LFP_analysis/Coh/Barney_FromEli/ElieSCode/figs/TwoMonkeys/2M_Movement'; % modify here for path
TargetDir='/Users/yuriria/Dropbox/FreiwaldLabMine/code/A_GitHub/Figure_3/Data_PPC_Fig3';
cd(TargetDir)
strG ='Granger';
PrintOn=0;
PlotG_freq=0;

%% Granger Lateral Areas to Medial ones

c= {'LM','ML','LL'}; % this defines lateral-medial, medial to lateral, lateral to lateral FLOWS

 lm_l=[];
 lm_t=[];
 
 ml_l=[];
 ml_t=[];
            
 ll_l=[];
 ll_t=[];


for i=1:length(c)

    if strcmp(c(i),'LM')
        strAA={ 'M1m-M3','PMVm-M3', 'S1-M3'}; % lateral to medial flow 
        strT = 'Latetal to Medial Areas';
        nF = 1;
        strLA={ 'LS M1->M3', 'LS PMV->M3', 'LS S1->M3'};
        strLA1={ 'T M1->M3', 'T PMVm->M3', 'T S1->M3'};

    elseif strcmp(c(i),'ML')
        strAA={ 'M3-M1m','M3-PMVm', 'M3-S1'}; % medial to lateral flow 
        strT = 'Medial to Lateral Areas';
        nF = 2;
        strLA={ 'LS M3->M1', 'LS M3->PMV', 'LS M3->S1'};
        strLA1={ 'T M3->M1', 'T M3->PMV', 'T M3->S1'};

    elseif strcmp (c(i),'LL')
        strAA={ 'M1m-PMVm', 'PMVm-M1m', 'M1m-S1','S1-M1m','PMVm-S1', 'S1-PMVm'}; % within lateral areas 
        strT = ' Within Lateral Areas';
        nF =3;
        %strLA={ 'LS M1->PMV', 'T M1->PMV', 'LS PMV-M1', 'T PMV->M1', 'LS M1->S1', 'T M1->S1', 'LS S1->M1', 'T S1->M1', 'LS PMV->S1', 'T PMV->S1',...
            %'LS S1->PMV', 'T S1-PMV'}

        strLA={ 'LS M1->PMV', 'LS PMV-M1', 'LS M1->S1', 'LS S1->M1','LS PMV->S1', 'LS S1->PMV'}
        strLA1={ 'T M1->PMV', 'T PMV->M1', 'T M1->S1', 'T S1->M1',  'T PMV->S1', 'T S1-PMV'}


    end

    lw=8;
    c1 = [0.9882    0.9255    0.0196; 0.9294    0.6941    0.1255; 1.0000    0.6902    0.0706; 0.7686    0.7216    0.0471;  0.6000    0.6000    0.1137; 0.7176    0.2745    1.0000]; % yellow
    c2 = [0.8510    0.3255    0.0980; 1 0 0;  0.6353    0.0784    0.1843; 0.5882    0.3490    0.3490; 0.9804    0.8039    0.8039;  1.0000    0.0745    0.6510]; % red

    for pair=1:length(strAA)
        load([strG,'_FX_',strAA{pair},'.mat']);
        var1; % mean Granger LS (?)
        var2; % mean Granger T
        stat.freq;
        idx = find(stat.freq>6 & stat.freq<31);
        freq = stat.freq(idx);

        % store the Freq and Magnitude of the Maximum Flow 
        [m1,ix1] = max(var1(idx));
        [m2,ix2] = max(var2(idx));

        if strcmp(c(i),'LM')==1
           
            lm_l(pair,:) = [m1 freq(ix1)];
            lm_t(pair,:) = [m2 freq(ix2)];

            disp('LM')
        elseif strcmp(c(i),'ML')==1
            ml_l(pair,:) = [m1 freq(ix1)];
            ml_t(pair,:) = [m2 freq(ix2)];
            disp('ML')

        elseif strcmp(c(i),'LL')==1
            i
            pair
            ll_l(pair,:) = [m1 freq(ix1)];
            ll_t(pair,:) = [m2 freq(ix2)];
            disp('LL')

        end
        
        %% plotting Granger as a function of Frequency
        if PlotG_freq==1
            figure(nF); hold on;
            plot(stat.freq(idx),var1(idx),'Color',  c1(pair,:), 'LineWidth', 5); % LS
            figure(nF+1000); hold on;
            plot(stat.freq(idx),var2(idx),'Color', c2(pair,:), 'LineWidth', 5)% Threats
        end
        clear var*
        clear stat*
    end


    if PlotG_freq==1
    axis tight
    figure(nF)
    legend(strLA, 'fontsize', 14)
    title([strT, ' LS'],'fontsize', 20);
    xlabel('Freq (Hz)', 'fontsize', 16)
    ylabel('Granger', 'fontsize', 16)

    figure(nF+1000)
    legend(strLA1, 'fontsize', 14)
    title([strT, ' T'],'fontsize', 20);
    xlabel('Freq (Hz)', 'fontsize', 16)
    ylabel('Granger', 'fontsize', 16)
    end
end

%% ----- Set Behavior beh='T' for threats, 'LS' for lipsmacks ------ %%
beh='LS'

% lateral to medial strLA={'M1->M3', 'PMV->M3', 'S1->M3'};
if strcmp(beh,'T')
    lm= lm_t;
elseif strcmp(beh,'LS')
    lm = lm_l;
end
lm;
m1_m3 = lm(1,1); %'M1->M3'
pmv_m3 = lm(2,1); % 'PMV->M3'
s1_m3 = lm(3,1); % 'S1->M3'

% medial to lateral strAA={'M3-M1m','M3-PMVm', 'M3-S1'}; % medial to lateral flow 
if strcmp(beh,'T')
    ml= ml_t;
elseif strcmp(beh,'LS')
    ml=ml_l;
end
ml;
m3_m1 = ml(1,1); %'M3-M1m'
m3_pmv = ml(2,1); % 'M3-PMVm'
m3_s1 = ml(3,1); % 'M3-S1'

% lateral to lateral areas
% strAA={ 'M1m-PMVm', 'PMVm-M1m', 'M1m-S1','S1-M1m','PMVm-S1', 'S1-PMVm'}; % within lateral areas 
if strcmp(beh,'T')
    ll=ll_t;
elseif strcmp(beh,'LS')
    ll=ll_l;
end
ll;

m1_pmv = ll(1,1); % M1m-PMVm
pmv_m1 = ll(2,1); % PMVm-M1m
m1_s1 = ll(3,1); % M1m-S1
s1_m1 = ll(4,1); % S1-M1m
pmv_s1 = ll(5,1); % PMV_S1
s1_pmv = ll(6,1); % S1_PMV 

% General Variables 
factor=1;
factorGW=10;
NodeSize=40; %30
MarkSize=35; %30 
thisFS= 0.0001; %invisible; 

labels={'M3','PMV','S1','M1'};
s = [ 1 2 3 4 1 2];
t = [ 2 3 4 1 3 4];
w=repmat(0.001,length(s),1);

M_NodeColor=[0.9804    0.4745    0.0275];
L_NodeColor = [0    0.4471    0.7412];

%% Medial --> Lateral strLA={'M3-M1m','M3-PMVm', 'M3-S1'};;
w(1)= m3_pmv; % 1-2, M3-PMV  --> here 
w(2)= 0.0001; % 2-3, S1-PMV 
w(3)= 0.0001; % 3-4 S1-M1 
w(4)= m3_m1; % 1-4 M3-M1  --> here 
w (5)= m3_s1; % 1-3 % S1-M3 --> here 
w(6)= 0.0001; % 2-4 M1-PMV
strT='Granger Medial --> Lateral'

figure(10); clf;
G = graph(s, t,w);
plot(G,'EdgeLabel',G.Edges.Weight);
G1 = digraph(s,t,w);
LWidths = factorGW*[G1.Edges.Weight]';

figure(11)
h=plot(G1,'Layout','force','EdgeLabel',G1.Edges.Weight,'LineWidth',LWidths,'NodeLabel',labels);
colormap turbo          % select color palette 
h.EdgeCData=LWidths;    % define edge colors
h.MarkerSize=MarkSize; % define node size
h.NodeFontSize=NodeSize;
h.NodeColor=[0.5 0.5 0.5];
h.EdgeFontSize=thisFS;
h.ArrowSize=0.001;
highlight(h,[1],'NodeColor', M_NodeColor)

axis off
caxis manual
caxis([0 10]);
colorbar;
c=colorbar; 
temp=0:.1:1;
temp2=num2cell(temp);
c.TickLabels=temp2;
c.FontSize=12;
c.Location='southoutside'
c.Label.String =[strT, ' ', beh]
c.Label.FontSize =18;
c.Visible='on'



%% Lateral to Medial strLA={'M1->M3', 'PMV->M3', 'S1->M3'};
w(1)= pmv_m3; % 1-2, M3-PMV  --> here 
w(2)= 0.0001; % 2-3, S1-PMV 
w(3)= 0.0001; % 3-4 S1-M1 
w(4)= m1_m3; % 1-4 M3-M1  --> here 
w(5)= s1_m3; % 1-3 % S1-M3 --> here 
w(6)= 0.0001; % 2-4 M1-PMV
strT='Granger Lateral --> Medial';

figure(20); clf;
G = graph(s, t,w);
plot(G,'EdgeLabel',G.Edges.Weight);
G1 = digraph(s,t,w);
LWidths = factorGW*[G1.Edges.Weight]';

figure(21)
h=plot(G1,'Layout','force','EdgeLabel',G1.Edges.Weight,'LineWidth',LWidths,'NodeLabel',labels);
colormap turbo          % select color palette 
h.EdgeCData=LWidths;    % define edge colors
h.MarkerSize=MarkSize; % define node size
h.NodeFontSize=NodeSize;
h.NodeColor=[0.5 0.5 0.5];
h.EdgeFontSize=thisFS;
h.ArrowSize= 0.001;
highlight(h,[2 3 4],'NodeColor', L_NodeColor);

axis off
caxis manual
caxis([0 10]);
colorbar;
c=colorbar; 
temp=0:.1:1;
temp2=num2cell(temp);
c.TickLabels=temp2;
c.FontSize=12;
c.Location='southoutside'
c.Label.String =[strT, ' ', beh];
c.Label.FontSize =18;
c.Visible='on';


%% Lateral to Lateral Areas
% strAA={ 'M1m-PMVm', 'PMVm-M1m', 'M1m-S1','S1-M1m','PMVm-S1', 'S1-PMVm'}; % within lateral areas 
%%  M1 and S1 as inputs into PMV 
w(1)= 0.0001; % 1-2, M3-PMV  
w(2)= s1_pmv; % 2-3, S1-PMV 
w(3)= m1_s1; % 3-4 S1-M1 ---> here 
w(4)= 0.0001; % 1-4 M3-M1  
w(5)= 0.0001; % 1-3 % S1-M3 
w(6)= m1_pmv ; % 2-4 M1-PMV --> here 

strT='Granger Lateral --> Lateral'

figure(30); clf;
G = graph(s, t,w);
plot(G,'EdgeLabel',G.Edges.Weight);
G1 = digraph(s,t,w);
LWidths = factorGW*[G1.Edges.Weight]';

figure(31)
h=plot(G1,'Layout','force','EdgeLabel',G1.Edges.Weight,'LineWidth',LWidths,'NodeLabel',labels);
colormap turbo          % select color palette 
h.EdgeCData=LWidths;    % define edge colors
h.MarkerSize=MarkSize; % define node size
h.NodeFontSize=NodeSize;
h.NodeColor=[0.5 0.5 0.5];
h.EdgeFontSize=thisFS;
h.ArrowSize= 0.001;
highlight(h,[4],'NodeColor', L_NodeColor)
highlight(h,[3],'NodeColor', [0.0745    0.6235    1.0000])

axis off
caxis manual
caxis([0 10]);
colorbar;
c=colorbar; 
temp=0:.1:1;
temp2=num2cell(temp);
c.TickLabels=temp2;
c.FontSize=12;
c.Location='southoutside'
c.Label.String =[strT, ' ', beh];
c.Label.FontSize =18;
c.Visible='on'


%%  PMV as inputs into M1 and S1 
w(1)= 0.0001; % 1-2, M3-PMV  
w(2)= pmv_s1; % 2-3, S1-PMV 
w(3)= s1_m1; % 3-4 S1-M1 ---> here 
w(4)= 0.0001; % 1-4 M3-M1  
w(5)= 0.0001; % 1-3 % S1-M3 
w(6)= pmv_m1 ; % 2-4 M1-PMV --> here 

strT='Granger Lateral --> Lateral'

figure(40); clf;
G = graph(s, t,w);
plot(G,'EdgeLabel',G.Edges.Weight);
G1 = digraph(s,t,w);
LWidths = factorGW*[G1.Edges.Weight]';

figure(41)
h=plot(G1,'Layout','force','EdgeLabel',G1.Edges.Weight,'LineWidth',LWidths,'NodeLabel',labels);
colormap turbo          % select color palette 
h.EdgeCData=LWidths;    % define edge colors
h.MarkerSize=MarkSize; % define node size
h.NodeFontSize=NodeSize;
h.NodeColor=[0.5 0.5 0.5];
h.EdgeFontSize=thisFS;
h.ArrowSize= 0.001;
highlight(h,[2],'NodeColor', L_NodeColor);
highlight(h,[3],'NodeColor', [0.0745    0.6235    1.0000]);

axis off
caxis manual
caxis([0 10]);
colorbar;
c=colorbar; 
temp=0:.1:1;
temp2=num2cell(temp);
c.TickLabels=temp2;
c.FontSize=12;
c.Location='southoutside'
c.Label.String =[strT, ' ', beh];
c.Label.FontSize =18;
c.Visible='on';



% store the figures to join with ProtoFig3_PPC_Granger
if PrintOn==1
    figure(11); gcf;
    savefig(['Fig3_Node_ML_', beh])
    print(['Fig3_Nodes_ML_', beh],'-dsvg')

    figure(21); gcf;
    savefig(['Fig3_Node_LM_', beh])
    print(['Fig3_Node_LM_', beh],'-dsvg')
 
    figure(31); gcf;
    savefig(['Fig3_Node_LL_M1_', beh])
    print(['Fig3_Node_LL_M1_', beh],'-dsvg')

    figure(41); gcf;
    savefig(['Fig3_Node_LL_PMV_', beh])
    print(['Fig3_Node_LL_PMV_', beh],'-dsvg')
end


