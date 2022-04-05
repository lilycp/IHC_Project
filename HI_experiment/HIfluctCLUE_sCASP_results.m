%%%%%%%%%%%%%%%%%%%%
%  Analysis of the fluct-CLUE model outputs, transformation into %correct
%  and calculation of SRTs. Includes plotting of the resulting predictions
%  in comparison to human data from Jorgensen et al. (2013)
%
%%%%%%%%%%%%%%%%%%%%%

clear
%% ----------  Experiment parameters and human data: ----------- %%

conditions = {'SSN'  'SAM' 'ISTS'};
SNRs = -25:5:15;

load('sbj_IDs_christiansen2012_2ears.mat')

names = {'NH' ;'HI1'; 'HI2'; 'HI3'; 'HI4'; 'HI5'; 'HI6'; 'HI7'; 'HI8';...
    'HI9';'HI10';'HI11';'HI12';'HI13'};

%%  ---------- Read values from results structure ---------- %%


load('res_fluct_CLUE_sCASP_HI.mat')

for j=1:length(res_fluct_sCASP) % For all the sentences
    
    dfinal_clue(j, :, :, :)=res_fluct_sCASP(j).dfinal;
    
end

%% ---------  Process data -------------------- %%

% Correlation index
d_mean = squeeze(nanmean(dfinal_clue, 1));


b = [-63.7644   12.3352]; % Fitted parameters

for s =1:size(d_mean, 1) %% Loop over subjects:
    
    % Transformation to percentage correct:
    Pcorrect_casp(s, :, :)=  100./(1 + exp(b(1)*squeeze(d_mean(s, : ,:)) + b(2)));
    
    % SRT estimation:
    [dSRT_casp_rect(s, :) SRTs_casp(s, :)] = dsrts_from_pc_mean(squeeze(Pcorrect_casp(s, :, :)), SNRs, 1:3);
    
end


%% NH model results:

SRTs_NH_pred =SRTs_casp(1, :);
% SRTs_NH_pred = load('srts_casp_fluct_NH.mat');  % Saved NH results 


%% Choice of the best ear for HI model data:

SRTs_L =SRTs_casp(2:14, :);
SRTs_R =SRTs_casp(15:end, :);

SRTs_BE = min(SRTs_L,SRTs_R  );

SRTs_casp_mean=nanmean(SRTs_BE, 1);
SRTs_casp_std=nanstd(SRTs_BE, 1);

%% Reference human data:

SRTs_HI=[ 1.200 -3.400	-7.000;...
    -1.100 	-3.100	-11.400;...
    0.000	-4.600	-10.600;...
    -0.100	-1.000	-5.300;...
    -0.100	-5.100	-8.700;...
    3.200	5.600	3.500;...
    0.400	0.800	-3.500;...
    -0.800	-4.300	-10.000;...
    8.200	9.800	9.000;... % Considered outlier
    -1.500	-5.500	-10.500;...
    -1.900	-2.800	-9.900;...
    2.200	3.300	1.700;... % Considered outlier
    2.700	2.400	-5.100];

SRs_NH=[-3.1353	-9	-18.3];

SRTs_HI_ave=mean(SRTs_HI, 1);
SRTs_HI_std=std(SRTs_HI, 1);


%% --------- Boxplots --------------- %%

cases ={'Human'; 'sCASP'};

cols = [0.7 0.7 0.7; 1 1 1];
figure
ha = tight_subplot(1,3,[.05 .01],[.1 .1],[.1 .1]);

axes(ha(1)) % Human data
H_all = [SRTs_HI(:, 1); SRTs_BE(:, 1)];
g = [zeros(length(SRTs_HI(:, 1)), 1); ones(length(SRTs_BE(:, 1)), 1)];
b1 = boxplot(H_all,g, 'Colors', [0 0 0]);
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects


h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),cols(j, :),'FaceAlpha',.5);
end



% Add NH data:

hold on

idx_stim = 1;
l(1) =   plot(2,SRTs_NH_pred(1) ,...
    'linestyle',    'none',...
    'linewidth',    1,...
    'color',        [0 0 0],...
    'marker',            'd',...
    'markerfacecolor', [0.7 0.7 0.7],...
    'markersize',       mark_size);



l(2) = plot(1,  SRTs_NH(1), ...
    'linestyle',    'none',...
    'linewidth',    1,...
    'color',           [ 0 0 0],...
    'marker',           'd',...
    'markerfacecolor',  [1 1 1],...
    'markersize',       mark_size); hold on





set(b1,{'linew'},{1})
ylim([-20 11])

le = legend('HI modelled', 'HI measured', 'NH modelled', 'NH measured');
set(le, 'box', 'off', 'location', 'southwest')
set(le, 'Fontsize', 12, 'Fontname', 'Times')

text(1.3,9.5,'SSN',  'Fontsize', 15, 'Fontname', 'Times', 'fontweight', 'bold')
set(gca, 'xtick', 1:2, 'xticklabel',cases)
% ylabel('SRT(dB SPL)')
set(gca, 'Fontsize', 15, 'Fontname', 'Times')


axes(ha(2)) % scasp

H_all = [SRTs_HI(:, 2); SRTs_BE(:, 2)];
% H_all = [scaspNH; scasp.SRTs_casp_v2'];
g = [zeros(length(SRTs_HI(:, 2)), 1); ones(length(SRTs_BE(:, 2)), 1)];
% g = [zeros(length(scaspNH), 1); ones(length(scasp.SRTs_casp_v2), 1)];
b3 = boxplot(H_all,g, 'Colors', [0 0 0]);
set(b3,{'linew'},{1})

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),cols(j, :),'FaceAlpha',.5);
end

hold on
plot(1,  SRTs_NH(2), ...
    'linestyle',    'none',...
    'linewidth',    1,...
    'color',           [0 0 0],...
    'marker',           'd',...
    'markerfacecolor',  [1 1 1],...
    'markersize',       mark_size); hold on


idx_stim = 1;
plot(2,SRTs_NH_pred(2) ,...
    'linestyle',    'none',...
    'linewidth',    1,...
    'color',        [0 0 0],...
    'marker',            'd',...
    'markerfacecolor', [0.7 0.7 0.7],...
    'markersize',       mark_size);


ylim([-20 11])
set(gca, 'xtick', 1:2, 'xticklabel',cases)
set(gca, 'ytick',[], 'yticklabel',[])
text(1.3,9.5,'SAM',  'Fontsize', 15, 'Fontname', 'Times', 'fontweight', 'bold')

set(gca, 'Fontsize', 15, 'Fontname', 'Times')

axes(ha(3)) %scorr
%
H_all = [SRTs_HI(:, 3); SRTs_BE(:, 3)];

g = [zeros(length(SRTs_HI(:, 3)), 1); ones(length(SRTs_BE(:, 3)), 1)];
b2 = boxplot(H_all,g, 'Colors', [0 0 0]);
set(b2,{'linew'},{1})

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),cols(j, :),'FaceAlpha',.5);
end

hold on
plot(1,  SRTs_NH(3), ...
    'linestyle',    'none',...
    'linewidth',    1,...
    'color',           [0 0 0],...
    'marker',           'd',...
    'markerfacecolor',  [1 1 1],...
    'markersize',       mark_size); hold on


idx_stim = 1;
plot(2,SRTs_NH_pred(3) ,...
    'linestyle',    'none',...
    'linewidth',    1,...
    'color',        [0 0 0],...
    'marker',            'd',...
    'markerfacecolor', [0.7 0.7 0.7],...
    'markersize',       mark_size);

ylim([-20 11])
set(gca, 'xtick', 1:2, 'xticklabel',cases)
set(gca, 'ytick',[], 'yticklabel',[])
text(1.3,9.5,'ISTS',  'Fontsize', 15, 'Fontname', 'Times', 'fontweight', 'bold')
set(gca, 'Fontsize', 15, 'Fontname', 'Times')



set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 18 10])
print('-dpng', '-r300',fullfile('Christiansen_group_boxplots.png'))

