%%%%%%%%%%%%%%%%%%%%
%  Analysis of the fluct-CLUE model outputs, transformation into %correct
%  and calculation of SRTs. Includes plotting of the resulting predictions
%  in comparison to human data from Jorgensen et al. (2013)
%
%  Helia Relaño Iborra, 2019
%
%%%%%%%%%%%%%%%%%%%%%

clear
%% ----------  Experiment parameters and human data: ----------- %%

conditions = {'SSN'  'SAM' 'ISTS'};
SNRs = -27:3:3;

% SRTs for CLUE Jørgensen et. al. 2013:
SRTs_Jetal2013 =  [-3.5 -9.5   -18.4]; % [SSN SAM ISTS]
SRTs_Jetal2013_std = [ 0.45 1.9   0.74 ];


%%  ---------- Read values from results structure ---------- %%

load('res_fluctCLUE_sCASP')

for j=1:length(res_fluct_sCASP) % For all the sentences
    
           dfinal_clue(j, :, :)=res_fluct_sCASP(j).dfinal;

end

%% ---------  Process data -------------------- %%

% Correlation index
d_mean = squeeze(nanmean(dfinal_clue, 1)); 
d_std = squeeze(nanstd(dfinal_clue, 1)); 

b = [ -64.9860   12.5050]; % Fitted parameters

% Transformation to percentage correct:
Pcorrect_casp=  100./(1 + exp(b(1)*d_mean' + b(2)));

% SRT estimation:
[dSRT_casp_rect SRTs_casp] = dsrts_from_pc_mean(Pcorrect_casp', SNRs, 1:3); 


%% ---------------  Accuracy metrics -------------- %%

MAE_casp = mean(abs(SRTs_casp -     SRTs_Jetal2013));
r_corr_casp = corr(SRTs_casp',    SRTs_Jetal2013', 'type', 'Pearson');
RMSE_casp= sqrt((sum(SRTs_casp- SRTs_Jetal2013).^2)/length(SRTs_Jetal2013));
 
%% ---------------- Plots ---------------- %%   
    
xmin = 0.5;
xmax = 3.5;
ymin = -30;
ymax = 5;
ytickmax = 5;
ytickmin = -30;
fnts = 12;
fig=figure;

sh = 500;
set(fig,'Position',[2*sh, 0.15*sh, 1.2*sh, 1*sh]);
 
    mark_sty = {'s','d','>','o'};
    mark_size = 10;
    mark_col(1,:) = [1 1 1]*0.7;
    mark_col(2,:) = [1 1 1]*0.4;
    mark_col(3,:) = [1 1 1]*0.6;
    mark_col(4,:) = [1 1 1]*0.01;

col=[0.67 .0 0.77; 0 0.66 0.47; 0 0 0 ];

    x = 1:3;
    offset = 0.04;
    
    h = plot(0,0,'color',[1 1 1]);
    
errorbar(x,  SRTs_Jetal2013, SRTs_Jetal2013_std,...
            'linestyle',    'none',...
            'linewidth',    1,...
            'color',           mark_col(4,:),...
            'marker',           char(mark_sty(1)),...
            'markerfacecolor',  [1 1 1],...
            'markersize',       mark_size); hold on          

 idx_stim = 4;    
    plot(x+offset,SRTs_casp,...
            'linestyle',    'none',...
            'linewidth',    1,...
            'color',          col(2,:),...
            'marker',            char(mark_sty(idx_stim)),...
            'markerfacecolor', col(2,:),...
            'markersize',       mark_size);hold on
                
    ylabel('SRT [dB]','FontSize',fnts, 'FontName', 'Colibri');
    ylim([ymin ymax])
    xlim([xmin xmax])

  le = legend('Human data','sCASP');
       set(le,'box','off','fontsize',12,'location','southwest');   
       
 text(2.5,0,{[' \rho_{sCASP} = ',num2str(r_corr_casp,2)],['MAE_{sCASP} = ',num2str(MAE_casp,2), ' dB']},'fontsize',12,'FontName', 'Colibri');

    set(gca,'xTick',1:4 ,'xTickLabel',conditions,'ytick',ytickmin:5:ytickmax,'FontSize',fnts,'yticklabel',ytickmin:5:ytickmax);
    set(gca, 'FontName', 'Colibri');
    set(gcf, 'Color', 'w');

set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 12 7])
print('-dpng',fullfile('figures', 'fluctCLUE_results.png' ))

