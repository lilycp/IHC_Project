%%%%%%%%%
%
% Retrieving fitting data from results structure and calculating the values
% of the free parametters of the mapping curve.
%
% Helia Relaño Iborra, 2020
%%%%%%

%clear
%%  ---------- Read values from results structure ---------- %%

load('res_CLUEfit_sCASP_65_me2')


for j=1:length(res_fit_sCASP) % For all the sentences

dfinal_clue(j, :)=res_fit_sCASP(j).dfinal;
       
end

%% mean outputs
d_mean_clue= squeeze(nanmean(dfinal_clue, 1));
dev_d= squeeze(nanstd(dfinal_clue, 1));

%% ----------- Fit logistic function ------------- %%

%Human Data
Pcorrect_human_CLUE = [0 8 35 71 90 100 ]; % Human data from Nielsen & Dau (2009)

% Define the logistic function:
fun = @(a,xdata) 100./(1 + exp(a(1)*xdata + a(2))); %logistic function

% Optimization parameters: 
pguess = [0 0];  %starting guess
xdataC = d_mean_clue; % X-values
ydataC = Pcorrect_human_CLUE; % Y-values

[fit_param,R,J,CovB,MSE] = nlinfit(xdataC,ydataC,fun,pguess); % non linear fitting

% Goodness of fit
ysim= fun(fit_param, xdataC);
rsq2 = 1 - sum(R.^2) / sum((ydataC- mean(ydataC)).^2);

%% ---------  Plot ----------- %%
 
colors = [0.67 .0 0.77; 0 0.66 0.47; 0 0 0 ];
linesize       = 5.5;      % Set thickness of plot lines.
axislabelsize  = 35;       % Set size of axis labels (text).
axisnumbersize = 30;       % Set size of axis numbers.


% "smooth" fitted function:
x = linspace(0, 1, 500); 
fit_functC= fun(fit_param, x);

figure
h(1) = scatter(d_mean_clue(:), Pcorrect_human_CLUE, 250, colors(2, :), 'o', 'filled');

hold on

h(2) = plot(x, fit_functC, '-', 'color', colors(3, :), 'linewidth', linesize);


text(0.5,92,{[' R^2 = ',num2str(rsq2,2)]},'FontSize', axisnumbersize,  'fontweight', 'bold','FontName', 'Colibri');

xlabel('Correlation Index','FontSize',16, 'FontName', 'Colibri')
ylabel('% correct','FontSize',16, 'FontName', 'Colibri')

le=legend(h, 'Nielsen & Dau (2009)', 'f_{CLUE}',...
          'Location', 'SouthEast','FontSize', axisnumbersize,  'fontweight', 'bold');
      
set(le, 'box', 'off', 'FontSize', axisnumbersize,  'fontweight', 'bold', 'location', 'southeast')

grid on
set(gca, 'Fontsize', 13)
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 12 7])
print('-dpng',fullfile('fitCLUE_curve.png' ))

