%%%%%
%
% Data analsis from Christiansen and Dau, 2012
%
%%%%%%

Subject_IDs = {'HI1','HI2','HI3','HI4','HI5','HI6','HI7','HI8','HI9','HI10','HI11','HI12','HI13'};

% Audiograms:

audG_L=[10 15 10 15 15 15 25 60 65 65 65 ;...
    20 30 35 35 35 45 50 60 65 55 55;...
    15 15 15 20 20 15 20 50 55 60 70 ;...
    15 25 40 50 45 50 45 50 55 50 65 ;...
    20 25 35 30 30 30 25 30 45 60 70 ;...
    20 35 55 65 65 60 55 55 60 65 70;...
    15 20 35 45 50 50 55 60 60 65 75;... 
    35 35 35 40 35 30 35 55 55 60 65;... 
    25 25 35 45 50 55 65 80 105 100 100;... % Considered outlier
    25 25 25 25 25 25 30 40 45 55 55;...
    20 25 30 40 40 55 50 50 60 85 90;... 
    25 40 55 50 55 55 55 80 80 70 75;... % Considered outlier
    25 35 35 45 45 45 50 70 80 80 80];

audG_R=[5 10 10 15 10 15 30 55 55 60 75;...
    20 35 35 35 40 45 45 50 55 50 50;...
    15 15 10 10 15 10 15 35 55 50 60;...
    20 25 40 60 60 55 60 55 50 40 55;...
    20 30 35 35 30 30 25 40 50 60 70;...
    25 45 55 65 55 55 50 55 60 55 70;...
    20 25 35 50 55 55 60 65 70 70 70;...
    25 30 35 40 35 30 30 50 50 55 60;...
    25 20 30 45 45 85 95 115 110 105 105;... % Considered outlier
    25 25 30 35 30 30 35 55 55 55 60;...
    35 40 40 50 50 65 65 60 75 95 100;...
    35 60 60 60 55 50 65 70 70 70 75;... % Considered outlier
    25 30 35 35 35 40 55 80 75 75 75];


freq=[125 250 500 750 1000 1500 2000 3000 4000 6000 8000];

%% Averaging 

% Remove outliers 
% audG_L([9, 12], :)= [];
% 
% audG_R([9, 12], :)= [];


% BE- selection

BE = min(audG_L, audG_R);

% Average audiogram

BE_ave = mean(BE);

both_ears = [audG_L; audG_R];

%% Ploting

figure,

plot(1:11, -audG_L', 'color', [0.85 0.85 0.85]);
hold on
plot(1:11, -audG_R', 'color', [0.85 0.85 0.85]);
hold on
plot(1:11, -mean(both_ears), 'color', 'k');

 set(gca,'xTick',1:11 ,'xTickLabel',freq,...
     'ytick', -120:10:0, 'yticklabel', 120:-10:0);
axis([0.2 11.5 -125 5])
grid on
set(gca,'GridLineStyle','--')
ref=refline(0, 0);
set(ref, 'color', 'k');
xlabel('Frequency (Hz)')
ylabel('Thresholds (dB HL)')

set(gca, 'fontname', 'times', 'fontsize', 12)

set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 12])
% print('-dpng','-r300',fullfile('figures', 'audiograms_christiansen_w_mean.png' ))
% print('-depsc','-r900',fullfile('figures', 'audiograms_christiansen_w_mean.eps' ), '-painters')

print('-dsvg','-r900',fullfile('audiograms_christiansen_w_mean.svg' ), '-painters')