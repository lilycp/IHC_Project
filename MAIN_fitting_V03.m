%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Modelling Inner Hair Cell Transduction
%                                 Main V03
%              Fitting function of cilia to conductance stage
%
%  References: 
%  Relaño-Iborra 2019: A speech-based computational auditory signal 
%                      processing and perception model
%  Sumner et al. 2002: A revised model of the IHC and auditory-nerve complex
%
%  Notes
%   - fitting of second stage (Ster Disp --> Conductance)
%   - Data from Palmer & Russell used to fit
%   - nlinfit 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

addpath(genpath('/Users/liliana/Documents/PhD/MATLAB/IHC_project'))

linesize       = 1.5;      % Set thickness of plot lines.
axislabelsize  = 10;       % Set size of axis labels (text).
axisnumbersize = 10;       % Set size of axis numbers.

%% PARAMETERS

fs = 44100;             % Sampling frequency
Ts = 1/fs;

Et = 100E-3;
Ek = -70.45E-3;
Omega = 40E-3;    %correction factor
Ekp = Ek+Et*Omega;
Gk = 18E-9;
Cm = 6E-12;

G0 = 1.974E-9;
Gmax = 8E-9;
s0 = 85E-9;
u0 = 7E-9;
s1 = 5E-9;
u1 = 7E-9;
Ga = G0-Gmax*1/((1+exp(u0/s0))*(1+exp(u1/s1)));

%% TEST
f_probe = [100 300 500 700 900 1000 2000 3000 4000 5000];   % input signal frequencies
param = [2.13e-3 0.1];

ratio = IHC_ACDC_V02(f_probe,param);

figure
loglog(f_probe, ratio,'linewidth', linesize);
grid on
%xlim([0.1e3 5e3])
ylim([0.1 10])
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

%% FITTING

load('data_palmer.mat')

handle = @IHC_ACDC_V02;
fun  = @(a,xdata) IHC_ACDC_V02(xdata,a);

%guess = [u0 s0 u1 s1];  % starting guess
guess = [9.45E-9 52.7E-9 63.1E-9 29.4E-9 12.7E-9 0.33E-9]; % Lopez

%guess = [9.45E-9 52.7E-9 63.1E-9 29.4E-9 12.7E-9 0.33E-9 2.13E-3 0.1];
%guess = [2.13E-3 0.1];

xdata = Data_01(:,1)*1e3;      % X values
ydata = Data_01(:,2);          % Y values

%opts = statset('nlinfit');
%opts.FunValCheck = 'off';

[param,~,~,~,~]   = nlinfit(xdata,ydata,fun,guess);  

x = linspace(100,4000,50);
fit_funct = fun(param, x);

figure
loglog(x/1e3, fit_funct, '-', 'LineWidth', linesize);
hold on
scatter(xdata/1e3,ydata,80, 'o', 'filled');
grid on
ylim([0.1 10])
xlabel('Frequency in kHz','FontSize',20) 
ylabel('AC/DC component','FontSize',20)
title('80dB')
yticks([0.1 0.3 1 3 ])
yticklabels({'0.1','0.3','1','3','10'})
xticks([0.1 1 3 10 ])
xticklabels({'0.1','1','3','10'})
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

%% Check Boltzmann function

dis_vs_cond = [-95.42483660130719, 0.016393442622952392
-76.79738562091504, 0.016393442622952392
-65.68627450980392, 1.7763568394002505e-15
-45.42483660130718, 0.016393442622952392
-36.92810457516339, 1.7763568394002505e-15
-17.320261437908485, 1.7763568394002505e-15
-13.071895424836583, 0.049180327868853624
-0.32679738562090677, 0.06557377049180424
7.8431372549019756, 0.44262295081967196
9.150326797385631, 1.196721311475411
18.62745098039217, 1.4918032786885247
31.37254901960786, 1.9836065573770494
36.9281045751634, 2.9836065573770494
36.9281045751634, 3.967213114754098
52.28758169934643, 4.049180327868852
58.49673202614382, 4.9016393442622945
64.37908496732028, 5.1803278688524586
77.77777777777779, 5.508196721311475
86.60130718954251, 6.0655737704918025
99.67320261437912, 6.475409836065573
106.86274509803924, 6.491803278688524
114.7058823529412, 6.868852459016393
124.18300653594771, 7.278688524590163
135.62091503267982, 7.5409836065573765];

dis_vs_cond(:,1) = dis_vs_cond(:,1)*1E-9;

x = linspace(min(dis_vs_cond(:,1))-1E-7,max(dis_vs_cond(:,1))+1E-7, 500);
a = param;

G = a(1)*1./((1+exp(-(x-a(2))/a(3))).*(1+exp(-(x-a(4))/a(5))))+ a(6) ;    
%G = Gmax*1./((1+exp(-(x-a(1))/a(2))).*(1+exp(-(x-a(3))/a(4)))); % + a(5);    

figure
scatter(dis_vs_cond(:,1)/1E-9, dis_vs_cond(:,2),80, 'o', 'filled');
grid on
hold on
plot(x/1E-9, G/1E-9, '-', 'LineWidth', linesize);

