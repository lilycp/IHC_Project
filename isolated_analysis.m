%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Modelling Inner Hair Cell Transduction
%                                 Main
%                          with Tone Signals
%
%  References: 
%  Relaño-Iborra 2019: A speech-based computational auditory signal 
%                      processing and perception model
%  Sumner et al. 2002: A revised model of the IHC and auditory-nerve complex
%
%  Notes
%   - tone signal inputs one freq & one level
%   - isolated analysis of each stage
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all, 
clear
clc

addpath(genpath('/Users/liliana/Documents/PhD/MATLAB/IHC_project'))

linesize       = 1.5;      % Set thickness of plot lines.
axislabelsize  = 10;       % Set size of axis labels (text).
axisnumbersize = 10;       % Set size of axis numbers.

%% PARAMETERS

% DRNL frequencies
flow=100;
fhigh= 8e3;
%fc = erbspace(flow, fhigh, 60); % choose centre frequencies
fc = 4e3;
drnlchannels = size(fc,2);

sbj='NH'; % DRNL profiles:'NH'/'OHCloss_5dB'/'OHCloss_10dB'/'OHCloss_20dB'/'OHCloss_30dB'/'HIx'

% Tone signal parameters
fs = 44100;             % Sampling frequency
Ts = 1/fs;
t_tot = 60e-3;          % Signal duration
Nsamp = t_tot*fs;       % Samples
ramp_dur = 5e-3;        % ramp duration in seconds
N_zeros = 20e-3*fs;     % zero-pad signal

t = (0:Nsamp-1)/fs;         % time axis
f = (0:Nsamp-1)*fs/Nsamp;   % freq axis

f_probe = 4000;

levels = 10;

amp = 20e-6.*10.^(levels./20); % Amplitude in Pa;

% Ramp up/downn
window = hanning(2*floor(fs*ramp_dur))'; 
w1 = window(1:ceil((length(window))/2)); 
w2 = window(ceil((length(window))/2)+1:end); 
w = [w1 ones(1,Nsamp-length(w1)-length(w2)) w2]; 

%% ------------- RUN MODEL ------------------ %%

x = amp.*sin(2*pi*f_probe*t); % Tone
x = [x.*w zeros(1,N_zeros)];
N_org= length(x);

b_hp = HeadphoneFilter(fs);     % calc headphone filtercoeffs
b_me = middleearfilter_v2(fs);  % calc middle ear filtercoeffs

inoutsig = filter(b_hp,1,x); % Outer-ear filterring
inoutsig = filter(b_me,1,inoutsig); % middle-ear-ear filterring

% DRNL filter
outsig = drnl_HI(inoutsig,fc,fs, sbj);

% IHC TRANSDUCTION 
tau = 2.13E-3;   %low-pass filter applied to drnl output
C = 0.1;

b = [0 tau*C];
a_filt = [tau 1];

[b_z,a_z] = bilinear(b,a_filt,pi);
ster_disp = filter(b_z,a_z,outsig,[],1); % stereocilia displacement        

% Conductance
Et = 100E-3;
Ek = -70.45E-3;
Omega = 40E-3;    % correction factor
Ekp = Ek+Et*Omega;
Gk = 18E-9;
Cm = 6E-12;
Gmax = 8E-9;
%a = [8E-9 7E-9 85E-9 7E-9 5E-9 0.215E-9]; %SOmner
a = [9.45E-9 52.7E-9 63.1E-9 29.4E-9 12.7E-9 0.33E-9]; % Lopez

G = a(1)*1./((1+exp(-(ster_disp-a(2))/a(3))).*(1+exp(-(ster_disp-a(4))/a(5)))) + a(6) ; 
G0 = a(1)*1./((1+exp(a(2)/a(3))).*(1+exp(a(4)/a(5)))) + a(6);

% Intracellular Potential
V_rest = (Gk*Ekp+G0*Et)./(G0+Gk);
V_now = V_rest;

V = zeros(1,Nsamp+N_zeros);
    for ii=1:Nsamp+N_zeros
            V_now = V_now + (-G(:,ii).*(V_now-Et)-...
                Gk*(V_now-Ekp))*Ts/Cm;
            V(:,ii)=V_now;
    end

% Reference - envelope extraction
signal = envextract(outsig,fs);

%% PLOT 

t_axis = (0:Nsamp+N_zeros-1)/fs;  % time vector

figure 
subplot(4,1,1)
plot(t_axis/1e-3, outsig/1e-3,'linewidth', linesize);
grid on
xlim([0 80])
xlabel('t in ms')
ylabel('BM velocity [mm/s]')
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

subplot(4,1,2)
plot(t_axis/1e-3, ster_disp/1e-9,'linewidth', linesize);
grid on
xlim([0 80])
xlabel('t in ms')
ylabel('cilia disp [nm]')
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

subplot(4,1,3)
plot(t_axis/1e-3, G/1e-9,'linewidth', linesize);
grid on
xlim([0 80])
xlabel('t in ms')
ylabel('Conductance [nS]')
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

subplot(4,1,4)
plot(t_axis/1e-3, V/1e-3,'linewidth', linesize);
grid on
xlim([0 80])
%ylim([-60 -30])
xlabel('t in ms')
ylabel('Potential [mV]')
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

figure 
plot(t_axis/1e-3, signal/1e-3,'linewidth', linesize);
grid on
xlim([0 80])
ylim([0 0.7])
xlabel('t in ms')
ylabel('Output [mMU]')
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

