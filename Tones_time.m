%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Modelling Inner Hair Cell Transduction
%                                 Main
%                           with Tone Signals
%
%  References: 
%  Relaño-Iborra 2019: A speech-based computational auditory signal 
%                      processing and perception model
%  Sumner et al. 2002: A revised model of the IHC and auditory-nerve complex
%
%  Notes
%   - tone signal inputs at specific level
%   - plots receptor potentials vs time for different frequencies
%   - plots AC/DC ratio of receptor potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
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
fc = erbspace(flow, fhigh, 60); % choose centre frequencies
%fc = [5000 7000];
drnlchannels = size(fc,2);

sbj='NH'; % DRNL profiles:'NH'/'OHCloss_5dB'/'OHCloss_10dB'/'OHCloss_20dB'/'OHCloss_30dB'/'HIx'

% Tone signal parameters
fs = 44100;             % Sampling frequency
t_tot = 60e-3;          % Signal duration
Nsamp = t_tot*fs;       % Samples
ramp_dur = 5e-3;        % ramp duration in seconds
N_zeros = 20e-3*fs;     % zero-pad signal

t = (0:Nsamp-1)/fs;         % time vector
f = (0:Nsamp-1)*fs/Nsamp;   % freq axis

f_probe = [100 300 500 700 900 1000 2000 3000 4000 5000];   % input signal frequencies
nprobs= length(f_probe);

level = 80;   % tone level
amp = 20e-6.*10.^(level./20); % Amplitude in Pa;

% Ramp up/downn
window = hanning(2*floor(fs*ramp_dur))'; 
w1 = window(1:ceil((length(window))/2)); 
w2 = window(ceil((length(window))/2)+1:end); 
w = [w1 ones(1,Nsamp-length(w1)-length(w2)) w2]; 
 

%% ------------- RUN MODEL ------------------ %%
output = zeros(nprobs,Nsamp+N_zeros,drnlchannels);
output_AC = zeros(nprobs,drnlchannels);
output_DC = zeros(nprobs,drnlchannels);

for r = 1:nprobs % Loop over tone frequencies
        
      x = amp*sin(2*pi*f_probe(r)*t); % Tone
      x = x.*w;
      x = [x zeros(1,N_zeros)];

      N_org= length(x);
    
      [~,outsig_nl,outsig_drnl,fc] = IHC_output(x, fs, fc,'drnl', N_org, 'NH');

      output(r,:,:) = outsig_nl;
      
      output_AC(r,:) = abs(max(squeeze(output(r,1000:2000,:)),[],1)-min(squeeze(output(r,1000:2000,:)),[],1));
      output_DC(r,:) = abs(min(squeeze(output(r,1000:2000,:)),[],1)+0.5*output_AC(r,:)-squeeze(output(r,1,:)).');
end

%% AC/DC RATIO

for r = 1:nprobs
    [BF_dif, BF_indx] = min(abs(fc-f_probe(r)));
    %BF_indx = 1;
    ACs(r) = output_AC(r,BF_indx);
    DCs(r) = output_DC(r,BF_indx);
end
ratio = ACs./DCs;

load('data_palmer.mat')

% PLOT 
figure
loglog(f_probe, ratio,'linewidth', linesize);
hold on
scatter(Data_01(:,1)*1e3,Data_01(:,2),80,'filled')
grid on
%xlim([0.1e3 5e3])
ylim([0.1 10])
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

%% Plot Freqs. -  Best frequency channel

t_axis = (0:Nsamp+N_zeros-1)/fs;         % time vector

i=1;
figure
for r = 1:nprobs
    subplot(5,2,r)
    [BF_dif, BF_indx] = min(abs(fc-f_probe(r)));
    plot(t_axis/1e-3, squeeze(output(r, :, BF_indx)/1e-3),'linewidth', linesize);
    hold on
    str2= sprintf('f = %.0f Hz',f_probe(r));
    legend(str2)
    grid on
    xlim([0 80])
    ylim([-70 -20])
    xlabel('t in ms')
    ylabel('Potential in mV')
    set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);
    i = i+1;
end
