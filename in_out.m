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
%   - tone signal inputs at different frequencies and levels
%   - input/output functions 
%   - Potential vs Level for different frequencies
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
fc = erbspace(flow, fhigh, 60); % choose centre frequencies
drnlchannels = size(fc,2);

sbj='NH'; % DRNL profiles:'NH'/'OHCloss_5dB'/'OHCloss_10dB'/'OHCloss_20dB'/'OHCloss_30dB'/'HIx'

% Tone signal parameters
fs = 44100;                 % Sampling frequency
t_tot = 1;                  % Signal duration
Nsamp = t_tot*fs;           % Samples
t = (0:Nsamp-1)/fs;         % time axis
f = (0:Nsamp-1)*fs/Nsamp;   % freq axis
N_zeros = 20e-3*fs;         % zero-pad signal

% Ramp up/downn
ramp_dur = 5e-3;        % ramp duration in seconds
window = hanning(2*floor(fs*ramp_dur))'; 
w1 = window(1:ceil((length(window))/2)); 
w2 = window(ceil((length(window))/2)+1:end); 
w = [w1 ones(1,Nsamp-length(w1)-length(w2)) w2]; 

f_probe = 100; % Tone frequency

levels = 10:10:90;  % tone levels
nlevels=length(levels);

amp = 20e-6.*10.^(levels./20); % Amplitude in Pa;

%% ------------- RUN MODEL ------------------ %%

for j = 1:nlevels  % Loop over levels
    
  x = amp(j).*sin(2*pi*f_probe*t); % Tone
  x = x.*w;
  x = [x zeros(1,N_zeros)];

  N_org= length(x);

  [~,outsig_nl,~,fc] = IHC_output(x, fs, fc,'drnl', N_org, 'NH');

   output(j,:,:) = outsig_nl;
      
   peak_pos(j,:) = max(outsig_nl);
   peak_neg(j,:) = min(outsig_nl);

end


[~, idx] = min(abs(fc-f_probe));


%% PLOT 
amp_axis = [-flip(amp) amp];
peak_axis = [flip(peak_neg(:,idx)); peak_pos(:,idx)];

figure
%plot(amp, peak_pos(:,idx),'linewidth', linesize);
plot(amp_axis, peak_axis,'-o','linewidth', linesize);
grid on
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

