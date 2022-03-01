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
%   - revised & summarised from Special Project
%   - tone signal inputs
%   - input/output functions 
%   - Potential vs Level for different frequencies
%   - Aim to reproduce figures from Patuzzi/Sellick 1983 (Fig. 1B & 6)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all, 
clear
clc

addpath(genpath('/Users/lpau/Desktop/PhD/MATLAB/IHC_project'))

%% PARAMETERS

% DRNL frequencies
flow=100;
fhigh= 8e3;
%fc = erbspace(flow, fhigh, 60); % choose centre frequencies
fc = [4e3 6e3];
drnlchannels = size(fc,2);

sbj='NH'; % DRNL profiles:'NH'/'OHCloss_5dB'/'OHCloss_10dB'/'OHCloss_20dB'/'OHCloss_30dB'/'HIx'

% Tone signal parameters
fs = 44100;             % Sampling frequency
t_tot = 60e-3;          % Signal duration
Nsamp = t_tot*fs;       % Samples
ramp_dur = 5e-3;        % ramp duration in seconds
N_zeros = 20e-3*fs;     % zero-pad signal

t = (0:Nsamp-1)/fs;         % time axis
f = (0:Nsamp-1)*fs/Nsamp;   % freq axis

%f_probe = [1 1.04 1.09 1.11]*fc(1);
f_probe = [1 0.89 0.72 0.39]*fc(1);
nprobs= length(f_probe);

levels = 10:10:100;  % tone levels
%levels = 65;
nlevels=length(levels);

amp = 20e-6.*10.^(levels./20); % Amplitude in Pa;

% Ramp up/downn
window = hanning(2*floor(fs*ramp_dur))'; 
w1 = window(1:ceil((length(window))/2)); 
w2 = window(ceil((length(window))/2)+1:end); 
w = [w1 ones(1,Nsamp-length(w1)-length(w2)) w2]; 

%% ------------- RUN MODEL ------------------ %%

% preallocation 
output = zeros(nlevels,nprobs,Nsamp+N_zeros);
output_AC = zeros(nlevels,nprobs,drnlchannels);
output_DC = zeros(nlevels,nprobs,drnlchannels);
output_BM = zeros(nlevels,nprobs,drnlchannels);

for r = 1:nprobs % Loop over tone frequencies
    for j = 1:nlevels  % Loop over levels
        
      x = amp(j).*sin(2*pi*f_probe(r)*t); % Tone
      x = x.*w;
      x = [x zeros(1,N_zeros)];

      N_org= length(x);
    
      [~,outsig_nl,outsig_BM,fc] = IHC_output(x, fs, fc,'drnl', N_org, 'NH');
        
      output(j,r,:) = outsig_nl(:,1);
      output_AC(j,r,:) = abs(max(outsig_nl(1000:2000,:),[],1)-min(outsig_nl(1000:2000,:),[],1));
      output_DC(j,r,:) = abs(min(outsig_nl(1000:2000,:),[],1)+0.5*squeeze(output_AC(j,r,:)).'-outsig_nl(1,:));

      output_BM(j,r,:) = rms(outsig_BM);

    end
end


%% Plot

linesize       = 2.5;      % Set thickness of plot lines.
axislabelsize  = 10;       % Set size of axis labels (text).
axisnumbersize = 10;       % Set size of axis numbers.

idx = 1;  % index of drnl channel investigated
DC = squeeze(output_DC(:, :, idx));
BM = squeeze(output_BM(:, :, idx));

figure
hold on
yyaxis left
plot(levels, DC/1e-3,'linewidth', linesize);
ylabel('DC receptor potential [mV]')
set(gca, 'YScale', 'log')
yyaxis right
plot(levels,BM/1e-3,'linewidth', linesize);
set(gca, 'YScale', 'log')
hold off

legend('1','1.04','1.09','1.11')
xlim([10 100])
ylabel('BM velocity [mm/s]')
xlabel('input level [dB]')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 10])
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);
