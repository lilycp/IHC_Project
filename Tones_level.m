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
%   - tone signal inputs
%   - AC & DC functions as a function of level
%   - Aim to reproduce figures from Dallos 1985 (Fig. 7)
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
fc = 1000;
drnlchannels = size(fc,2);

sbj='NH'; % DRNL profiles:'NH'/'OHCloss_5dB'/'OHCloss_10dB'/'OHCloss_20dB'/'OHCloss_30dB'/'HIx'

% Tone signal parameters
fs = 44100;             % Sampling frequency
t_tot = 60e-3;          % Signal duration
Nsamp = t_tot*fs;       % Samples
ramp_dur = 5e-3;        % ramp duration in seconds
N_zeros = 60e-3*fs;     % zero-pad signal

t = (0:Nsamp-1)/fs;         % time axis
f = (0:Nsamp-1)*fs/Nsamp;   % freq axis

%f_probe = [1 1.04 1.09 1.11]*fc(1);
f_probe = [1000 800 1500];
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
output_BM = zeros(nlevels,nprobs,Nsamp+N_zeros);

for r = 1:nprobs % Loop over tone frequencies
    for j = 1:nlevels  % Loop over levels
        
      x = amp(j).*sin(2*pi*f_probe(r)*t); % Tone
      x = x.*w;
      x = [x zeros(1,N_zeros)];

      N_org= length(x);
    
      [~,outsig_nl,outsig_BM,fc] = IHC_output(x, fs, fc,'drnl', N_org, 'NH');
        
      output(j,r,:) = outsig_nl(:,1);
      output_AC(j,r,:) = abs(max(outsig_nl(1500:2500,:),[],1)-min(outsig_nl(1500:2500,:),[],1));
      output_DC(j,r,:) = abs(min(outsig_nl(1500:2500,:),[],1)+0.5*squeeze(output_AC(j,r,:)).'-outsig_nl(1,:));

      output_BM(j,r,:) = outsig_BM; %rms(outsig_BM);

    end
end


%% Plot

linesize       = 2.5;      % Set thickness of plot lines.
axislabelsize  = 10;       % Set size of axis labels (text).
axisnumbersize = 10;       % Set size of axis numbers.

idx = 1;  % index of drnl channel investigated
DC = squeeze(output_DC(:, :, idx));
AC = squeeze(output_AC(:, :, idx));
BM = squeeze(output_BM(:, :, idx));

figure
subplot(1,2,1)
semilogy(levels, DC/1e-3,'linewidth', linesize);
ylabel('DC receptor potential [mV]')
xlabel('SPL [dB]')
v = [47 1e-4;72 1e-4;47 1e2;72 1e2];
f = [1 2 4 3];
patch('Faces',f,'Vertices',v,'FaceColor','red','FaceAlpha',.3);
legend('BF','below BF','above BF')
ylim([1e-4 1e2])
grid on
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

subplot(1,2,2)
semilogy(levels, AC/1e-3,'linewidth', linesize);
ylabel('AC receptor potential [mV]')
xlabel('SPL [dB]')
v = [47 1e-4;72 1e-4;47 1e2;72 1e2];
f = [1 2 4 3];
patch('Faces',f,'Vertices',v,'FaceColor','red','FaceAlpha',.3);
ylim([1e-4 1e2])
grid on
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);
