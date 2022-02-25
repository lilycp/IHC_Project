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

close all, clear all, clc

addpath(genpath('/Users/liliana/Documents/PhD/MATLAB/IHC_project'))

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

f_probe=[300 600 1e3 2e3 4e3 6e3]; % Tone frequency
nprobs= length(f_probe);

levels = 10:10:100;  % tone levels
nlevels=length(levels);

amp = 20e-6.*10.^(levels./20); % Amplitude in Pa;

%% ------------- RUN MODEL ------------------ %%

for r = 1:nprobs % Loop over tone frequencies
    for j = 1:nlevels  % Loop over levels
        
      x = amp(j).*sin(2*pi*f_probe(r)*t); % Tone
      x = x.*w;
      x = [x zeros(1,N_zeros)];

      N_org= length(x);
    
      [outsig_lin,outsig_nl,outsig,fc] = IHC_output(x, fs, fc,'drnl', N_org, 'NH');

      output_OHC(j,r,:) = rms(outsig);

      outsig_lin_fft = fft(outsig_lin,[],1);
      output_lin(j,r,:)  = outsig_lin_fft(1,:);

      outsig_nl_fft = fft(outsig_nl,[],1);
      output_nl(j,r,:)  = outsig_nl_fft(1,:);

    end
end

%% Plot Freqs. -  Best frequency channel

linesize       = 1.5;      % Set thickness of plot lines.
axislabelsize  = 10;       % Set size of axis labels (text).
axisnumbersize = 10;       % Set size of axis numbers.

% colfreq=linspecer(6,'qualitative');
lin_s = {':'; '-.'; '-.'; '-'; '--'; '-'; '-.'; '-';'--'; '-'; '-'};
lin_w = [1, 1, 1 , 1,1, 2 ,2 ,2, 2,2, 2];
col_s= [.4 .4 .4;  .4 .4 .4; 0 0 0;...
    .4 .4 .4;  0 0 0; 0.3 0.3 0.3;...
    0.2 .2 .2; 0.4 0.4 0.4; 0 0 0;...
    0 0 0; 0.4 0.4 0.4];

i=1;

%IHC trans
figure
subplot(3,1,3)
for r = 1:nprobs
    [BF_dif, BF_indx] = min(abs(fc-f_probe(r)));
    plot(amp, squeeze(output_nl(:, r, BF_indx)),...
    'linewidth', 5,...
    'linestyle', lin_s{r},...
    'markerfacecolor', col_s(r, :));
    hold on
    str2{i}= sprintf('Probe freq. = %.0f Hz',f_probe(r));
    i = i+1;
end
grid on
ylabel('output level [MU]')
xlabel('input level [Pa]')
t = title('IHC non-lin');
set(t, 'fontweight', 'bold')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 10])
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

subplot(3,1,2)
for r = 1:nprobs
    [BF_dif, BF_indx] = min(abs(fc-f_probe(r)));
    plot(amp, squeeze(output_lin(:, r, BF_indx)),...
    'linewidth', 5,...
    'linestyle', lin_s{r},...
    'markerfacecolor', col_s(r, :));
    hold on
    str2{i}= sprintf('Probe freq. = %.0f Hz',f_probe(r));
    i = i+1;
end
grid on
ylabel('output level [MU]')
%le = legend(str2{:});
%set(le,'FontSize', 20,  'fontweight', 'bold','Location','northwest','box', 'off')
t = title('IHC lin');
set(t, 'fontweight', 'bold')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 10])
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

subplot(3,1,1)
for r = 1:nprobs
    [BF_dif, BF_indx] = min(abs(fc-f_probe(r)));
    plot(amp, squeeze(output_OHC(:, r, BF_indx)),...
    'linewidth', 5,...
    'linestyle', lin_s{r},...
    'markerfacecolor', col_s(r, :));
    hold on
    str2{i}= sprintf('Probe freq. = %.0f Hz',f_probe(r));
    i = i+1;
end
grid on
ylabel('output level [MU]')
%le = legend(str2{:});
%set(le,'FontSize', 20,  'fontweight', 'bold','Location','northwest','box', 'off')
t = title('OHC');
set(t, 'fontweight', 'bold')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 10])
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);