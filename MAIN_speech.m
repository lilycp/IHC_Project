%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Modelling Inner Hair Cell Transduction
%                                 Main
%                       Speech Signal Generation
%
%  References: 
%  Relaño-Iborra 2019: A speech-based computational auditory signal 
%                      processing and perception model
%  Sumner et al. 2002: A revised model of the IHC and auditory-nerve complex
%
%  Version 01
%   - revised & summarised from Special Project
%   - speech signal inputs
%   - plots outputs after OHC (DRNL) and IHC stage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all, clear all, clc

addpath(genpath('/Users/liliana/Documents/PhD/MATLAB/IHC_project'))

%% PARAMETERS

% DRNL frequencies
flow =100;
fhigh = 8e3;
fc = erbspace(flow, fhigh, 60);
drnlchannels = size(fc,2);

Nmod = 12; % number of mod filter frequencies

sbj='NH'; % DRNL profiles:'NH'/'OHCloss_5dB'/'OHCloss_10dB'/'OHCloss_20dB'/'OHCloss_30dB'/'HIx'

load single_150_SentArray22kHz_varLength % CLUE speech
sentenceFileLevel = -26;  % Average sentence level

Nsentences = 7; % Number of sentences to test

speechSPL = 60;   % Presentation level of the speech
Nspeech   = 27000;    % truncate speech signals to all be same size);
Nspl = length(speechSPL);
Pref= 20e-6;                    % Transformation to Pascals
fs = 22050;                     % Sampling freq

%% ------------- RUN MODEL ------------------ %%

% initialise
IHC_lin_av = 0;
IHC_nl_av = 0;
OHC_av = 0;

% Loop over sentences
   % parfor q = 1:Nsentences  % Parallel procesing
for q = 1:Nsentences % Normal processing

    % Load sentences
    speech  = sentenceArray{q};
    speech = speech(1:Nspeech);
    N_org= length(speech);     % Calculate length of original sentence
    % speech = [speech; speech]; % Prepane the same sentence
    
    speech = Pref*speech*(1/rms(speech))*10^((speechSPL)/20); % Set speech level
    N = length(speech);  % Overall speech size
                
    [outsig_lin,outsig_nl,outsig,fc] = IHC_output(speech, fs, flow, fhigh,'drnl', N_org, 'NH');

    IHC_lin_av = IHC_lin_av + outsig_lin;  % average over sentences
    IHC_nl_av  = IHC_nl_av + outsig_nl;
    OHC_av     = OHC_av + outsig;


  disp(['sentence nr: ' num2str(q) ' ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);  
end % end loop over sentences

IHC_lin_av = IHC_lin_av/Nsentences;
IHC_nl_av  = IHC_nl_av/Nsentences;
OHC_av     = OHC_av/Nsentences;


%% PLOT 
linesize       = 1.5;      % Set thickness of plot lines.
axislabelsize  = 10;       % Set size of axis labels (text).
axisnumbersize = 10;       % Set size of axis numbers.

t = (0:1:N-1)/fs; % time axis

figure
subplot(3,1,1)
plot(t,outsig(:,26))
ylabel('DRNL','Interpreter','Latex')
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);
grid on
xlim([0 t(end)])

subplot(3,1,2)
plot(t,outsig_lin(:,26))
ylabel('IHC lin','Interpreter','Latex')
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);
grid on
xlim([0 t(end)])

subplot(3,1,3)
plot(t,outsig_nl(:,26))
ylabel('IHC non-lin','Interpreter','Latex')
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);
grid on
xlim([0 t(end)])
