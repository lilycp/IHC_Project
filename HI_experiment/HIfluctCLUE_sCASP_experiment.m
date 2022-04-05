%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script runs the experiment of additive noise as reportend in
% Jorgensen et al (2013) with the sCASP model
%
%  v1.0 - August, 2017. Helia Relaño Iborra
%  v2.0 - March, 2019. Helia Relaño Iborra
%  v3.0 - Spetember, 2020. Helia Relaño Iborra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['start Christiansen2012: ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]) % Report experiment start
%% -------- Initialization ---------- %%

load single_150_SentArray22kHz_varLength %C LUE material
sentenceFileLevel = -26; % Average sentence level

noise_names = {'SSN_CLUE_22kHz.wav','SSN_MOD_CLUE.wav','ISTS_eq.wav'}; % Noise files
conditions = {'SSN'  'SAM' 'ISTS'};

fs = 22050; % sampling frequency

Pref=20e-6;

speechSPL = 80; % Presentation level of speech

templateSPL = 65; % PResentation level of the NH template

SNRs = -25:5:15; % SNRs to test

NSpeechsamples=25; % Number of sentences


load('sbj_IDs_christiansen2012_2ears.mat')
subjs = ['NH' sbj_IDs_christiansen2012_2ears]; % NH + HI subjets

Nsubjs = length(subjs);

%% ------------- Experiment ------------------ %%
%% Loop over sentences
% parfor q = 1:NSpeechsRamples % Parallel processing
for q = 1:NtSpeechsamples % Linear processing
    warning('off', 'all')
    
    %% Preallocate memory for variables:
    
    dfinal=zeros([Nsubjs , length(SNRs),length(conditions)]);
    dint=cell([Nsubjs , length(SNRs),length(SNRs)]);
    dsegments=cell([Nsubjs, length(SNRs),length(SNRs)]);
    
    noise_scaled= cell(Nsubjs ,length(noise_names), length(SNRs));
    noise_glob= cell(Nsubjs , length(noise_names), length(SNRs));
    
    %% Load sentence
    x = sentenceArray{q}';
    N_org= length(x);
    x = Pref*speech*(1/rms(x))*10^((speechSPL)/20);  % the level of the sentence is set such that the long-term RMS of all
    % sentences are presented at a level of 65 dB SPL.
    x=[x x]; % Prepan sentence
    clean_speech=Pref*(x/rms(x)*10^((templateSPL)/20));
    
    Ts = 1/fs;
    T = length(x)/fs;
    t = 0:Ts:T;
    t = t(1:end-1);
    N = length(t);
    
    %% Loop over subjects
    for s = 1:Nsubjs
        
        %% Loop over noise types (conditions)
        for n = 1:length(noise_names)
            
            [ tmp fs_tmp ] = audioread(noise_names{n});  % load the noise files
            if fs_tmp ~= fs
                noise_scaled{n} = resample(tmp,fs,fs_tmp);
            else
                noise_scaled{n} = tmp;
            end
            
            Nsegments = floor(length(noise_scaled{n})/N);
            
            
            %% Loop over SNRs
            for k = 1:length(SNRs)
                
                %% Set the noise level
                startIdx = randi(Nsegments-2 ,1)*N;
                noise = noise_scaled{n}(startIdx:startIdx+N -1); %Select random noise segment
                noise = Pref*(noise/rms(noise)*10^((speechSPL-SNRs(k))/20)); %Set noise level
                
                if size(noise) ~= size(x)
                    noise = noise';
                end
                test = noise + x;
                
                
                %% Run the model:
                
                %Preprocessing:
                [td_int_rep_speech, fc, td_cmod_cf] = sCASP_preprocessing(clean_speech', fs, 100, 8000,'drnl', N_org, 'NH'); % NH clean template
                [td_int_rep_mix, fc, td_mod_cf] =  sCASP_preprocessing(test', fs, 100, 8000,'drnl', N_org, subjs{s}); % HI noisy target
                
                % Back end:
                tmp = sCASP_backend(td_int_rep_speech, td_int_rep_mix, fs, fc, td_cmod_cf);
                %Save the metric
                dfinal(s,k, n)=tmp.dfinal;
                dint{s,k, n}=tmp.dint;
                dsegments{s,k,n}=tmp.dsegments;
                
                
                
            end  %% End loop over SNRs
            
        end  %% End loop over noise types
        
    end %% End loop over listener profile
    
    disp(['sentence nr: ' num2str(q) ' ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
    
    %% Write variables:
    
    res_fluct_sCASP(q).dfinal = dfinal;
    res_fluct_sCASP(q).dint = dint;
    res_fluct_sCASP(q).dsegments = dsegments;
    
    
end %% End loop over sentences


%% ------------- Save data ------------------ %%
current  = pwd;
name2=strcat( current, '/', 'res_fluct_CLUE_sCASP_HI',  '.mat');
save(name2, 'res_fluct_sCASP' , '-v7.3')

