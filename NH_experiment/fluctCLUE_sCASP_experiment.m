%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script runs the experiment of additive noise as reportend in
% Jorgensen et al (2013) with the sCASP model
%
%  v1.0 - August, 2017. Helia Relaño Iborra
%  v2.0 - March, 2019. Helia Relaño Iborra
%  v3.0 - September, 2020. Helia Relaño Iborra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['start Jetal2013: ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]) % Report experiment start
%% -------- Initialization ---------- %%

load single_150_SentArray22kHz_varLength %C LUE material
sentenceFileLevel = -26; % Average sentence level

noise_names = {'SSN_CLUE_22kHz.wav','SSN_MOD_CLUE.wav','ISTS_eq.wav'}; % Noise files
conditions = {'SSN'  'SAM' 'ISTS'};

fs = 22050; % sampling frequency
Pref=20e-6; % level setting

speechSPL = 65; % Presentation level of speech

SNRs = -27:3:3; % SNRs to test

NSpeechsamples=2; % Number of sentences

%% ------------- Experiment ------------------ %%
%% Loop over sentences
% parfor q = 1:NSpeechsamples % Parallel processing
for q = 1:NSpeechsamples % Linear processing
    warning('off', 'all')
    
    %% Preallocate memory for variables:
    dfinal=zeros([length(SNRs),length(conditions)]);
    
    dint=cell([length(SNRs),length(SNRs)]);
    dsegments=cell([length(SNRs),length(SNRs)]);
    
    noise_scaled= cell(length(noise_names), length(SNRs));
    noise_glob= cell(length(noise_names), length(SNRs));
    
    %% Load sentence
    x = sentenceArray{q}';
    N_org= length(x);
    x = Pref*x*(1/rms(x))*10^((speechSPL)/20); % Set speech level
                                                      % sentences are presented at a level of 65 dB SPL.
    x=[x x]; % Prepan sentence
    clean_speech=x;
    
    Ts = 1/fs;
    T = length(x)/fs;
    t = 0:Ts:T;
    t = t(1:end-1);
    N = length(t);
    
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
                        
            [td_int_rep_speech, fc, td_cmod_cf] = sCASP_preprocessing(x', fs, 100, 8000,'drnl', N_org, 'NH'); % Template
            [td_int_rep_mix, fc, td_mod_cf] =  sCASP_preprocessing(test', fs, 100, 8000,'drnl', N_org, 'NH'); % Target
            
                    
            tmp = sCASP_backend(td_int_rep_speech, td_int_rep_mix, fs, fc, td_cmod_cf); % Back-end
            
            %Save the metric
            dfinal(k, n)=tmp.dfinal;
            dint{k, n}=tmp.dint;
            dsegments{k,n}=tmp.dsegments;
                        
            
        end %% End loop over SNRs
       
    end %% End loop over noise types
    
    disp(['sentence nr: ' num2str(q) ' ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
    
    %% Write variables:
     
    res_fluct_sCASP(q).dfinal = dfinal;
    res_fluct_sCASP(q).dint = dint;
    res_fluct_sCASP(q).dsegments = dsegments;
       
    end %% End loop over sentences

%% ------------- Save data ------------------ %%
current  = pwd;
name2=strcat( current, 'res_fluctCLUE_sCASP', '.mat');
save(name2, 'res_fluct_sCASP' )

