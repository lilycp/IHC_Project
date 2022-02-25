%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This scrip finds the free parameters of the logistic function that
% relates the output of the model (correlation coefficient d) with the
% percentage of correct answers. The data used for this fitting is taken
% from Nielsen and Dau (2009). A least squares analysis is used to fit
% the function to the data.
%
% v1.0 - August, 2017. Helia Relaño Iborra
% v2.0 - March, 2019. Helia Relaño Iborra
% v3.0 - September, 2020. Helia Relaño Iborra
%%%%%%%%%%%%%%%%%%%%%%%%%


%%  -------------- Initialization ------------ %%

Pcorrect_human = [0 8 35 71 90 100 ]; % Human data
SNRs= -8:2:2; % tested SNRs

load single_150_SentArray22kHz_varLength % CLUE speech
sentenceFileLevel = -26;  % Average sentence level

Nsentences=5; % Number of sentences to test

noise_name = 'SSN_CLUE_22kHz.wav'; % SSN Noise
speechSPL = 65; % Presentation level of the speech

Pref= 20e-6; % Transformation to Pascals

fs = 22050; % Sampling freq
res_fit_sCASP = struct([]); %Initialize results structure

%% ------------- Experiment ------------------ %%
%% Loop over sentences
% parfor q=1:Nsentences  % Parallel procesing
for q=1:Nsentences % Normal processing
    
    warning('off', 'all')
    %%  Preallocate memory for output metric:
    
    
    dfinal=zeros([1,length(SNRs)]);
    dint=cell([1,length(SNRs)]);
    dsegments=cell([1,length(SNRs)]);
    
    
    %% Load sentences
    speech  = sentenceArray{q};
    N_org= length(speech);  % Calculate length of original sentence
    speech = [speech; speech]; % Prepane the same sentence
    
    speech = Pref*speech*(1/rms(speech))*10^((speechSPL)/20); % Set speech level
    N = length(speech);  % Overall speech size
    
       
        
    %% Loop over SNRs
    for n=1:length(SNRs)
        
        %% Define noise token:
        noise = audioread(noise_name);
        Nsegments = floor(length(noise)/N); % number of possible random segments
        
        startIdx = randi(Nsegments-2 ,1)*N; % random initial index
        noise = noise(startIdx:startIdx+N -1)'; % select random segment from the noise file
        
        noise = Pref*(noise./rms(noise)*10^((speechSPL-SNRs(n))/20)); % sets the level of the noise signal
        
        if size(noise) ~= size(speech)
            noise = noise'; % Flips noise signal if needed
        end
        
        %% Build mixture:
        test = noise + speech; % Gerating the noisy mixture
        
        %% Run the model:
        
        %Preprocessing:
                [td_int_rep_speech, fc, td_cmod_cf] = sCASP_preprocessing(speech, fs, 100, 8000,'GT', N_org, 'NH');
                [td_int_rep_mix, fc, td_mod_cf] = sCASP_preprocessing(test, fs, 100, 8000,'drnl', N_org, 'NH');
        
              
        %Back-end:        
        tmp = sCASP_backend(td_int_rep_speech, td_int_rep_mix, fs, fc, td_cmod_cf);
        %Save the metric
        dfinal(n)=tmp.dfinal;
        dint{n}=tmp.dint;
        dsegments{n}=tmp.dsegments;
        
        
        
        
        
    end %% End of loop over SNRs
    disp(['sentence nr: ' num2str(q) ' ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
        
    %% Write into results structure:
    res_fit_sCASP(q).dfinal = dfinal;
    res_fit_sCASP(q).dint = dint;
    res_fit_sCASP(q).dsegments = dsegments;
    
end %% End loop over sentences

%% Save results:
current  = pwd;
name2=strcat( current, '/', 'res_CLUEfit_sCASP_5',  '.mat');
save(name2, 'res_fit_sCASP' )
