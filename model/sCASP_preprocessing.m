function [outsig_mf, fc, infpar] = sCASP_preprocessing(inoutsig, fs, flow, fhigh,BMtype, N_org, sbj, offset)
% The code is based on previous versions of authors: Torsten Dau, Morten
% Løve Jepsen, Boris Kowalesky and Peter L. Soendergaard
%
% Notes:
%
% The model has been optimized to work with speech signals, and the
% preprocesing and variable names follow this principle. The model is
% also designed to work with broadband signals. In order to avoid undesired
% onset enhancements in the adaptation loops, the model expects to recive a
% prepaned signal to initialize them.
%
%  Inputs:
%
%       * intsig_ss  :  clean speech template signal
%       * insig_ssnn :  noisy speech target signal
%       * fs         :  Sampling frequency
%       * flow       :  lowest center frequency of auditory filterbank
%       * fhigh      :  highest center frequency of auditory filterbank
%       * BMtype     :  model design: 'GT' gammatone (PEMO) or 'drnl'(CASP)
%       * N_org      :  length of original sentence
%       * sbj        :  subject profile for drnl definition
%
%  Outputs:
%       * out           : correlation metric structure inlcuding:
%            .dint      : correlation values for each modulation band
%            .dsegments : correlation values from each time window and mod. band.
%            .dfinal    : final (averaged) correlation
%
% REFERENCES:
%
%   Jepsen, M. L., Ewert, S. D., & Dau, T. (2008). A computational model
%   of human auditory signal processing and perception. Journal of the
%   Acoustical Society of America, 124(1), 422–438.
%
%   Relaño-Iborra, H., Zaar, J., & Dau, T. (2019). A speech-based computational
%   auditory signal processing and perception model (sCASP). The Journal of the
%   Acoustical Society of America, 146(5), 3306–3317.
%
%   -------------------------------
% 
%    v1.0 August 2017, Helia Rela?o Iborra:
% 
%    - Added input N_org (number of time samples in the original sentence)
%    and added cutting of prepanned sentence after adaptation loops.
%    - Changed the maximum modulation frequency for consistency with
%    Jepsen et al. 2008.
% 
%   -------------------------------
%   
%    v2.0 October 2018, Helia Relano Iborra:
%
%    - Added internal noise and removed minimum level check.
%    -------------------------------
%
%    v3.0 December 2018, Helia Relano Iborra:
%
%    - Changed drnl parameters for consistency with Jepsen (2008)
%    - Changed middle ear filter parameters from Lopez-Poveda & Meddis
%    (2001)
%    -------------------------------
%
%    v4.0 March 2019, Helia Relano Iborra:
%
%    - Adjust internal noise levels for NH free field thresholds
%    -------------------------------

%% find the center frequencies used in the filterbank, 1 ERB spacing
 if strcmp(BMtype,'GT')  
    %% ------- Case of linear auditory processing (gammatone)
    fc = erbspacebw(flow, fhigh, 1/2, 1000);
     % Calculate filter coefficients for the gammatone filter bank.
     [gt_b, gt_a]=gammatone(fc, fs);
     
     %% Apply the Gammatone filterbank
     outsig = 2*real(filterbank(gt_b,gt_a,inoutsig));
     
 else 
    %%  -------  Case of DRNL
    fc = erbspace(flow, fhigh, 60);
     nchannels = size(fc,2);
       
     %% Outer- and middle-ear filtering
     
    b_hp = HeadphoneFilter(fs); % calc headphone filtercoeffs
    b_me = middleearfilter_v2(fs);% calc middle ear filtercoeffs

    inoutsig = filter(b_hp,1,inoutsig); % Outer-ear filterring
    inoutsig = filter(b_me,1,inoutsig); % middle-ear-ear filterring
    
    %% DRNL filter
    outsig = zeros(length(inoutsig),nchannels); % channels are in colums
    for n = 1:nchannels
        outsig(:,n) = drnl_HI(inoutsig',fc(n),fs, sbj);
    end
 end

  
  %% Add noise to channels (freq specific):
    
offset = -3.7; % Correction for broadband levels:

int_noise_lvl = [-68.15;-73.05;-74.75;-75.25;-77.45;-77.75;-77.65;-78.15;...
    -79.15;-79.75;-79.25;-78.35;-77.75;-77.35;-77.05;-76.45;-75.85;-75.37;...
    -74.85;-74.25;-73.65;-73.65;-74.15;-75.35;-74.65;-75.15;-76.15;-73.05;...
    -74.55;-75.45;-76.45;-77.35;-77.95;-80.35;-79.75;-78.15;-78.85;-79.55;...
    -82.25;-79.25;-82.75;-81.45;-83.05;-85.05;-87.45;-90.15;-92.55;-91.55;...
    -90.85;-90.15;-89.45;-88.65;-89.65;-90.55;-90.95;-92.95;-94.45;-96.15;...
    -97.85;-99.65] + offset;  %% Internal noise (monochannel)
  
  
for m =1:size(outsig, 2)
int_noise = wgn(size(outsig, 1),1, int_noise_lvl(m));
outsig(:, m) = outsig(:, m) + int_noise;
end    
    

%% 'Haircell' envelope extraction
outsig = envextract(outsig,fs);

% IHCloss stage
if strcmp(sbj,'NH') == 0
         namestr = ['IHCloss_',sbj,'.mat']; load(namestr);
         for n = 1:nchannels
    [CFdiff,lookupnum] = min(abs(([50:25:10000]) - fc(n)));
    outsig(:, n) = outsig(:, n) .* out(lookupnum);
         end
end
 

if strcmp(BMtype,'GT')==0  % --- only for use with DRNL
    outsig = outsig * 10^(50/20); % Gain to compensate for the Outer/middle ear attenuation
  
%% Expansion (Auditory nerve) 
    outsig = outsig.^2; 
    
end  % --- end of only DRNL path
 
%% Non-linear adaptation loops:

    outsig = adaptloop(outsig, fs,10,2e-7); % lowest level in minlim.mat = 2e-7
       
%% Modulation processing:

    % set lowest mf as constant value. The multiplication by 0 is just an easy
    % way to get an array of zeros of the correct size.
    mflow = fc.*0;

    mfhigh= 1500; % Set maximum modulation center freq. to 1.5 kHz (v1.0 - HRI)
      
    [infpar, outsig_mf] = mfbtd(outsig,mflow,mfhigh,1,fs);	% MFB incl 150 LP
    outsig_mf = mfbtdpp(outsig_mf,infpar,fs);
   
     outsig_mf= outsig_mf((N_org+1):end, :, :);  % Removing prepanned sentence
end


