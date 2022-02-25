function [outsig_lin, outsig_nl,outsig,fc] = IHC_output(inoutsig, fs, fc,BMtype, N_org, sbj)
% 
% Notes:
%
% The model is heavily based on sCASP_preprocessing (Relaño-Iborra) and
% has been truncated to produce the output signal after the IHC stage of
% the auditory pathway. 
%
%  Inputs:
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
%   -------------------------------------------------------------------  %

 if strcmp(BMtype,'GT')  % Case of linear auditory processing (gammatone)   
     
     %fc = erbspacebw(flow, fhigh, 1/2, 1000);  % centre frequencies (1ERB spacing)
     nchannels = size(fc,2);
     [gt_b, gt_a]=gammatone(fc, fs); % filter coefficients
     
     % Apply the Gammatone filterbank
     outsig = 2*real(filterbank(gt_b,gt_a,inoutsig'));
     
 else % Case of DRNL
 
     %fc = erbspace(flow, fhigh, 60);
     nchannels = size(fc,2);
       
     % Outer- and middle-ear filtering
     b_hp = HeadphoneFilter(fs);     % calc headphone filtercoeffs
     b_me = middleearfilter_v2(fs);  % calc middle ear filtercoeffs

     inoutsig = filter(b_hp,1,inoutsig); % Outer-ear filterring
     inoutsig = filter(b_me,1,inoutsig); % middle-ear-ear filterring
    
    % DRNL filter
    outsig = zeros(length(inoutsig),nchannels); % channels are in colums
    for n = 1:nchannels
        outsig(:,n) = drnl_HI(inoutsig,fc(n),fs, sbj);
    end
 end
  
% INNER HAIR CELL STAGE
outsig_lin = envextract(outsig,fs);  % envelope extraction with Half-wave rectification and LP
outsig_nl  = IHC_transduction(outsig,fs,nchannels);   % non-linear IHC model


end


