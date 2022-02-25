function [ratio] = IHC_ACDC_V02(f_probe,param)
% [ratio] = IHC_ACDC(f_probe,a)
%       Obtains the AC/DC ratio of the IHC receptor potential given the
%       probe frequencies. 
%
%   Input:
%       f_probe  - Frequencies of the input signals
%          a     - Parameters of the Boltzmann function 
%   Output:
%       ratio - AC/DC ratio of the receptor potential 


% PARAMETERS
flow=100;
fhigh= 8e3;
fc = erbspace(flow, fhigh, 60); % choose centre frequencies DRNL
nchannels = size(fc,2);

sbj='NH'; % DRNL profiles:'NH'/'OHCloss_5dB'/'OHCloss_10dB'/'OHCloss_20dB'/'OHCloss_30dB'/'HIx'

% Tone signal parameters
fs = 44100;             % Sampling frequency
Ts = 1/fs;
t_tot = 60e-3;          % Signal duration
Nsamp = t_tot*fs;       % Samples
ramp_dur = 5e-3;        % ramp duration in seconds
N_zeros = 20e-3*fs;     % zero-pad signal

t = (0:Nsamp-1)/fs;     % time vector

level = 80;   % tone level
amp = 20e-6.*10.^(level./20); % Amplitude in Pa;

% Ramp up/downn
window = hanning(2*floor(fs*ramp_dur))'; 
w1 = window(1:ceil((length(window))/2)); 
w2 = window(ceil((length(window))/2)+1:end); 
w = [w1 ones(1,Nsamp-length(w1)-length(w2)) w2]; 

nprobs= length(f_probe);

output = zeros(nprobs,Nsamp+N_zeros,nchannels);
output_AC = zeros(nprobs,nchannels);
output_DC = zeros(nprobs,nchannels);

for r = 1:nprobs % Loop over tone frequencies
        
      x = amp*sin(2*pi*f_probe(r)*t); % Tone
      x = x.*w;
      x = [x zeros(1,N_zeros)];
     
      b_hp = HeadphoneFilter(fs);     % calc headphone filtercoeffs
      b_me = middleearfilter_v2(fs);  % calc middle ear filtercoeffs

      inoutsig = filter(b_hp,1,x); % Outer-ear filterring
      inoutsig = filter(b_me,1,inoutsig); % middle-ear-ear filterring
    
     % DRNL filter
      outsig = zeros(length(inoutsig),nchannels); % channels are in colums
    for n = 1:nchannels
        outsig(:,n) = drnl_HI(inoutsig,fc(n),fs, sbj);
    end
    
      % IHC TRANSDUCTION 
       tau = 2.13e-3;
       C = 0.1;

       b = [0 tau*C];
       a_filt = [tau 1];
    
       [b_z,a_z] = bilinear(b,a_filt,pi);
       ster_disp = filter(b_z,a_z,outsig,[],1); % stereocilia displacement  
               
       % Conductance
        N = size(outsig,1); %length of signal

        Et = 100E-3;
        Ek = -70.45E-3;
        Omega = 40E-3;    % correction factor
        Ekp = Ek+Et*Omega;
        Gk = 18E-9;
        Cm = 6E-12;
        
        %a = [9.45E-9 52.7E-9 63.1E-9 29.4E-9 12.7E-9 0.33E-9]; % Lopez
        a = param;
        
        G = a(1)*1./((1+exp(-(ster_disp-a(2))/a(3))).*(1+exp(-(ster_disp-a(4))/a(5)))) + a(6);    
        G0 = a(1)*1./((1+exp(a(2)/a(3))).*(1+exp(a(4)/a(5)))) + a(6);
        
        % Intracellular Potential
        V_rest = (Gk*Ekp+G0*Et)./(G0+Gk);
        V_now = V_rest;
        
        V = zeros(N,nchannels);
        for ii=1:N
            V_now = V_now + (-G(ii,:).*(V_now-Et)-...
                Gk*(V_now-Ekp))*Ts/Cm;
            V(ii,:)=V_now;
         end


      output(r,:,:) = V;
      
      output_AC(r,:) = abs(max(squeeze(output(r,1000:2000,:)),[],1)-min(squeeze(output(r,1000:2000,:)),[],1));
      output_DC(r,:) = abs(min(squeeze(output(r,1000:2000,:)),[],1)+0.5*output_AC(r,:)-squeeze(output(r,1,:)).');

end

ACs = zeros(nprobs,1);
DCs = zeros(nprobs,1);
for r = 1:nprobs
    [~, BF_indx] = min(abs(fc-f_probe(r)));
    %BF_indx = 1;
    ACs(r) = output_AC(r,BF_indx);
    DCs(r) = output_DC(r,BF_indx);
end

ratio = ACs./DCs;

end

