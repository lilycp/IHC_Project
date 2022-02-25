function [ratio] = IHC_ACDC(f_probe,param)
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
fc = 4000; % DRNL frequency

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
ster_disp = zeros(nprobs,N_zeros+Nsamp);
for r = 1:nprobs % Loop over tone frequencies
        
      x = amp*sin(2*pi*f_probe(r)*t); % Tone
      x = x.*w;
      x = [x zeros(1,N_zeros)];
     
      b_hp = HeadphoneFilter(fs);     % calc headphone filtercoeffs
      b_me = middleearfilter_v2(fs);  % calc middle ear filtercoeffs

      inoutsig = filter(b_hp,1,x); % Outer-ear filterring
      inoutsig = filter(b_me,1,inoutsig); % middle-ear-ear filterring
    
     % DRNL filter
      outsig = drnl_HI(inoutsig,fc,fs, sbj);
    
      % IHC TRANSDUCTION 
       %tau = 2.13E-3;   %low-pass filter applied to drnl output
       tau = param(1);
       %C = 0.1;
       C = param(2);

       b = [0 tau*C];
       a_filt = [tau 1];
    
       [b_z,a_z] = bilinear(b,a_filt,pi);
       ster_disp(r,:) = filter(b_z,a_z,outsig,[],1); % stereocilia displacement        
end

% Conductance
nfreq = size(ster_disp,1);
Nsamp = size(ster_disp,2);

Et = 100E-3;
Ek = -70.45E-3;
Omega = 40E-3;    % correction factor
Ekp = Ek+Et*Omega;
Gk = 18E-9;
Cm = 6E-12;

a = [9.45E-9 52.7E-9 63.1E-9 29.4E-9 12.7E-9 0.33E-9]; % Lopez

G = a(1)*1./((1+exp(-(ster_disp-a(2))/a(3))).*(1+exp(-(ster_disp-a(4))/a(5)))) + a(6);    
G0 = a(1)*1./((1+exp(a(2)/a(3))).*(1+exp(a(4)/a(5)))) + a(6);

% Intracellular Potential
V_rest = (Gk*Ekp+G0*Et)./(G0+Gk);
V_now = V_rest;

V = zeros(nfreq,Nsamp);
    for ii=1:Nsamp
            V_now = V_now + (-G(:,ii).*(V_now-Et)-...
                Gk*(V_now-Ekp))*Ts/Cm;
            V(:,ii)=V_now;
    end

AC = abs(max(V(:,1000:2000),[],2)-min(V(:,1000:2000),[],2));
DC = abs(min(V(:,1000:2000),[],2)+0.5*AC-V(:,1));

ratio = AC./DC;

end

