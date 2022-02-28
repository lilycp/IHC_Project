function [in_out] = IHC_IO(amp,param)
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

%amp = 20e-6.*10.^(levels./20); % Amplitude in Pa;
nlevels= length(amp);

f_probe = 600;

% Ramp up/downn
window = hanning(2*floor(fs*ramp_dur))'; 
w1 = window(1:ceil((length(window))/2)); 
w2 = window(ceil((length(window))/2)+1:end); 
w = [w1 ones(1,Nsamp-length(w1)-length(w2)) w2]; 

for j = 1:nlevels % Loop over tone frequencies
        
      x = amp(j)*sin(2*pi*f_probe*t); % Tone
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
       %tau = param(1);
       C = 0.0251;
       %C = param(2);

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

   peak_pos(j,:) = max(V);
   peak_neg(j,:) = min(V);

end

idx = 48;
M = max(peak_pos(:,idx)) - V_rest;

% I/O Function
%in_out = [flip(peak_neg(:,idx)); peak_pos(:,idx)] - V_rest;
%in_out = [peak_pos(:,idx)] - V_rest;
in_out = [peak_pos(1:4,idx); peak_neg(5:end,idx)] - V_rest;
in_out = in_out/M;

end

