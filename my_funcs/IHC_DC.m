function [DC] = IHC_DC(levels,param)
% [DC] = IHC_DC(levels,a)
%       Obtains the DC component of the IHC receptor potential given 
%       different levels. 
%
%   Input:
%        levels  - Levels of the input signals
%          a     - Parameters of the Boltzmann function 
%   Output:
%       DC - DC component of the receptor potential 


% PARAMETERS
fc = [18000 20000]; % DRNL frequency
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

f_probe = [7e3 18e3];
nprobs= length(f_probe);

nlevels=length(levels);
amp = 20e-6.*10.^(levels./20); % Amplitude in Pa;

% Ramp up/downn
window = hanning(2*floor(fs*ramp_dur))'; 
w1 = window(1:ceil((length(window))/2)); 
w2 = window(ceil((length(window))/2)+1:end); 
w = [w1 ones(1,Nsamp-length(w1)-length(w2)) w2]; 

%% ------------- RUN MODEL ------------------ %%

% preallocation 
output_AC = zeros(nlevels,nprobs,nchannels);
output_DC = zeros(nlevels,nprobs,nchannels);

for r = 1:nprobs % Loop over tone frequencies
    for j = 1:nlevels  % Loop over levels
        
      x = amp(j).*sin(2*pi*f_probe(r)*t); % Tone
      x = x.*w;
      x = [x zeros(1,N_zeros)];
      N_org= length(x);
 
     % Outer- and middle-ear filtering
     b_hp = HeadphoneFilter(fs);     % calc headphone filtercoeffs
     b_me = middleearfilter_v2(fs);  % calc middle ear filtercoeffs

     inoutsig = filter(b_hp,1,x); % Outer-ear filterring
     inoutsig = filter(b_me,1,inoutsig); % middle-ear-ear filterring
    
    % DRNL filter
    outsig = zeros(N_org,nchannels); % channels are in colums
    for n = 1:nchannels
        outsig(:,n) = drnl_HI(inoutsig,fc(n),fs, sbj);
    end

    % Parameters of high-pass filter
%     tau = 2.13E-3;
%     C = 0.1;
    tau = param(1);
    C = param(2);

    b = [0 tau*C];
    a_filt = [tau 1];

    [b_z,a_z] = bilinear(b,a_filt,pi);
    ster_disp = filter(b_z,a_z,outsig,[],1); % stereocilia displacement
    
    Et = 100E-3;
    Ek = -70.45E-3;
    Omega = 40E-3;    %correction factor
    Ekp = Ek+Et*Omega;
    Gk = 18E-9;
    Cm = 6E-12;
    a = [8E-9 7E-9 85E-9 7E-9 5E-9 1.215E-9];
    
    G = a(1)*1./((1+exp(-(ster_disp-a(2))/a(3))).*(1+exp(-(ster_disp-a(4))/a(5)))) + a(6) ; 
    G0 = a(1)*1./((1+exp(a(2)/a(3))).*(1+exp(a(4)/a(5)))) + a(6);

    % Intracellular Potential
    V_rest = (Gk*Ekp+G0*Et)./(G0+Gk);
    V_now = V_rest;
    
    V = zeros(N_org,nchannels);
        for ii=1:N_org
                V_now = V_now + (-G(ii,:).*(V_now-Et)-...
                    Gk*(V_now-Ekp))*Ts/Cm;
                V(ii,:)=V_now;
        end

      outsig_nl = V;
        
      output_AC(j,r,:) = abs(max(outsig_nl(1000:2000,:),[],1)-min(outsig_nl(1000:2000,:),[],1));
      output_DC(j,r,:) = abs(min(outsig_nl(1000:2000,:),[],1)+0.5*squeeze(output_AC(j,r,:)).'-outsig_nl(1,:));

    end
end

DC = squeeze(output_DC(:, 2, 1))/1e-3;

end

