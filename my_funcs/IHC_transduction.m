function [V] = IHC_transduction(outsig,fs,nchannels)
% [V] = IHC_transduction(outsig,fs,nchannels)
%       Obtains the output of the inner hair cell transduction stage. This
%       is obtained via three main stages: 
%            1. Low-pass filter applied to BM velocity 
%            2. physiologically fitted function to convert to Conductance 
%            3. electrical circuit analogy to obtain receptor potential
%       (Finally a low-pass filter is applied.) -- currently not
%   Input:
%       outsig - Output of the OHC stage
%         fs   - sampling frequency
%    nchannles - number of DRNL channels
%
%   Output:
%       V - Receptor Potential of the inner hair cells [ModelUnits]

N = size(outsig,1); %length of signal
Ts = 1/fs;

% Parameters of high-pass filter
tau = 2.13E-3;
C = 0.0251;

b = [0 tau*C];
a = [tau 1];

[b_z,a_z] = bilinear(b,a,pi);
ster_disp = filter(b_z,a_z,outsig,[],1); % stereocilia displacement

Et = 100E-3;
Ek = -70.45E-3;
Omega = 40E-3;    %correction factor
Ekp = Ek+Et*Omega;
Gk = 18E-9;
Cm = 6E-12;

%a = [8E-9 7E-9 85E-9 7E-9 5E-9 1.2E-9]; %SOmner
a = [9.45E-9 52.7E-9 63.1E-9 29.4E-9 12.7E-9 0.33E-9]; % Lopez
%a = [9.45E-9 52.7E-9 63.1E-9 29.4E-9 12.7E-9 0.33E-9];

G = a(1)*1./((1+exp(-(ster_disp-a(2))/a(3))).*(1+exp(-(ster_disp-a(4))/a(5)))) + a(6) ; 
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

% Low-pass filtering
% f_cutoff = 1000;
% [b, a] = butter(2, f_cutoff*2/fs);
% V = filter(b,a, V);

end

