function [V] = IHC_parameters(outsig,fs,nchannels,C)
% [V] = IHC_parameters(outsig,fs,nchannels,C)
%       Obtains the output of the inner hair cell transduction stage. This
%       is obtained via three main stages: 
%            1. Low-pass filter applied to BM velocity 
%            2. physiologically fitted function to convert to Conductance 
%            3. electrical circuit analogy to obtain receptor potential
%       This function also allows for a manipulation of the model
%       parameters. Here these are inputs
%   Input:
%       outsig - Output of the OHC stage
%         fs   - sampling frequency
%    nchannles - number of DRNL channels
%         C    - gain of the low-pass filter applied to BM velocity
%
%   Output:
%       V - Receptor Potential of the inner hair cells [ModelUnits]

N = size(outsig,1); %length of signal
Ts = 1/fs;

% Parameters of high-pass filter
tau = 2.13E-3;
%C = 0.0251; 
 %C = 0.1;

b = [0 tau*C];
a = [tau 1];

[b_z,a_z] = bilinear(b,a,pi);
ster_disp = filter(b_z,a_z,outsig,[],1); % stereocilia displacement

G0 = 1.974E-9;
Gmax = 8E-9;  % Gmax = 9.45E-9;
s0 = 85E-9;   % s0 = 63.1E-9;
u0 = 7E-9;    % u0 = 52.7E-9;
s1 = 5E-9;    % s1 = 12.7E-9;
u1 = 7E-9;    % u1 = 29.4E-9;

Et = 100E-3;
Ek = -70.45E-3;
Omega = 40E-3;    %correction factor
Ekp = Ek+Et*Omega;
Gk = 18E-9;
Cm = 6E-12;

Ga = G0-Gmax*1/((1+exp(u0/s0))*(1+exp(u1/s1))); 
% Ga = 0.33E-9;
% G0 = Ga + Gmax*1/((1+exp(u0/s0))*(1+exp(u1/s1)));  

G = Gmax*1./((1+exp(-(ster_disp-u0)/s0)).*(1+exp(-(ster_disp-u1)/s1))) + Ga;

% Intracellular Potential
V_rest = (Gk*Ekp+G0*Et)./(G0+Gk);
V_now = V_rest;

V = zeros(N,nchannels);
    for ii=1:N
            V_now = V_now + (-G(ii,:).*(V_now-Et)-...
                Gk*(V_now-Ekp))*Ts/Cm;
            V(ii,:)=V_now;
    end

end

