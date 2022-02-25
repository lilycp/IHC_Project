%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Modelling Inner Hair Cell Transduction
%                                 Main
%                           Shamma reproduction
%
%  References: 
%  Relaño-Iborra 2019: A speech-based computational auditory signal 
%                      processing and perception model
%  Sumner et al. 2002: A revised model of the IHC and auditory-nerve complex
%  Shamma et al. 1986: A biophysical model of cochlear processing: 
%                       Intensity dependence of pure tone responses 
%
%  Version 01
%   - revised & summarised IHC transduction from Special Project
%   - Aim to reproduce Fig.5 in Shamma 1986
%   - sinusoidal stereocilia displacement directly inputted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all, clear all, clc

addpath(genpath('/Users/liliana/Documents/PhD/MATLAB/IHC_project'))

%% IHC TRANSDUCTION STAGE

% sinusoidal stereocilia displacement
freq = [20 200 500 1000 1500 2000];
w0_s = 2*pi*freq;
u = [0.05e-6 0.1e-6 0.3e-6 0.5e-6 0.8e-6];

fs = 44100;           % Sampling frequency
Ts = 1/fs;
t_tot = 100e-3;       % Signal duration
N = t_tot*fs;         % Samples
t = (0:N-1)/fs;       % time 

V_mat = zeros(N,length(u),length(freq));

params = 'Sumner';  % 1/0 depending on values used

for ff = 1:length(freq)
    f = freq(ff);
    w0 = w0_s(ff);

    for uu = 1:length(u)
    
    u_0 = u(uu);
    ster_disp = u_0*sin(w0*t);

switch params
    case 'Shamma'
        Et = 100E-3;
        Ek = -84E-3;
        Omega = 40E-3;    %correction factor
        Ekp = Ek+Et*Omega;
        G0 = 0.43e-8;
        Gmax = 0.15e-8;
        Gk = 1.078e-8;
        Cm = 1e-12;

        G = Gmax./(1+exp(-(xdata-a(2))/a(3)).*(1+exp(-(xdata-a(4))/a(5)))) + a(6);
    case 'Sumner'
        fs = 22050;    %sampling freq
        Ts = 1/fs;
        tau = 2.13E-3;
        C = 0.0251;
        
        G0 = 1.974E-9;
        Gmax = 8E-9;
        s0 = 85E-9;
        u0 = 7E-9;
        s1 = 5E-9;
        u1 = 7E-9;

        Et = 100E-3;
        Ek = -70.45E-3;
        Omega = 40E-3;    %correction factor
        Ekp = Ek+Et*Omega;
        Gk = 18E-9;
        Cm = 6E-12;
        
        Ga = G0-Gmax*1/(1+exp(u0/s0)*(1+exp(u1/s1)));       

        G = Gmax*1./(1+exp(-(ster_disp-u0)/s0).*(1+exp(-(ster_disp-u1)/s1))) + Ga;
        %G = Gmax*1./(1+exp(-(ster_disp-u0)/s0)) + Ga;

    case 'me'
        % boltzmann fit 
        fun_b   = @(a,xdata) a(1)*1./(1+exp(-(xdata-a(2))/a(3)).*(1+exp(-(xdata-a(4))/a(5)))) + a(6);
        a_b = [1.20594060102378e-08,8.69665159243074e-08,9.24125737653338e-08,3.08250738729188e-08,1.27205958125243e-08,-2.66214810891231e-11];
        G = fun_b(a_b, ster_disp); % apical conductance
    
        Et = 100E-3;
        Ek = -70.45E-3;
        Omega = 40E-3;    %correction factor
        Ekp = Ek+Et*Omega;
        G0 = a_b(6)+a_b(1)/((1+exp(a_b(2)/a_b(3)))*(1+exp(a_b(4)/a_b(5))));
        Gmax = a_b(1);
        Gk = 18E-9;
        Cm = 6E-12;
    end

    % Intracellular Potential
    V_rest = (Gk*Ekp+G0*Et)./(G0+Gk);
    V_now = V_rest;
    
    V = zeros(N,1);
    V(1) = V_rest;
        for ii=1:N-1
                V_now = V_now + (-G(ii).*(V_now-Et)-...
                    Gk*(V_now-Ekp))*Ts/Cm;
                V(ii+1)=V_now;
        end
    
        V_mat(:,uu,ff) = V;
        %%%%%%% CHECK THIS %%%%%%%%%%5
%     minim = V(1,:);
%     V = V-repmat(minim,N,1);
    
    % Low-pass filtering
    % f_cutoff = 1000;
    % [b, a] = butter(2, f_cutoff*2/fs);
    % V = filter(b,a, V);

    end
end


%% Plot output

linesize       = 1.5;      % Set thickness of plot lines.
axislabelsize  = 10;       % Set size of axis labels (text).
axisnumbersize = 10;       % Set size of axis numbers.

figure
plot(t/1e-3, squeeze(V_mat(:,:,2))/1e-3,'linewidth', linesize);
grid on
xlim([0 5])
xlabel('t in ms')
ylabel('Hair Cell Potential in mV')
legend('0.05$\mu$', '0.1$\mu$', '0.3$\mu$', '0.5$\mu$', '0.8$\mu$','Interpreter','latex','FontSize', axisnumbersize,  'fontweight', 'bold');
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

%%  peak-to-peak, AC-values

AC = squeeze(abs(max(V_mat,[],1)-min(V_mat,[],1)));

figure
plot(freq, AC/1e-3,'linewidth', linesize);
grid on
xlabel('frequency in Hz')
ylabel('Hair Cell Potential in mV')
legend('0.05$\mu$', '0.1$\mu$', '0.3$\mu$', '0.5$\mu$', '0.8$\mu$','Interpreter','latex','FontSize', axisnumbersize,  'fontweight', 'bold');
set(gca, 'FontSize', axisnumbersize,  'fontweight', 'bold', 'linewidth', linesize);

