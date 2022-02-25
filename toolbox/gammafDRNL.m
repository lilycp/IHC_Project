
% gammafDRNL.m - digital gammatone filter transfer function 
%
% Usage: [b,a] = gammaf1(fO,beta,n,type,fs)
%
% f0 = center frequency in Hz
% beta = 1/time constant
% n = gammatone filter order
% type = 1: bilinear transform, 2: impulse invariance
% fs = sampling rate in Hz
%
%
% b = [b0,b1, ... ,bN] = numerator coefficients
% a = [1,a1, ... ,aN] = denominator coefficients


function [b,a]=gammafDRNL(f0,beta,n,type,fs);

switch type
case 1
w0=2*pi*f0/fs;		
W0=tan(w0/2);		%prewarping
Beta=tan(pi*beta/fs);

tmp=1+1/Beta-i*W0/Beta;

b=[1, 1]/tmp;
a=[tmp, 1-1/Beta-i*W0/Beta]/tmp;

case 2

b=1-exp(-2*pi*beta/fs);
a=[1, -exp(-(2*pi*beta + i*2*pi*f0)/fs)];

case 3
O = 2*pi*f0/fs;
OI = 2*pi*beta/fs;
alfa = -exp(-OI)*cos(O);

b1 = 2*alfa;
b2 = exp(-2*OI);
a0 = abs((1+b1*cos(O)-i*b1*sin(O)+b2*cos(2*O)-i*b2*sin(2*O))/(1+alfa*cos(O)-i*alfa*sin(O)));
a1 = alfa*a0;
b=[a0, a1];
a=[1, b1, b2];
    
end
 
b2=1;
a2=1;

for j=1:n
	b2=conv(b,b2);
	a2=conv(a,a2);
end

b=b2;
a=a2;



