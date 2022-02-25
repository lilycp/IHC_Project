function [b,a]=gammatone(fc,fs,n,betain);
%GAMMATONE  Gammatone filter coefficients
%   Usage: [b,a] = gammatone(fc,fs,n,beta);
%          [b,a] = gammatone(fc,fs,n);
%          [b,a] = gammatone(fc,fs);
%
%   Input parameters:
%      fc    -  center frequency in Hz.
%      fs    -  sampling rate in Hz.
%      n     -  filter order.
%      beta  -  bandwidth of the filter.
%
%   Output parameters:
%      b     -  nominator coefficients.
%      a     -  denominator coefficients.
%
%   GAMMATONE(fc,fs,n,beta) computes the filter coefficients of a digital
%   gammatone filter with center frequency fc, bandwidth beta, order n and
%   sampling rate fs.
%
%   GAMMATONE(fc,fs,n) will do the same but choose a filter bandwidth
%   according to Glasberg and Moore (1990). The order n can only be 2 or 4.
%
%   GAMMATONE(fc,fs) will do as above for a 4th order filter.
%
%   If fc is a vector, each entry of fc is considered as one center
%   frequency, and the corresponding coefficients are returned as row
%   vectors in the output.
%
%   The inpulse response of the gammatone filter is given by
% 
%     g(t) = a*t^(n-1)*cos(2*pi*fc*t)*exp(-2*pi*beta*t)
% 
%   To create the filter coefficients of a 1-erb spaced filter bank using
%   gammatone filters use the following construction
% 
%     [b,a] = gammatone(fs,erbscalebw(flow,fhigh));
% 
%   REFERENCES:
%     A. Aertsen and P. Johannesma. Spectro-temporal receptive fields of
%     auditory neurons in the grassfrog. I. Characterization of tonal and
%     natural stimuli. Biol. Cybern, 38:223-234, 1980.
%     
%     B. Glasberg and B. Moore. Derivation of auditory filter shapes from
%     notched-noise data. Hearing Research, 47(1-2):103, 1990.

% Copyright (C) 2009 CAHR.
% This file is part of CASP version 0.01
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
%   AUTHOR : Stephan Ewert, Peter L. Soendergaard
%   Last changed on $Date: 2009-01-14 17:46:50 +0100 (ons, 14 jan 2009) $
%   Last change occured in $Rev: 1 $

% ------ Checking of input parameters ---------
  
% betain is used as input parameter in order not to mask the beta function.
  
error(nargchk(2,4,nargin));

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

if ~isnumeric(fc) || ~isvector(fc) || any(fc<0) || any(fc>fs/2)
  error(['%s: fc must be a vector of positive values that are less than half ' ...
         'the sampling rate.'],upper(mfilename));
end;

if nargin==4
  if ~isnumeric(betain) || ~isscalar(betain) || betain<=0
    error('%s: beta must be a positive scalar.',upper(mfilename));
  end;
end;

if nargin>2
  if ~isnumeric(n) || ~isscalar(n) || n<=0 || fix(n)~=n
    error('%s: n must be a positive, integer scalar.',upper(mfilename));
  end;
end;


% ------ Computation --------------------------

if nargin==2
  % Choose a 4th order filter.
  n=4;
end;

if nargin<4
  % Determine the correct multiplier for beta depending on the filter
  % order.
  switch(n)
   case 2
    betamul =  0.637;
   case 4
    betamul = 1.0183;
   otherwise
    error(['GAMMATONE: Default value for beta can only be computed for 2nd ' ...
           'and 4th order filters.']);
  end;
end;

nchannels = length(fc);

b=zeros(nchannels,1);
a=zeros(nchannels,n+1);

for ii = 1:nchannels

  if nargin<4
    % Determine betain parameter.
    betain = betamul*audfiltbw(fc(ii));
  end;
  
  btmp=1-exp(-2*pi*betain/fs);
  atmp=[1, -exp(-(2*pi*betain + i*2*pi*fc(ii))/fs)];
  
  b2=1;
  a2=1;
  
  for jj=1:n
    b2=conv(btmp,b2);
    a2=conv(atmp,a2);
  end
  
  % Place the result (a row vector) in the output matrices.
  b(ii,:)=b2;
  a(ii,:)=a2;

end;


