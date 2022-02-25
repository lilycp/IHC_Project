function [y,bw] = audspace(scale,flow,fhigh,n)
%AUDSPACE  Equidistantly spaced points on auditory scale
%   Usage: y=audspace(scale,flow,fhigh,n);
%
%   AUDSPACE(scale,flow,fhigh,n) computes a vector of length n containing values
%   equistantly scaled on the selected auditory scale between flow and fhigh. All
%   frequencies are specified in Hz.
%
%   See the help on FREQTOAUD to get a list of the supported values of the
%   scale parameter.
%  
%   [y,bw]=AUDSPACE(...) does the same but outputs the bandwidth between
%   each sample measured on the selected scale.
%  
% 
%   See also: freqtoaud, audspacebw, audfiltbw
% 
%   REFERENCES:
%     B. Glasberg and B. Moore. Derivation of auditory filter shapes from
%     notched-noise data. Hearing Research, 47(1-2):103, 1990.
%     

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
  
%   AUTHOR : Peter L. Soendergaard
%   Last changed on $Date: 2009-01-14 17:46:50 +0100 (ons, 14 jan 2009) $
%   Last change occured in $Rev: 1 $
  
% ------ Checking of input parameters ---------
  
error(nargchk(4,4,nargin));
  
% Default parameters.

if ~isnumeric(flow) || ~isscalar(flow) || flow<0
  error('%s: flow must be a non-negative scalar.',upper(mfilename));
end;

if ~isnumeric(fhigh) || ~isscalar(fhigh) || fhigh<0
  error('%s: fhigh must be a non-negative scalar.',upper(mfilename));
end;

if ~isnumeric(n) || ~isscalar(n) || n<=0 || fix(n)~=n
  error('%s: n must be a positive, integer scalar.',upper(mfilename));
end;

if flow>fhigh
  error('%s: flow must be less than or equal to fhigh.',upper(mfilename));
end;


% ------ Computation --------------------------

audlimits = freqtoaud(scale,[flow,fhigh]);

y = audtofreq(scale,linspace(audlimits(1),audlimits(2),n));  

bw=(audlimits(2)-audlimits(1))/(n-1);
