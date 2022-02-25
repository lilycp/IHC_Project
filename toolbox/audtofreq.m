function freq = audtofreq(scale,aud);
%AUDTOFREQ  Converts auditory units to frequency (Hz)
%   Usage: freq = audtofreq(aud);
%  
%   AUDTOFREQ(scale,aud) converts values on the selected auditory scale to
%   values on the frequency scale measured in Hz.
%
%   See the help on FREQTOAUD to get a list of the supported values of the
%   scale parameter. 
% 
%   See also: freqtoaud, audspace, audfiltbw

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

%   AUTHOR: Peter L. Soendergaard
%   Last changed on $Date: 2009-01-14 17:46:50 +0100 (ons, 14 jan 2009) $
%   Last change occured in $Rev: 1 $

% ------ Checking of input parameters ---------

error(nargchk(2,2,nargin));

if ~isnumeric(aud) ||  all(aud(:)<0)
  error('%s: aud must be a non-negative number.',upper(mfilename));
end;

if ~ischar(scale)
  error('%s: the scale must be denoted by a character string.',upper(mfilename))
end;

% ------ Computation --------------------------
  
switch(lower(scale))
 case 'mel'
  freq = 700*(exp(aud/1127.01048)-1);
 case 'erb'
  freq = 228.8455*(exp(aud/9.265)-1);
 case 'bark'
  % This one was found through http://www.ling.su.se/STAFF/hartmut/bark.htm
  freq = 1960./(26.81./(aud+0.53)-1);
 case 'erb83'
  freq = 14363./(1-exp((aud-43.0)/11.7))-14675;
 otherwise
  error(['%s: unknown auditory scale: %s. Please see the help for a list ' ...
         'of supported scales.'],upper(mfilename),scale);
end;