function aud = freqtoaud(scale,freq);
%FREQTOAUD  Converts frequencies (Hz) to auditory scale units.
%   Usage: aud = freqtoaud(scale,freq);
%
%   FREQTOAUD(scale,freq) converts values on frequency scale (measured in Hz) to
%   values on the selected auditory scale. The value of the parameter
%   scale determines the auditory scale:
%
%     'mel'  -  The mel scale is a perceptual scale of pitches judged by
%               listeners to be equal in distance from one another. The
%               reference point between this scale and normal frequency
%               measurement is defined by equating a 1000 Hz tone, 40 dB above
%               the listener's threshold, with a pitch of 1000 mels.
%               The mel-scale is defined in Stevens et. al (1937).
%
%     'bark'  - The bark-scale is originally defined in Zwicker (1961). A
%               distance of 1 on the bark scale is known as a critical
%               band. The implementation provided in this function is
%               described in Traunmuller (1990).
%
%     'erb'   - A distance of 1 erb is equal to the equivalent rectangular
%               bandwidth of the auditory filters at that point on the
%               frequency scale. The scale is normalized such that 0 erbs
%               corresponds to 0 Hz. The width of the auditory filters were
%               determined by a notched-noise experiment. The erb scale is
%               defined in Glasberg and Moore (1990).
%
%     'erb83' - This is the original defintion of the erb scale given in
%               Moore. et al. (1983).
% 
% 
%   See also: freqtoaud, audspace, audfiltbw
% 
%   REFERENCES:
%     B. Glasberg and B. Moore. Derivation of auditory filter shapes from
%     notched-noise data. Hearing Research, 47(1-2):103, 1990.
%     
%     B. Moore and B. Glasberg. Suggested formulae for calculating
%     auditory-filter bandwidths and excitation patterns. The Journal of the
%     Acoustical Society of America, 74:750, 1983.
%     
%     S. Stevens, J. Volkmann, and E. Newman. A scale for the measurement of
%     the psychological magnitude pitch. The Journal of the Acoustical
%     Society of America, 8:185, 1937.
%     
%     H. Traunmüller. Analytical expressions for the tonotopic sensory scale.
%     The Journal of the Acoustical Society of America, 88:97, 1990.
%     
%     E. Zwicker. Subdivision of the audible frequency range into critical
%     bands (frequenzgruppen). The Journal of the Acoustical Society of
%     America, 33(2):248-248, 1961. [1]http ]
%     
%     Henvisninger
%     
%     1. http://link.aip.org/link/?JAS/33/248/1

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

if ~isnumeric(freq) ||  all(freq(:)<0)
  error('%s: freq must be a non-negative number.',upper(mfilename));
end;

if ~ischar(scale)
  error('%s: the scale must be denoted by a character string.',upper(mfilename))
end;

% ------ Computation --------------------------

switch(lower(scale))
 case 'mel'
  aud = 1127.01048*log(1+freq/700);
 case 'erb'
  aud = 9.265*log(1+freq/228.8455);
 case 'bark'
  % The bark scale seems to have several different approximations available.
  
  % This one was found through http://www.ling.su.se/STAFF/hartmut/bark.htm
  aud = (26.81./(1+1960./freq))-0.53;
  
  % The one below was found on Wikipedia.
  %aud = 13*atan(0.00076*freq)+3.5*atan((freq/7500).^2);  
 case 'erb83'
  aud = 11.17*log((freq+312)./(freq+14675))+43.0;
 otherwise
  error(['%s: unknown auditory scale: %s. Please see the help for a list ' ...
         'of supported scales.'],upper(mfilename),scale);
end;