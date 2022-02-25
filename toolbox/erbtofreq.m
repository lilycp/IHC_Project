function freq = erbtofreq(erb);
%ERBTOFREQ  Converts erb units to frequency (Hz)
%   Usage: freq = erbtofreq(erb);
%  
%   This is a wrapper around audtofreq that selects the erb-scale. Please
%   see the help on AUDTOFREQ for more information.
% 
%   See also: audtofreq, freqtoaud

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
  
freq = audtofreq('erb',erb);


