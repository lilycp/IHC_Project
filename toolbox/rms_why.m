function y = rms_why(x)
%RMS RMS value of signal
%   Usage: y = rms(x)
%
%   RMS(x) computes the RMS (Root Mean Square) value of a finte sampled
%   signal sampled at a uniform sampling rate.
%
%   The RMS value of a signal x of length N is computed by
% 
%                     1  M
%      rms(x) = sqrt( - sum x(n)^2 )
%                     N n=1
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
  
error(nargchk(1,1,nargin));

if ~isnumeric(x) || ~isvector(x)
  error('RMS: Input must be a vector.');
end;

% ------ Computation --------------------------

% It is better to use 'norm' instead of explicitly summing the squares, as
% norm (hopefully) attempts to avoid numerical overflow. Matlab does not
% have 'sumsq'.
y = norm(x)/sqrt(length(x));
