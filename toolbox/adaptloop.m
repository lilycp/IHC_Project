function inoutsig = adaptloop(inoutsig,fs,limit,minlvl,tau);
%ADAPTLOOP   Adaptation loops.
%   Usage: outsig = adaptloop(insig,fs,limit,minlvl,tau);
%          outsig = adaptloop(insig,fs,limit,minlvl);
%          outsig = adaptloop(insig,fs,limit);
%          outsig = adaptloop(insig,fs);
%
%   ADAPTLOOP(insig,fs,limit,minlvl,tau) applies non-linear adaptation to an
%   input signal insig sampled at a sampling frequency of fs Hz. limit is used
%   to limit the overshoot of the output, minlvl determines the lowest
%   audible threshhold of the signal and tau are time constants involved
%   in the adaptation loops. The number of adaptation loops is determined
%   by the length of tau.
%
%   ADAPTLOOP(insig,fs,limit,minlvl) does as above, but uses the values for
%   tau determined in Dau. 1996.
%
%   ADAPTLOOP(insig,fs,limit) does as above with a minimum threshhold minlvl
%   equal to 1e-5.
%
%   ADAPTLOOP(insig,fs) does as above with an overshoot limit of limit=10.
% 
%   REFERENCES:
%     T. Dau, D. Püschel, and A. Kohlrausch. A quantitative model of the
%     effective signal processing in the auditory system. I. Model structure.
%     The Journal of the Acoustical Society of America, 99(6):3615-3622,
%     1996a.
%     
%     D. Püschel. Prinzipien der zeitlichen Analyse beim Hören. PhD thesis,
%     Universität Göttingen, 1988.

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

% Copyright (c) 1999 - 2004 Stephan Ewert. All rights reserved.

%   AUTHOR : Stephan Ewert, Morten L. Jepsen, Peter L. Soendergaard
%   Last changed on $Date: 2009-01-14 17:46:50 +0100 (ons, 14 jan 2009) $
%   Last change occured in $Rev: 1 $

% ------ Checking of input parameters and default parameters ---------

error(nargchk(2,5,nargin));
  
% Default parameters for tau measured in seconds.
if nargin<5
  tau=[0.005 0.050 0.129 0.253 0.500];
else
  if ~isnumeric(tau) || ~isvector(tau) || tau<=0
    error('%s: tau must be a vector with positive values.',upper(mfilename));
  end;
end;

if nargin<4
  minlvl =1e-5;
else
  if ~isnumeric(minlvl) || ~isscalar(minlvl) || minlvl<=0
    error('%s: minlvl must be a positive scalar.',upper(mfilename));
  end;
end;

if nargin<3
  limit = 10;
else
  if ~isnumeric(limit) || ~isscalar(limit) || limit<=0
    error('%s: "limit" must be a positive scalar.',upper(mfilename));
  end;  
end;

% -------- Computation ------------------

inoutsig=comp_adaptloop(inoutsig,fs,limit,minlvl);
