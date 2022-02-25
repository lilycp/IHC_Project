function out=filterbank(b,a,in,hopsize)
%FILTERBANK  Wrapper around filter to multiple filters
%   Usage: out=filterbank(b,a,in);
%          out=filterbank(b,a,in,hopsize);
%
%   FILTERBANK(b,a,in) filters the input signal with the filters
%   described in a and b.
%
%   FILTERBANK(b,a,in,hopsize) does the same, but only outputs every
%   hopsize sample in the time domain.
%
%   If a and b are matrices then each row corresponds to a subband
%   channel.
%
%   If f is a matrix then filtering is applied along the columns.
%
%   If both f, a and b are matrices the output will be
%   3-dimensional. First dimension is time, second dimension is frequency
%   channel third dimension corresponds to the column of the input signal.

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

error(nargchk(3,4,nargin));

if nargin==3
  hopsize=1;
end;

nchannels=size(b,1);

L=size(in,1);
W=size(in,2);

outlen=ceil(L/hopsize);

out=zeros(outlen,nchannels,W);

for ii=1:nchannels
  res = filter(b(ii,:),a(ii,:),in);
  res = res(1:hopsize:L,:);  
  out(:,ii,:) = reshape(res,outlen,1,W);
end;