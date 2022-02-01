function varargout = lin2db(varargin)
% DESCRIPTION y = lin2db(x1, x2, ...)
%  Converts linear to dB, no warning if input is zero.
%  The function accepts any number of arguments.
%  An error is issued for values less than zero.
% INPUT
%  x1  -- Real, Positive matrix
%  x2  -- Real, Positive matrix
%  x3  -- Real, ...
% OUTPUT
%  y1 -- As 1st input, but coverted into dB.
%  y2 -- As 2nd input, but coverted into dB.
%  y3 -- As 3rd input, but ...
% TRY
%  lin2db(0), lin2db([1 1]), lin2db([0 1],ones(2,2,2))
% SEE ALSO
%  db2lin

% by Magnus Almgren 990301 revised 050425

for i=1:nargin
if any(varargin{i}(:)<0)
  error('Input argument is less than zero')
end

% supress warning for log of 0 
warning off
varargout{i} = 10*log10(varargin{i}+0); % +0 enabler for operation on logical variables
warning on
end
% $Id: lin2db.m 28 2013-10-07 18:45:03Z maciel $
