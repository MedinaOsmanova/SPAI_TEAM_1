function [t0_pos,s0_pos,t0_neg,s0_neg] = crossing_V7(S,t,level,imeth)

% [ind,t0,s0,t0close,s0close] = crossing_V6(S,t,level,imeth,slope_sign) % older format

% CROSSING find the crossings of a given level of a signal
%   ind = CROSSING(S) returns an index vector ind, the signal
%   S crosses zero at ind or at between ind and ind+1
%   [ind,t0] = CROSSING(S,t) additionally returns a time
%   vector t0 of the zero crossings of the signal S. The crossing
%   times are linearly interpolated between the given times t
%   [ind,t0] = CROSSING(S,t,level) returns the crossings of the
%   given level instead of the zero crossings
%   ind = CROSSING(S,[],level) as above but without time interpolation
%   [ind,t0] = CROSSING(S,t,level,par) allows additional parameters
%   par = {'none'|'linear'}.
%	With interpolation turned off (par = 'none') this function always
%	returns the value left of the zero (the data point thats nearest
%   to the zero AND smaller than the zero crossing).
%
%	[ind,t0,s0] = ... also returns the data vector corresponding to 
%	the t0 values.
%
%	[ind,t0,s0,t0close,s0close] additionally returns the data points
%	closest to a zero crossing in the arrays t0close and s0close.
%
%	This version has been revised incorporating the good and valuable
%	bugfixes given by users on Matlabcentral. Special thanks to
%	Howard Fishman, Christian Rothleitner, Jonathan Kellogg, and
%	Zach Lewis for their input. 

% Steffen Brueckner, 2002-09-25
% Steffen Brueckner, 2007-08-27		revised version

% Copyright (c) Steffen Brueckner, 2002-2007
% brueckner@sbrs.net

% M Noe
% added positive or negative slope condition

% check the number of input arguments
error(nargchk(1,4,nargin));

% check the time vector input for consistency
if nargin < 2 | isempty(t)
	% if no time vector is given, use the index vector as time
    t = 1:length(S);
elseif length(t) ~= length(S)
	% if S and t are not of the same length, throw an error
    error('t and S must be of identical length!');    
end

% check the level input
if nargin < 3
	% set standard value 0, if level is not given
    level = 0;
end

% check interpolation method input
if nargin < 4
    imeth = 'linear';
end


% make row vectors
t = t(:)';
S = S(:)';

% always search for zeros. So if we want the crossing of 
% any other threshold value "level", we subtract it from
% the values and search for zeros.
S   = S - level;

% first look for exact zeros
ind0 = find( S == 0 ); 

% then look for zero crossings between data points
S1 = S(1:end-1) .* S(2:end);
ind1 = find( S1 < 0 );

% bring exact zeros and "in-between" zeros together 
ind = sort([ind0 ind1]);

% and pick the associated time values
t0 = t(ind); 
s0 = S(ind);

if ~isempty(ind)

    if strcmp(imeth,'linear')
        % linear interpolation of crossing
        for ii=1:length(t0)
            %if abs(S(ind(ii))) > eps(S(ind(ii)))    % MATLAB V7 et +
            if abs(S(ind(ii))) > eps*abs(S(ind(ii)))    % MATLAB V6 et -    EPS * ABS(X)

                % interpolate only when data point is not already zero
                NUM = (t(ind(ii)+1) - t(ind(ii)));
                DEN = (S(ind(ii)+1) - S(ind(ii)));
                slope =  NUM / DEN;
                slope_sign(ii) = sign(slope);
                t0(ii) = t0(ii) - S(ind(ii)) * slope;
                s0(ii) = level;
            end
        end
    end

    % extract the positive slope crossing points 
    ind_pos = find(sign(slope_sign)>0);
    t0_pos = t0(ind_pos);
    s0_pos = s0(ind_pos);

    % extract the negative slope crossing points 
    ind_neg = find(sign(slope_sign)<0);
    t0_neg = t0(ind_neg);
    s0_neg = s0(ind_neg);

else

    % empty output
    ind_pos = [];
    t0_pos = [];
    s0_pos = [];

    % extract the negative slope crossing points 
    ind_neg = [];
    t0_neg = [];
    s0_neg = [];
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Addition:
% % Some people like to get the data points closest to the zero crossing,
% % so we return these as well
% [CC,II] = min(abs([S(ind-1) ; S(ind) ; S(ind+1)]),[],1); 
% ind2 = ind + (II-2); %update indices 
% 
% t0close = t(ind2);
% s0close = S(ind2);
end

