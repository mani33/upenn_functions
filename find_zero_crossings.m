function [ind,t0] = find_zero_crossings(S,t)

% make row vectors
t = t(:)';
S = S(:)';

% first look for exact zeros
ind0 = find( S == 0 ); 

% then look for zero crossings between data points
S1 = S(1:end-1) .* S(2:end);
ind1 = find( S1 < 0 );

% bring exact zeros and "in-between" zeros together 
ind = sort([ind0 ind1]);

% and pick the associated time values
n = length(ind1);
tt0 = nan(1,n);

for i = 1:n
    ai = [0 1] + ind1(i);
    s = S(ai);
    tt = t(ai);
    % Use slope intercept formula from highschool and solve for
    % interpolated t0    
    m = (s(2)-s(1))/(tt(2)-tt(1));
    b = s(1) - m*tt(1);
    tt0(i) = -b/m;
end

t0 = sort([t(ind0) tt0]);
if isempty(t0)
    disp('Could not find zero-crossing')
end