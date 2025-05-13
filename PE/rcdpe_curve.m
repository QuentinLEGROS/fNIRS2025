function H = rcdpe_curve(x,d,mmax)
% 191015 M
% Refined Composite MPE curve vs m
% Creates all the MPE values for all scales up to mmax
%
% Input:
%       x:= original signal column vector
%       d:= embedded dimension
%       mmax := maximum time scale to analyze
%
% Output:
%       H:= column vector with the clasical MPE values for x
%
% Notes:
%       - Depends on function pe.m and cms_m1 in loop.
%       - Can I do it better without loops?
%
%% Computations
H = zeros(mmax,1);

for m=1:mmax
    xcoarse = resample(x,1,m);
    H(m) = rcpe(xcoarse,d);
end
