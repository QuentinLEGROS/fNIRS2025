function [Out_PE, p]  = rcpe(x,d)
% Refined Composite multiscale
% 150520L
% New method
% This works from the beginning of PE. It does not invoke pe() function
% 
% Refined Composite Multiscaling - 1º Moment (Average)
%
%   xcoarse:= Coarse-grained signal, already decomposed in all the composite versions.
%        The number of columns represent the scale of conversion.
%        Does not allow multiple signals.
%   d:= embedded dimension
%
% Output
%   Out_PE := Scalar value with the Refined Composite MPE value of the 
%   decomposed signal.
%% Calculations
% Vector orientation
n=size(x,1);
niv = n-(d-1);
nsig = size(x,2);
df=factorial(d);

%% Segmentation
% rescales x such that each column is a pattern

B=repmat(permute(x,[3,1,2]),d,1);
C = bsxfun( @plus, bsxfun(@plus,bsxfun(@plus,repmat((0:(d-1))',1,1,nsig),0:n)*d,(1:d)'), permute(n*d*[0:(nsig-1)],[1,3,2]));
A=B(C(:,1:niv,:));

%% Sorting Order
% Change the signal into pattern labels

[~,iA]=sort(A);     % Labels patterns
[~,iv]=sort(iA);     % Labels patterns

%% Histogram       

b = sum(bsxfun(@times,iv, (10.^((d-1):-1:0))' ));   %Convert patterns to decimal integers
c = b(:);
a = transpose( histc(c,unique(c)));                 % Counts number of patterns
a = sort(a,'descend');
if (df-size(a,2)) ~= 0
    a = [a,repmat(0,1,df-size(a,2))];
end

mn = length(c);

%% Computation of entropy based on the histogram of the motifs

p = a/mn;
P = p .* log(p);
P(isnan(P))=0;
Out_PE = -sum(P,2) / log(factorial(d));

% Number of forbidden patterns:= factorial(d) - size(a,2) - #zeros in each line
