function [Out_PE, p,bu,patternS] = pe(x,d)
%Out_PE = pe(x,m,t)
% Calculate the permutation entropy (PE)
% This version calculates the NORMALIZED PE
%
% Input:    X: time series in vector form. if X is a matrix, each column is a different time series vector;
%           d: order of permuation entropy

% Output:
%           Out_PE:  Permutation Entropy
%
%Ref: 1)C. Bandt, and B. Pompe. "Permutation entropy: a natural complexity measure for time series." Physical review letters 88.17 (2002).
% 2) H. Azami and J. Escudero, “Improved Multiscale Permutation Entropy for Biomedical Signal Analysis: Interpretation and Application to Electroencephalogram Signals”,  Biomedical Signal Processing and Control , 2015.     %
%
% 190405V
% New features:
% Normalization
% Faster (no cycle) algorithm x10000
% 2D Matrix-capable

%% Calculations
% Vector orientation
n=length(x);
niv = n-(d-1);

%% Segmentation
% rescales x such that each column is a pattern
C = repmat((1:d)',1,n)+ones(d,1)*(0:n-1);
A=x(reshape(C(:,1:niv),1,d*(niv)));
A=reshape(A,d,niv);


%% Sorting Order
% Change the signal into pattern labels

[~,iv]=sort(A);     % Labels patterns
%% Histogram       
b= (10.^((d-1):-1:0))*iv;
bu=unique(b);
a=histc(b,bu);
%% Computation of entropy based on the histogram of the motifs
%a=a(a~=0);
p = a/niv;
P = p .* log(p);
P(isnan(P))=0;
Out_PE = -sum(P,2) / log(factorial(d));
[~,v]=max(p);
patternS=bu(v);

% [p,ip]=sort(p);
% bu=bu(ip);
[bu,ip]=sort(bu);
p=p(ip);
end