function [OP, loc,bu,p]=Pattern_prob(x,d)
% this function determines the most probable OP  within the signal x and
% its localization

% Input:    x: time series in vector form. 
%           d: order of ordinal pattern

% Output:
%           OP:  most probable ordinal pattern
%           loc: ilts localization within x
%

%% Calculations
% Vector size
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
%% Computation of probability based on the histogram of the motifs
p = a/niv;
[~,v]=max(p);
OP=bu(v);
%% Localisation of the OP
loc=find(b==bu(v));
end