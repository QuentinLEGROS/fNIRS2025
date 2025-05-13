function prsa=PRSAfnirs(x,n,L,loc,Fs)

% this function determines the PRSA filter of order 2*L+1

% Input:    x: time series in vector form. 
%           N: length of x   


%% Localisation of window of length 2L+1, centererd on OP
loc=loc(find(loc>L & loc<n-L))-L;
Ln=2*L+1;
niv = n-(Ln-1);

%% Segmentation
% rescales x such that each column is a window
C = repmat((1:Ln)',1,n)+ones(Ln,1)*(0:n-1);
A=x(reshape(C(:,1:niv),1,Ln*(niv)));
A=reshape(A,Ln,niv);

prsa=mean(A(:,loc)'); prsa=prsa-mean(prsa); % PRSA signal

%% normalization
[H,f] = freqz(prsa,1,1024,Fs);    %  frequency response of the filter
scale = max(abs(H));             %  scaling factor
prsa = prsa/scale;
end