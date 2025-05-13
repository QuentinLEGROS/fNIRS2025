function [Slop] = SlopEn(Sig, varargin) 
% SlopEn  estimates the slope entropy of a univariate data sequence.
%
%   [Slop] = SlopEn(Sig) 
% 
%   Returns the slope entropy ('Slop') estimates for embedding dimensions
%   [2, ..., 'm'] of the data sequence ('Sig') using the default parameters: 
%   embedding dimension = 2, time delay = 1, angular thresholds = [5 45], 
%   logarithm = base 2 
% 
%   [Slop] = SlopEn(Sig, name, value, ...)
% 
%   Returns the slope entropy (``Slop``) estimate of the data sequence ('Sig')  
%   using the specified name/value pair arguments:
% 
%       * 'm'     - Embedding Dimension, an integer > 1. 
%         'SlopEn' returns estimates for each dimension [2, ..., 'm']
%       * 'tau'   - Time Delay, a positive integer
%       * 'Lvls'  - Angular thresolds, a vector of monotonically increasing
%         values in the range [0 90] degrees.
%       * ``Logx'  - Logarithm base, a positive scalar (enter 0 for natural log)
%       * 'Norm'  - Normalisation of ``Slop`` value, a boolean:
%               *  [false]  no normalisation
%               *  [true]   normalises w.r.t. the number of patterns found (default)
%
%   See also:
%       PhasEn, GridEn, MSEn, CoSiEn, SampEn, ApEn
%   
%   References:
%     [1] David Cuesta-Frau,
%           "Slope Entropy: A New Time Series Complexity Estimator Based on
%           Both Symbolic Patterns and Amplitude Information." 
%           Entropy 
%           21.12 (2019): 1167.
% 


Sig = Sig./sum(Sig);

narginchk(1,11)
Sig = squeeze(Sig);

p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
Chk3 = @(x) isnumeric(x) && isscalar(x) && (x > 1) && (mod(x,1)==0);
Chk4 = @(x) isnumeric(x) && (length(x) > 1) && (min(x)>0) && (max(x)<90);
addRequired(p,'Sig',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addParameter(p,'m',2,Chk3);
addParameter(p,'tau',1,Chk);
addParameter(p,'Lvls',[5 45],Chk4);
addParameter(p,'Logx',2,Chk2);
addParameter(p,'Norm',true,@(x) islogical(x));
parse(p,Sig,varargin{:})
m = p.Results.m; tau = p.Results.tau;  Lvls = p.Results.Lvls;
Logx = p.Results.Logx;  Norm = p.Results.Norm; 

if Logx == 0
    Logx = exp(1);
end
m = m-1;
Tx = atand(Sig(1+tau:end)-Sig(1:end-tau));
N = length(Tx);
Sx = zeros(N,m);
Symbx = zeros(size(Tx));
Slop = zeros(1,m);
Lvls = sort(Lvls,'ascend');
% Lvls = [0.1 4.5];

for q = 2:length(Lvls)
    Symbx(Tx<= Lvls(q) & Tx> Lvls(q-1)) = q-1;
    Symbx(Tx>=-Lvls(q) & Tx<-Lvls(q-1)) = -(q-1);
    
    if q == length(Lvls)
        Symbx(Tx> Lvls(q)) = q;
        Symbx(Tx<-Lvls(q)) = -q;
    end
end

for k = 1:m
    Sx(1:N-k+1,k) = Symbx(k:N);
    [~,~,Locs] = unique(Sx(1:N-k+1,1:k),'rows');
    
    if Norm
        p = accumarray(Locs,1)/(N-k+1);
        if round(sum(p))~= 1
            warning('Potential Error: Some permutations not accounted for!')
        end
    else
        p = accumarray(Locs,1)/numel(accumarray(Locs,1));
    end
    
    Slop(k) = -sum(p.*(log(p)/log(Logx)));
    clear Locs p
end

end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub