function [Slop] = SlopEn(Sig, varargin)
    % SlopEn estimates the slope entropy of a univariate data sequence.
    %
    % [Slop] = SlopEn(Sig)
    % Returns the slope entropy (``Slop``) estimates for embedding dimensions 
    % [2, ..., ``m``] of the data sequence (``Sig``) using the default parameters:
    % embedding dimension = 2, time delay = 1, angular thresholds = [5 45], 
    % logarithm = base 2
    %
    % [Slop] = SlopEn(Sig, name, value, ...)
    % Returns the slope entropy (``Slop``) estimate of the data sequence (``Sig``) 
    % using the specified name/value pair arguments:
    %
    %   * ``m``     - Embedding Dimension, an integer > 1.
    %   * ``tau``   - Time Delay, a positive integer.
    %   * ``Lvls``  - Angular thresholds, a vector of monotonically increasing
    %     values in the range [0 90] degrees.
    %   * ``Logx``  - Logarithm base, a positive scalar (enter 0 for natural log).
    %   * ``Norm``  - Normalisation of ``Slop`` value, a boolean (default: true).
    %
    % References:
    %   [1] David Cuesta-Frau, "Slope Entropy: A New Time Series Complexity Estimator
    %   Based on Both Symbolic Patterns and Amplitude Information." Entropy 21.12 (2019): 1167.
    
    % Parse the inputs
    narginchk(1, 11);  % Validate number of input arguments
    Sig = squeeze(Sig);  % Ensure Sig is a vector

    p = inputParser;
    addRequired(p, 'Sig', @(x) isnumeric(x) && isvector(x) && (length(x) > 10));
    addParameter(p, 'm', 2, @(x) isnumeric(x) && isscalar(x) && (x > 1) && mod(x, 1) == 0);
    addParameter(p, 'tau', 1, @(x) isnumeric(x) && isscalar(x) && x > 0 && mod(x, 1) == 0);
    addParameter(p, 'Lvls', [5 45], @(x) isnumeric(x) && isvector(x) && all(x >= 0) && all(x <= 90));
    addParameter(p, 'Logx', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Norm', true, @(x) islogical(x));
    
    parse(p, Sig, varargin{:});
    
    % Extract parameters
    m = p.Results.m - 1;  % Decrement m to start from 2, as per the description
    tau = p.Results.tau;
    Lvls = p.Results.Lvls;
    Logx = p.Results.Logx;
    Norm = p.Results.Norm;

    if Logx == 0
        Logx = exp(1);  % Use natural logarithm if Logx is set to 0
    end

    % Calculate the slopes (derivative-like) of the signal
    Tx = atand(Sig(1 + tau:end) - Sig(1:end - tau));
    N = length(Tx);
    
    % Initialize variables
    Sx = zeros(N, m);
    Symbx = zeros(size(Tx));
    Slop = zeros(1, m);
    
    Lvls = sort(Lvls, 'ascend');  % Ensure angular thresholds are sorted in ascending order
    
    % Create symbolic representation of slopes
    for q = 2:length(Lvls)
        % Symbolic mapping of slopes based on thresholds
        Symbx(Tx <= Lvls(q) & Tx > Lvls(q - 1)) = q - 1;
        Symbx(Tx >= -Lvls(q) & Tx < -Lvls(q - 1)) = -(q - 1);

        % Assign extreme values for the last threshold
        if q == length(Lvls)
            Symbx(Tx > Lvls(q)) = q;
            Symbx(Tx < -Lvls(q)) = -q;
        end
    end

    % Compute slope entropy for each embedding dimension
    for k = 1:m
        Sx(1:N - k + 1, k) = Symbx(k:N);  % Assign the symbolic values to Sx matrix
        [~, ~, Locs] = unique(Sx(1:N - k + 1, 1:k), 'rows');  % Find unique symbolic patterns
        
        % Calculate the probability distribution of patterns
        if Norm
            p = accumarray(Locs, 1) / (N - k + 1);
            if round(sum(p)) ~= 1
                warning('Potential Error: Some permutations not accounted for!');
            end
        else
            p = accumarray(Locs, 1) / numel(accumarray(Locs, 1));
        end

        % Compute slope entropy for the current embedding dimension
        Slop(k) = -sum(p .* (log(p) / log(Logx)));
        clear Locs p;
    end
    
    % Display the final result
    disp('Slope Entropy:');
    disp(Slop);
end
