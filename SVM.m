clear all
close all
clc


% Path to nireject-benchmark (to download)
basePath = fullfile('../Code_fNIRS/nireject-benchmark', 'R22');

% Load the CSV file (contains labels, files name, etc)
csvFilePath = fullfile(basePath, 'R22.csv');
if ~isfile(csvFilePath)
    error('Fichier CSV introuvable : %s', csvFilePath);
end
csvData = readtable(csvFilePath);

% Does "decision" column exists (labels)?
if ~ismember('decision', csvData.Properties.VariableNames)
    error('Colonne "decision" absente du fichier CSV : %s', csvFilePath);
end

% Initialisation
X = []; % Training data
X0 = {};  % Training data without padding
Y = []; % Label

% Signals max length for padding
maxLength = 0;

% Upgradable : brose with a stepp of 44 -> 44 signals per patient
% A patient has only signals of equal size.

for i = 1:44:height(csvData) % Browse data
    fileName = csvData.filename{i}; % current file name
    matFilePath = fullfile(basePath, fileName); % current path
    if isfile(matFilePath) 
        matData = load(matFilePath); 
        try 
            signalLength = size(matData.Data.Raw.Hbo, 1);
            if signalLength > maxLength
                maxLength = signalLength; % Store max length
            end
        catch
            continue;
        end
    end
    
end

% Browse each raw of CSV file to load the signals
for i = 1:height(csvData)
    fileName = csvData.filename{i}; % file name
    decision = csvData.decision(i); % associated label
    channelIdx = csvData.channel(i);  % Considered channel

    % Label vérification
    if isempty(decision) || ~ismember(decision, [0, 1])
        warning('Label invalide à la ligne %d du fichier CSV %s', i, csvFilePath);
        continue;
    end

    matFilePath = fullfile(basePath, fileName);
    if ~isfile(matFilePath)
        warning('Fichier introuvable : %s', matFilePath);
        continue;
    end

    matData = load(matFilePath);
    try
        signalHbo = matData.Data.Raw.Hbo(:, channelIdx); % Load Hbo
        signalHbR = matData.Data.Raw.HbR(:, channelIdx); % Load HbR

        combinedSignal = signalHbo';

        % Padding
        paddedSignal = [combinedSignal, zeros(1, 2 * maxLength - length(combinedSignal))];

        X = [X; paddedSignal];  % Add signals
        Y = [Y; decision];
        
        X0{end+1} = combinedSignal;  % Add signals


    catch ME
        warning('Erreur lors de la lecture des signaux dans %s: %s', matFilePath, ME.message);
        continue;
    end
end

% Verify dimensions
if isempty(X) || isempty(Y)
    error('Les données X ou Y sont vides. Vérifiez les fichiers ou les labels.');
end

Y = categorical(Y);

% Separate test and training data
cv = cvpartition(Y, 'HoldOut', 0.2);
XTrain = X(training(cv), :); % 80% data -> training data
YTrain = Y(training(cv)); % 80% data -> training labels
XTest = X(test(cv), :); % 20% data -> test data
YTest = Y(test(cv)); % 20% data -> test labels

% Training SVM
tic
svmModel = fitcsvm(XTrain, YTrain, 'KernelFunction', 'linear');
toc
% Evaluation
predictions = predict(svmModel, XTest);
accuracy = mean(predictions == YTest);
fprintf('Précision du SVM : %.2f%%\n', accuracy * 100);

save('svmModel.mat', 'svmModel', 'X', 'Y', 'X0', 'maxLength');

