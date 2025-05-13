clear all
close all
clc

folder = '../Code_fNIRS/';
%% required paths 
addpath(folder);
addpath(strcat([folder 'data']));
addpath(strcat([folder 'Homer2cc']));
addpath(strcat([folder 'PE']));
addpath(strcat([folder 'PRSA']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load trained SVM and data
load('svmModel.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Load test data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separate test and training data
cv = cvpartition(Y, 'HoldOut', 0.2);
XTrain = X(training(cv), :); % 80% data -> training data
YTrain = Y(training(cv)); % 80% data -> training labels
XTest = X(test(cv), :); % 20% data -> test data
YTest = Y(test(cv)); % 20% data -> test labels
M = size(XTest,1);
N = size(XTest,2);
MM = size(X0,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Method parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Paramètres pour la méthode de Permutation Entropy
mpe = 2;% size embedding PE
mse = 15; %2; % size embedding SloPE
m = 3; %2; % size embedding
lag = 1; % sampling delay
ns = 1; % Normalisation (if applicable)
c = 6; % number of class for Dispersion Entropy
mmax = 10;  
%% Débruitage PRSA
d=4;% OP length
L=32;% PRSA length

%% Random Forest Parameters
numTrees = 100; % Tree number

%% Initialisation des méthodes
methods_name = {'PRSA',...
                'PE',...
                'SlopEn',...
                'rcdPE',...
                'SVM',...
                'RandomForest'
                };

methods_to_use = [1 2 3 4 6];  

nb_methods = length(methods_to_use);
MAE_out = zeros(nb_methods,1);

%% Entraînement du Random Forest
tic
rfModel = TreeBagger(numTrees, XTrain, YTrain, 'Method', 'classification');
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            SVM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
svm_predictions = predict(svmModel, XTest);
toc
% MAE_out(5) = mean(double(svm_predictions) ~= (double(YTest)-1)); % Taux d'erreur du Random Forest
MAE_out(5) = mean(svm_predictions ~= YTest); % Taux d'erreur (cohérent avec MAE)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                        Random Forest
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
rf_predictions = predict(rfModel, XTest);
toc
MAE_out(6) = mean(str2double(rf_predictions) ~= (double(YTest)-1));

%% Init
perm_entropy = zeros(MM,1);
slope_entropy = zeros(MM,mse-1);
rcdpe_entropy = zeros(mmax,MM);
Y = double(Y)-1;

for ind_met = 1:length(methods_to_use)
    MAE_tmp = zeros(MM,1);
    bd = zeros(MM,1);
    for n=1:MM
        x = (X0{n})';

        switch(methods_to_use(ind_met))
            case 1  %% PRSA 
                    try
                        % tic
                        [~, loc,~,~] = Pattern_prob(x, d);
                        prsa = PRSAfnirs(x, length(x), L, loc, 1);
                        if any(isnan(prsa))
                            MAE_tmp(n) = 1;
                            continue;
                        end
                        prsa_var = var(prsa);
                        prsa_threshold = 0.05;
                        if prsa_var >= prsa_threshold
                            bd(n) = 1;
                        end
                        % t_prsa(n) = toc;
                    catch ME
                        MAE_tmp(n) = 1;
                        continue;
                    end
                
            case 2  %% PE
                % tic
                perm_entropy(n) = pe(x, mpe);
                if perm_entropy(n) <= 0.7
                    bd(n) = 1;
                end
                % t_pe(n) = toc;
                
            case 3  %% SlopEn
                % tic
                slope_entropy(n,:) = SlopEn(x, 'm', mse, 'tau', lag);
                if mean(slope_entropy(n,:))>=0.7
                    bd(n)=1;
                end
                % t_spe(n) = toc;
            case 4 %% rcdPE 
                % tic
                rcdpe_entropy(:,n) = rcdpe_curve(x, m, mmax);
                if mean(rcdpe_entropy(:,n))<=0.7
                    bd(n)=1;
                end
                % t_rpe(n) = toc;
        end
        MAE_tmp(n) = abs(bd(n)-Y(n));
    end
    MAE_out(ind_met) = mean(MAE_tmp);
end

% Precision
P = 1-MAE_out


% mean(t_prsa)
% mean(t_pe)
% mean(t_spe)
% mean(t_rpe)
