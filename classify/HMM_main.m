% Main function of HMM learning
load('predictors&label.mat')
addpath('./HMM')
%% Data process
subId = cell2mat(subjectId);
% % 4-class AW/REM/N1/N2N3
% n_class = 4;
% label = labelmat;
% label(label == 5) = 4;

% 3-class AW/REM/N1~N3
n_class = 3;
label = labelmat;
label(label == 4 | label == 5) = 3;


%% Training
sub = 2;

    %% Descretize features
    % [symbols_train,max_idx,nmi_list] = Feature2sym(n_label,label(subId ~= sub),predictors(subId ~= sub));
       K = 10;
       symbols_train = kmeans(predictors(subId ~= sub),K);
% Generate initial parameters 
[TRANS_ems, EMIS_ems] = Generate_Para_in_ems(n_class, K, symbols_train,label(subId ~= sub));
% 4-class
%ini_ems = [0.6004    0.1321    0.1435    0.1241];
ini_ems = [0.6004    0.1321    0.1435];

% E-M alg training
[TRANS_F, EMIS_F, ini_F, i] = Training_Scaled(TRANS_ems, EMIS_ems , ini_ems, symbols_train');
%% 

%% Generate test symbols
mean_list = [];
predictors_train = predictors(subId ~=sub,:);
for i = 1:K
mean_list = [mean_list; mean(predictors_train(symbols_train == i,:),1)];
end

symbols_test = [];
for i = find(subId == sub)'
[~,d] = sort(vecnorm(ones([K,1])*predictors(i,:)-mean_list,2,2),'ascend');
symbols_test = [symbols_test;d(1)] ;
end
%% Test
Re = Viterbi_Scaled(TRANS_F, EMIS_F, ini_F, symbols_test', label(subId == sub)');
