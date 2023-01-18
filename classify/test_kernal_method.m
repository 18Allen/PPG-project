%% Load the features & labels & subject partition
load('predictors&label.mat') % N_sub,labelmat,predictors,subjectId
%% Create Partition data for cross-validation
%cvp = cvpartition(cell2mat(subjectId),'Kfold',N_sub);
%% 3/4 class
% 3 class
N_class = 3;
label = labelmat;
label(labelmat == 5 | labelmat == 4) = 3;

% 4 class
% N_class = 4;
% label = labelmat;
% label(labelmat == 4) = 3;
% label(labelmat == 5) = 4;

%%  Set class index
    sub = N_sub;
    t = label; 
    for ss = 1:N_class
        disp(length(find(t==ss)))
    end
    t = full(ind2vec(double(t(:)')))';
%% Fit the model
rng('default')
% s1 = RandStream('mlfg6331_64','Seed',2718)
%Mdl = fitrkernel(predictors,label,'CVPartition',cvp,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus'))
template = templateKernel('Learner','svm','BlockSize',2e3,'Verbose',0,'KernelScale',5);
% One of four BinaryLearners
%              ClassNames: [-1 1]
%                   Learner: 'svm'
%   NumExpansionDimensions: 4096
%               KernelScale: 5.1315

% Mdl = fitcecoc(predictors,label, ...
%         'Learners', template, ...
%         'Coding', 'onevsone', ...
%         'ClassNames', (1:size(t,2)).',...
%         'CVPartition',cvp,...
%         'Verbose',1);
%%
% load('KernelModel_3class.mat') 
%%
cm = cell([sub,1]);
matsub = cell2mat(subjectId);
for i = 1:sub
    [pred, validationScores] = predict(Mdl.Trained{i}, predictors(matsub == i, :));
    [~, cm{i}, ~, ~] = confusion(t(matsub==i, :)', validationScores');
    disp(sum(diag(cm{i}))/sum(cm{i}(:)));
end
%%
% %
% % Mdl= fitcecoc(predictors,label,'Learner','kernel','OptimizeHyperparameters','auto','HyperparameterOptimizationOptions', ...
% %     struct('CVPartition',cvp,'Verbose',1,'UseParallel',true,'ShowPlots',false,'AcquisitionFunctionName','expected-improvement-plus'));
y = t;
[~, response] = max(y, [], 2);

% create leave-one-out cross validation partition
cvp = struct;
cvp.NumObservations = size(response, 1);
cvp.testSize = zeros(1, sub);
cvp.trainSize = zeros(1, sub);
cvp.testStore = cell(1, sub);
cvp.trainStore = cell(1, sub);
for i = 1:sub
    cvp.testStore{i} = cell2mat(subjectId) == i;
    cvp.testSize(i) = sum(cvp.testStore{i});
    cvp.trainSize(i) = cvp.NumObservations - cvp.testSize(i);
    cvp.trainStore{i} = ~cvp.testStore{i};
end

cm = cell(sub,1);
storage = cell(sub,1);
Target = cell(sub,1);
for i = 1:sub
    storage{i,1} = [];
    disp(num2str(i));
    tic;
    Mdl = fitcecoc(predictors(cvp.trainStore{i}, :), response(cvp.trainStore{i}, :), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'Verbose',1,...
        'ClassNames', (1:size(t,2)).');%,...
%         'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions', ...
%      struct('Verbose',1,'UseParallel',false,'ShowPlots',false,...
%      'MaxObjectiveEvaluations',400,'AcquisitionFunctionName','expected-improvement-plus'));   
    
    [pred, validationScores]   = predict(Mdl, predictors(cvp.testStore{i}, :));
    Target{i} = pred;
    [~, cm{i}, ~, ~] = confusion(y(cvp.testStore{i}, :)', validationScores');
    toc

    disp(sum(diag(cm{i}))/sum(cm{i}(:)));
end

%% statistics
   
    ConfMat = zeros(size(t,2));
    for i = 1:sub
        ConfMat = ConfMat + cm{i};
    end
%     conf{tt} = ConfMat;
    
    acc = sum(diag(ConfMat))/sum(ConfMat(:));
    stat = zeros(3,size(t,2));
    for i = 1:size(t,2)
        stat(1,i) = ConfMat(i,i)/sum(ConfMat(:,i));
        stat(2,i) = ConfMat(i,i)/sum(ConfMat(i,:));
        stat(3,i) = 2*stat(1,i)*stat(2,i)/(stat(1,i)+stat(2,i));
    end
    macroF1 = sum(stat(3,:))/size(stat,2);
    disp(['Accuracy = ',num2str(acc),'; macroF1 = ',num2str(macroF1)]);