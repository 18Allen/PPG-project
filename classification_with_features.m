clear all
%% 
load('./data/allsubPPG_fs200_test.mat');
load('./data/allsubLabel_test.mat');
load('./data/features&labels.mat')
addpath('./lib')
%%
label = PPG_label;
data = features;
N_sub = size(data,1);
Channels = cell2mat(data);


n_class = 5;
if n_class == 4
% Use N1N2 as light, N3 as deep
for i=1:N_sub
    idx = label{i} == 4;
    label{i}(idx) = 3;
    idx =  label{i} == 5;
    label{i}(idx) = 4;
end
end
%% organization and indexing
Info.PID = [1:size(data,1)]';

t = cell(size(data, 1), 1);
t_save = t;
subjectId = cell(size(Info.PID));
for i = 1:size(label, 1)
    
    %id = find(label{i} > 0;
    % Sometime some class will be absent. "t" will have probelm
%     if isempty(find(label{i} == 1))
%         label{i} = [label{i}; 1];
%     end
%     if isempty(find(label{i} == 2))
%         label{i} = [label{i}; 2];
%     end
%     if isempty(find(label{i} == 3))
%         label{i} = [label{i}; 3];
%     end
%     if isempty(find(label{i} == 4))
%         label{i} = [label{i}; 4];
%     end
%     if isempty(find(label{i} == 5))
%         label{i} = [label{i}; 5];
%     end
    
    t{i} = full(ind2vec(double(label{i}(:)')))';
    %t{i} = t{i}(id, :);
    subjectId{i} = repelem(Info.PID(i), length(t{i}))';
    
    t_save{i} = t{i};
end


%% Load features to COM
COM = Channels;
%% n-class SVM
sub = size(data,1); % leave-one-out validatoin, this number should be #subject
predictors = COM;
y = cell2mat(t);
[~, response] = max(y, [], 2);

% create leave-one-out cross validation partition
cvp = struct;
cvp.NumObservations = size(response, 1); % the number of epoch
cvp.testSize = zeros(1, sub);
cvp.trainSize = zeros(1, sub);
cvp.testStore = cell(1, sub);
cvp.trainStore = cell(1, sub);

% create epoch-picking list
idx = cell(size(data,1),n_class);

for i = 1:sub
    cvp.testStore{i} = cell2mat(subjectId) == i;
    cvp.testSize(i) = sum(cvp.testStore{i});
    cvp.trainSize(i) = cvp.NumObservations - cvp.testSize(i);
    cvp.trainStore{i} = ~cvp.testStore{i};
end



template = templateSVM('KernelFunction', 'gaussian', ...
    'PolynomialOrder', [], 'KernelScale', 10, ...
    'BoxConstraint', 1, 'Standardize', true);

cm = cell(sub, 1);
classname = [1:n_class]';
for i = 1:sub
    disp(num2str(i))
    tic;
    Mdl = fitcecoc(predictors(cvp.trainStore{i}, :), response(cvp.trainStore{i}, :), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames',classname);%[1; 2; 3; 4; 5]);
    
    [tpp, validationScores] = predict(Mdl, predictors(cvp.testStore{i}, :));
    
    [~, cm{i}, ~, ~] = confusion(y(cvp.testStore{i}, :)', validationScores');
    toc
    
%     % Make index matrix of which epoch from the subject to pick
%     % Pass by the subject with < 50 valid epoch
%     if size(validationScores,1) < 50
%         for j = 1:5
%             idx{i,j} = [];
%         end
%         continue
%     end
%     % For those subjects with >= 50 valid epoch
%     bounds = sort(validationScores,1);
%     % First, pick the best 20 (Because it select w.r.t to score,
%     % so the number may not be exactly 20.
%     bounds = bounds(end-50,:);
%     for j = 1:5
%         %idx{i,j} = [idx{i,j},find(validationScores(:,j) >= bounds(j))'];
%         idx{i,j} = find(validationScores(:,j) >= bounds(j))';
%         idx{i,j} = idx{i,j}(1:50); %Second, pick 15 epoch out of 20;
%     end
end

% Calculating confusion matrix
ConfMat = zeros(n_class);
for i = 1:sub
    ConfMat = ConfMat + cm{i};
end
acc = sum(diag(ConfMat)) / sum(ConfMat(:));

disp(['Accuracy = ' num2str(acc)])
 %% Pick epoch
% selection = cell(5,4);
% for i = 1:5 % stage
%     for k=1:4 % channels .1,2: EEG; 3: ECG; 4: PPG
%         for j = 1:sub
%             % Only pick the last 30 seconds from each epoch
%             if k~=4
%                 selection{i,k} = [selection{i,k}; copy_channels{j,k}(end-(30*200 -1):end,idx{j,i})'];
%             else
%                 selection{i,k} = [selection{i,k}; copy_channels{j,k}(:,idx{j,i})'];
%             end
%        end
%     end
% end
% save('selection_A5','selection','ConfMat');
%% metrics

%ConfMat = ConfMat_fin;

SUM=sum(ConfMat,2);
nonzero_idx=find(SUM~=0);
normalizer=zeros(n_class,1);
normalizer(nonzero_idx)=SUM(nonzero_idx).^(-1);
matHMM=diag(normalizer)*ConfMat;
normalized_confusion_matrix = matHMM;

SUM=sum(ConfMat,1);
nonzero_idx=find(SUM~=0);
normalizer=zeros(n_class,1);
normalizer(nonzero_idx)=SUM(nonzero_idx).^(-1);
normalized_sensitivity_matrix=ConfMat*diag(normalizer);

recall = diag(normalized_confusion_matrix);
precision = diag(normalized_sensitivity_matrix);

F1_score = 2*(recall.*precision)./(recall+precision);
Macro_F1 = mean(F1_score);

TOTAL_EPOCH = sum(sum(ConfMat));
ACC = sum(diag(ConfMat))/TOTAL_EPOCH;
EA = sum(sum(ConfMat,1).*sum(transpose(ConfMat),1))/TOTAL_EPOCH^2;
kappa = (ACC-EA)/(1-EA);
if n_class == 5
output = cell(8, 9);
output(1, 2:end) = {'Predict-W', 'Predict-REM', 'Predict-N1', 'Predict-N2', 'Predict-N3', 'PR', 'RE', 'F1'};
output(2:6, 1) = {'Target-W', 'Target-REM', 'Target-N1', 'Target-N2', 'Target-N3'};
output(2:6, 2:6) = num2cell(ConfMat);
output(2:6, 7) = num2cell(precision);
output(2:6, 8) = num2cell(recall);
output(2:6, 9) = num2cell(F1_score);
output(8, 1:3) = {['Accuracy: ' num2str(ACC)], ['Macro F1: ' num2str(Macro_F1)], ['Kappa: ' num2str(kappa)]};
elseif n_class ==4
output = cell(7, 8);
output(1, 2:end) = {'Predict-W', 'Predict-REM', 'Predict-N1,N2', 'Predict-N3', 'PR', 'RE', 'F1'};
output(2:5, 1) = {'Target-W', 'Target-REM', 'Target-N1,N2', 'Target-N3'};
output(2:5, 2:5) = num2cell(ConfMat);
output(2:5, 6) = num2cell(precision);
output(2:5, 7) = num2cell(recall);
output(2:5, 8) = num2cell(F1_score);
output(7, 1:3) = {['Accuracy: ' num2str(ACC)], ['Macro F1: ' num2str(Macro_F1)], ['Kappa: ' num2str(kappa)]};
end
time = clock;
xlswrite(strcat('./','matrix-',num2str(time(4)), num2str(time(5)), '.xls'), output);