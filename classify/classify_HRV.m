clear
%%
%load('./data/features_labels.mat'); features1 = features; % For 270 sec, 69 features
load('./data/features&labels_time_resample.mat'); features1 = features; % For 270 sec, 74 features
load('./data/features_labels_MSE.mat'); features2 = features; % For 510 sec, 20 features
load('./data/features_labels_DFA.mat'); features3 = features; % For 330 sec, 7 features
% load('./toNengTai_HRV/features&labels_MMA.mat'); features4 = features; % For 510 sec, 54 features
load('./data/features_labels_PDFA.mat'); features5 = features;
N_sub = length(features);
load('./data/allsubLabel.mat');

%%
for sub = 1:N_sub
    N = size(allsubLabel{sub},1);
    features{sub} = cell(N,1);
    bad = 0;
    for ep = 9:N-9+1
        % nth beat index feature
        %qual = features1{sub}(ep,69);
        %tmp = [features1{sub}(ep,:) features2{sub}(ep,:) features3{sub}(ep, [2 4:end]) features5{sub}(ep,:)];% features4{sub}(ep-8,:)];
%         tmp = [features4{sub}(ep-8,:)];
        
        % time resample index feature
        qual = features1{sub}(ep,74);
        features1{sub}(ep,61:67) = old_features1{sub}(ep,61:67);
        tmp = [features1{sub}(ep,:) features2{sub}(ep,:) features3{sub}(ep, [2 4:end]) features5{sub}(ep,:)];% features4{sub}(ep-8,:)];
        if sum(isnan(tmp))==0 && qual>0.41
            features{sub}{ep} = tmp;
        else
            bad = bad + 1;
        end
    end
    % features{sub} = cell2mat(features{sub});
    disp([num2str(sub),':',num2str(bad)]);
end

%% 
% conf = cell(10,1);
% for tt = 1:10
    label = allsubLabel;
    
    alldata_tmp = cell(N_sub,1);
    for sub = 1:N_sub
        count = 0;
        alldata_tmp{sub} = cell(size(allsubLabel{sub},1),1);
        for i = 1:size(allsubLabel{sub},1)
            if isempty(features{sub}{i}) % || isempty(alldata{sub}{i})
                count = count+1;
                label{sub,1}(i,1) = -1;
            else
                alldata_tmp{sub}{i} = [features{sub}{i}];% alldata{sub}{i}];
            end
        end
        disp(['rubbish: ',num2str(count)]);
        alldata_tmp{sub} = cell2mat(alldata_tmp{sub});
        label{sub} = label{sub}(label{sub,1}>0);

        %% data upsample or downsample
        % SMOTE scheme
%         for stage = [1]
%             ls = find(label{sub}==stage);
%             if length(ls)>5
%                 [~,~,Xn,~] = smote(alldata_tmp{sub}(ls,:), 1);
%             elseif length(ls)>2
%                 [~,~,Xn,~] = smote(alldata_tmp{sub}(ls,:), 1, length(ls)-1);
%             else
%                 Xn = [];
%             end
%             alldata_tmp{sub} = [alldata_tmp{sub}; Xn];
%             label{sub} = [label{sub}; stage*ones(size(Xn,1),1)];
%         end
% 
%         for stage = [1 5]
%             ls = find(label{sub}==stage);
%             if length(ls)>5
%                 [~,~,Xn,~] = smote(alldata_tmp{sub}(ls,:), 2);
%             elseif length(ls)>2
%                 [~,~,Xn,~] = smote(alldata_tmp{sub}(ls,:), 1, length(ls)-1);
%             else
%                 Xn = [];
%             end
%             alldata_tmp{sub} = [alldata_tmp{sub}; Xn];
%             label{sub} = [label{sub}; stage*ones(size(Xn,1),1)];
%         end
        
        % sampling on N1 and N2 scheme
        if length(find(label{sub}==4))>50
            dws = randsample(find(label{sub}==4), 50); rest = find(label{sub}~=4);
            label{sub,1} = label{sub,1}([dws; rest], :);
            alldata_tmp{sub,1} = alldata_tmp{sub,1}([dws; rest],:);
        end
        if length(find(label{sub}==3))>50
            dws = randsample(find(label{sub}==3), 50); rest = find(label{sub}~=3);
            label{sub,1} = label{sub,1}([dws; rest], :);
            alldata_tmp{sub,1} = alldata_tmp{sub,1}([dws; rest],:);
        end
        % subjectId{sub} = sub*ones(size(alldata_tmp{sub},1),1);
    end
    %%
    % subjectId = subjectId(~any(cellfun('isempty', subjectId),2),:);
    alldata_tmp = alldata_tmp(~any(cellfun('isempty', alldata_tmp),2),:);
    label = label(~any(cellfun('isempty', label),2),:);
    t = cell2mat(label); 
    for ss = 1:5
        disp(length(find(t==ss)))
    end
    t = full(ind2vec(double(t(:)')))';
%%
    % subjectId for different number of folds
    N_sub = 8;
    subjectId = cell(N_sub, 1);
    for sub = 1:length(subjectId)
        subjectId{sub} = sub*ones(size(t,1)/N_sub, 1);
    end
    labelmat = cell2mat(label);
    labeltmp = buffer(labelmat, size(t,1)/N_sub);
    label = cell(N_sub,1);
    for l = 1:length(label)
        label{l} = labeltmp(:,l);
    end
    
    %% 5-class SVM
    sub = N_sub;
    predictors = cell2mat(alldata_tmp);
    predictors = zscore(predictors);
    % predictors = bsxfun(@minus, predictors, mean(predictors, 1));
    % [~, predictors] = pca(predictors, 'NumComponents', 300);


%     for i = 1:size(labelmat,1)
%         if labelmat(i,1) == 4; labelmat(i,1) = 3; end
%         if labelmat(i,1) == 5; labelmat(i,1) = 4; end
%     end
%% Print to csv for LGBM use 
%     csvwrite('./csvFiles/HRV_all.csv', predictors);
%     csvwrite('./csvFiles/HRV_all_label.csv', labelmat);
%     csvwrite('./csvFiles/subjectId.csv', cell2mat(subjectId));
%% 3/4 class svm
%     t(:,3) = t(:,3) + t(:,4);
%     t(:,4) = [];

    t(:,3) = t(:,3) + t(:,4) + t(:,5);
    t(:,[4 5]) = [];
    
    %% Classification
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
            
    template = templateSVM('KernelFunction', 'gaussian', ...
        'PolynomialOrder', [], 'KernelScale', 19, ...
        'BoxConstraint', 1, 'Standardize', true);
    
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
            'ClassNames', (1:size(t,2)).');   
        [pred, validationScores] = predict(Mdl, predictors(cvp.testStore{i}, :));
        Target{i} = pred;
        [~, cm{i}, ~, ~] = confusion(y(cvp.testStore{i}, :)', validationScores');
        toc
    
%         for j = 1:length(label{i})
%             if label{i}(j,1) == pred(j,1)
%                 storage{i,1} = [storage{i,1}; j validationScores(j,:)];
%             end
%         end
        disp(sum(diag(cm{i}))/sum(cm{i}(:)));
    end
    
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
% end