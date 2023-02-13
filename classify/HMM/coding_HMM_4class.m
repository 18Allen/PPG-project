function [coded_training_feature,coded_testing_feature] = coding_HMM(input,inputGT,testing_V)
% input: training feature
% inputGT: the sleep stage corresponding to the training feature
% testing_V: testing feature
% AW:1/REM:2/N1:3/N2N3:4


N = length(inputGT);
detail_GT = zeros(N,1);

for stage = 1:5
temp_inputGT = inputGT;
disconnected_interval = find(inputGT==stage);
starting = [disconnected_interval(1); disconnected_interval(1+find(diff(disconnected_interval)~=1))];

disconnected_interval = find(flipud(temp_inputGT)==stage);
ending = [disconnected_interval(1); disconnected_interval(1+find(diff(disconnected_interval)~=1))];
ending = flipud(N-ending+1);
if numel(starting)~=numel(ending)
    numel(starting)
    numel(ending)
    error('error')
end
segment_cell = cell(numel(starting),1);
for i = 1:numel(starting)
    segment_cell{i} = [starting(i):ending(i)]';
end

acc1_box = [];
acc2_box = [];
acc3_box = [];
    for i = 1:numel(starting)
       current_cell = segment_cell{i};
       size_cell = numel(current_cell); 
       if (size_cell>=5)
       acc1_box = [acc1_box; current_cell(1:1)];
       acc2_box = [acc2_box; current_cell(2:end-1)];
       acc3_box = [acc3_box; current_cell(end:end)];
       end
       if (size_cell==4)
       acc1_box = [acc1_box; current_cell(1:1)];
       acc2_box = [acc2_box; current_cell(2:3)];
       acc3_box = [acc3_box; current_cell(4:4)];
       end
       if (size_cell==3)
       acc1_box = [acc1_box; current_cell(1)];
       acc2_box = [acc2_box; current_cell(2)];
       acc3_box = [acc3_box; current_cell(3)];
       end
       if (size_cell==2)
       acc1_box = [acc1_box; current_cell(1)];
       acc3_box = [acc3_box; current_cell(end)];
       end
    end
    detail_GT(acc1_box) = (stage-1)*3+1;
    detail_GT(acc2_box) = (stage-1)*3+2;
    detail_GT(acc3_box) = (stage-1)*3+3;
end

noniso_index = find(detail_GT>0);
iso_index = find(detail_GT==0);

template = templateSVM('KernelFunction', 'gaussian', ...
    'PolynomialOrder', [], 'KernelScale', 10, ...
    'BoxConstraint', 1, 'Standardize', true);

    input_tem = Vectorize(input);
testing_V_tem = Vectorize(testing_V);

Sigma_input = cov(input_tem);
 Sigma_test = cov(testing_V_tem);

[U_i, ~] = eigs(Sigma_input, 50);
[U_t, ~] = eigs(Sigma_test, 50);

input = input_tem * U_i;
testing_V = testing_V_tem * U_t;

%%
if 1    
    Mdl = fitcecoc(input(noniso_index,:), detail_GT(noniso_index), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', [1:15]');
    
    [~, validationScores] = predict(Mdl, input(iso_index,:));
    [~,detail_GT(iso_index)] = max(validationScores,[],2);
end

if 0
A = find(inputGT==1);    
Mdl_Awake = fitcecoc(input(intersect(noniso_index,A),:), detail_GT(intersect(noniso_index,A)), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', [1:3]');
[~, validationScores] = predict(Mdl_Awake, input(intersect(iso_index,A),:));
[~,detail_GT(intersect(iso_index,A))] = max(validationScores,[],2);     


A = find(inputGT==2);    
Mdl_REM = fitcecoc(input(intersect(noniso_index,A),:), detail_GT(intersect(noniso_index,A)), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', [4:6]'); 
[~, validationScores] = predict(Mdl_REM, input(intersect(iso_index,A),:));
[~,detail_GT(intersect(iso_index,A))] = max(validationScores,[],2);         
detail_GT(intersect(iso_index,A)) = detail_GT(intersect(iso_index,A))+3;

A = find(inputGT==3);    
Mdl_N1 = fitcecoc(input(intersect(noniso_index,A),:), detail_GT(intersect(noniso_index,A)), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', [7:9]'); 
[~, validationScores] = predict(Mdl_N1, input(intersect(iso_index,A),:));
[~,detail_GT(intersect(iso_index,A))] = max(validationScores,[],2);  
detail_GT(intersect(iso_index,A)) = detail_GT(intersect(iso_index,A))+6;

    
A = find(inputGT==4);    
Mdl_N2 = fitcecoc(input(intersect(noniso_index,A),:), detail_GT(intersect(noniso_index,A)), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', [10:12]'); 
[~, validationScores] = predict(Mdl_N2, input(intersect(iso_index,A),:));
[~,detail_GT(intersect(iso_index,A))] = max(validationScores,[],2);  
detail_GT(intersect(iso_index,A)) = detail_GT(intersect(iso_index,A))+9;    
    
    
A = find(inputGT==5);    
Mdl_N3 = fitcecoc(input(intersect(noniso_index,A),:), detail_GT(intersect(noniso_index,A)), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', [13:15]');  
[~, validationScores] = predict(Mdl_N3, input(intersect(iso_index,A),:));
[~,detail_GT(intersect(iso_index,A))] = max(validationScores,[],2);  
detail_GT(intersect(iso_index,A)) = detail_GT(intersect(iso_index,A))+12;      
    
end
%%
    
    Md2 = fitcecoc(input(noniso_index,:), inputGT(noniso_index), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', [1:5]');
    
    [~, validationScores] = predict(Md2, testing_V);
    [~,testing_GT] = max(validationScores,[],2);

%%       
        A = find(detail_GT>=1 & detail_GT<=3);
        Md3_Awake = fitcecoc(input(A,:), detail_GT(A), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', [1:3]');
    
        A = find(detail_GT>=4 & detail_GT<=6);
        Md3_REM = fitcecoc(input(A,:), detail_GT(A), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', [4:6]');
    
    
        A = find(detail_GT>=7 & detail_GT<=9);
        Md3_N1 = fitcecoc(input(A,:), detail_GT(A), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', [7:9]');
    
         A = find(detail_GT>=10 & detail_GT<=12);
         Md3_N2 = fitcecoc(input(A,:), detail_GT(A), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', [10:12]');
    
        A = find(detail_GT>=13 & detail_GT<=15);
        Md3_N3 = fitcecoc(input(A,:), detail_GT(A), ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', [13:15]');
   
    

    
   coded_testing_V2 = zeros(length(testing_GT),1); 
   for sleep_stage = 1:5
       A = find(testing_GT==sleep_stage);
       
       
       if (sleep_stage==1)
        [~, validationScores] = predict(Md3_Awake, testing_V(A,:));
        [~,coded_testing_V2(A)] = max(validationScores,[],2);
       elseif (sleep_stage==2)
        [~, validationScores] = predict(Md3_REM, testing_V(A,:));
        [~,coded_testing_V2(A)] = max(validationScores,[],2);
        coded_testing_V2(A) = 3+coded_testing_V2(A);
       elseif (sleep_stage==3)
        [~, validationScores] = predict(Md3_N1, testing_V(A,:));
        [~,coded_testing_V2(A)] = max(validationScores,[],2);
        coded_testing_V2(A) = 6+coded_testing_V2(A);
       elseif (sleep_stage==4)
        [~, validationScores] = predict(Md3_N2, testing_V(A,:));
        [~,coded_testing_V2(A)] = max(validationScores,[],2);
        coded_testing_V2(A) = 9+coded_testing_V2(A);
       elseif (sleep_stage==5)
        [~, validationScores] = predict(Md3_N3, testing_V(A,:));
        [~,coded_testing_V2(A)] = max(validationScores,[],2);
        coded_testing_V2(A) = 12+coded_testing_V2(A);
       else
           error('error');
       end
   end
   
    coded_training_feature = detail_GT;   
    coded_testing_feature = coded_testing_V2;
    
end


