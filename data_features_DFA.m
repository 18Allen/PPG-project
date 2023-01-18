%% Workspace Hygiene
clear;
%% Load data, tools
% mkdir './txt_data';
% for i = 1:30
%     writematrix(RRI_data{i}, ['sub_',num2str(i),'.txt'], 'Delimiter', 'space');
% end
%%
% For ECG
% load('./data/ECG/allsubECG_fs200.mat');
% load('./data/ECG/allsubLabel.mat');
% load('./data/ECG/processed_data.mat');
% For PPG
load('./data/allsubPPG_fs200.mat');
load('./data/allsubLabel.mat');
load('./data/processed_data.mat');
% addpath('../SST_TF_analysis/TF_anaylsis');
addpath('./lib')
%% Data parsing--------------------------------------
%% Set parameters
N_sub = size(allsubLabel,1);
PPG_data = cell(size(allsubLabel));
PPG_label = cell(size(allsubLabel));
PPG_label_index = cell(size(allsubLabel));

fs = 200;
len_epoch = 330; %in second, 5 min
len_orig = 30;

% basic parameters for STFT 
sampling_rate = 200; 

% Loop parameters
n_class = 5;
stats = cell([n_class,3]);

% HRV
upsampling_rate = 500; % For PPG_peak_detection

% Generated feature
features = cell([n_class,1]);
n_features = 7;
% 1~3 DFA, 4~7 WDFA
% (See 'data_features_270s.m') 1~44 Traditional time, 45~54 Traditional freq (2 HFpole features left),
% 55 ApEn, 56 Higuchi fractal dimension, 57~58 teager energy

% Temp function
% extract = cell([n_class,1]);
%% Parse the data to epoch interval
IHR_list = cell(size(allsubLabel));
RRI_list = cell(size(allsubLabel));
% WDFA_list = cell(size(allsubLabel));
% WDFA_curves = cell(size(allsubLabel)); %four combination n and sigma
RRI_res_list = cell(size(allsubLabel));
s = round(len_epoch/len_orig);
for i =1:N_sub
    i
    slabel = ceil(s/2);
    % PPG
    PPG_data{i} = buffer(allsubPPG{i},fs*len_epoch,fs*(len_epoch-len_orig))';
    % PPG_data{i} = buffer(allsubECG{i},fs*len_epoch,fs*(len_epoch-len_orig))';
    PPG_data{i} = PPG_data{i}(s:end,:);
    %Label
    
    PPG_label{i} = allsubLabel{i};
    %Index (Unused = 0, Used = 1)
    index = ones(size(allsubLabel{i}));
    index([1:slabel-1 end-slabel+2:end]) = 0;
    PPG_label_index{i} = index;
     
    %IHR
    IHR_list{i} = buffer(IHR_data{i},len_epoch*4,(len_epoch-len_orig)*4)';
    IHR_list{i} = IHR_list{i}(s:end,:);
    %TBD temp. In correspond to IHR(end-12:end) =0 in data_parsing.m
    %IHR_list{i}(1,1:12) = 1;IHR_list{i}(end,end-12:end) = 1;
    
%     % Get WDFA curve from RRI_res_data
    RRI_res = RRI_res_data{i};
    WDFA_curves{i} = cell([4,1]);
    WDFA_curves{i}{1} = WDFA_fun(RRI_res,30,30,1);
    WDFA_curves{i}{2} = WDFA_fun(RRI_res,30,90,1);
    WDFA_curves{i}{3} = WDFA_fun(RRI_res,90,30,1);
    WDFA_curves{i}{4} = WDFA_fun(RRI_res,90,90,1);
    
    WDFA_list{i} = cell([4,1]);
    for j = 1:4
        WDFA_list{i}{j} = buffer(WDFA_curves{i}{j},len_epoch,(len_epoch-len_orig))';
        WDFA_list{i}{j} = WDFA_list{i}{j}(s:end,:);
    end
    
    % RRI_res_list for teager energy
    RRI_res_list{i} = buffer(RRI_res,len_epoch,(len_epoch-len_orig))';
    RRI_res_list{i} = RRI_res_list{i}(s:end,:);
    
    % RRI
    RRI_list{i} = cell(size(PPG_label{i}));
    for j = 1:length(PPG_label{i})
        idx = find(locs_data{i} > (-(len_epoch-len_orig)/2 + (j-1)*len_orig)*upsampling_rate & locs_data{i} <= ((len_epoch-len_orig)/2+(j)*len_orig)*upsampling_rate);
        RRI_list{i}{j} = RRI_data{i}(idx(1:end-1));
        
        % Apply change if less than 300 RRI in an epoch
%         if length(RRI_list{i}{j}) < 300
%             RRI_list{i}{j} = interp1(1:length(RRI_list{i}{j}),RRI_list{i}{j},...
%                 [linspace(1,(300-length(RRI_list{i}{j}))+1,(300-length(RRI_list{i}{j}))*2+1 ),(300-length(RRI_list{i}{j}))+2:length(RRI_list{i}{j})]);
%         end 
    end
    RRI_list{i} = RRI_list{i}(slabel:end-slabel+1); % notice: start from slabel
end

%%  Analysis
for l = 1:N_sub
    features{l} = zeros([size(allsubLabel{l},1), n_features]);

%     qual = Quality{l};
    for m = find(PPG_label_index{l} > 0).'
        % m
%         if length(RRI_list{l}{m})<165 %|| ~(qual(m-floor(s/2))<prctile(qual,95))
%             features{l}(m,:) = nan;
%             continue;
%         end
        RRI = RRI_list{l}{m-slabel+1};
    
        %% DFA features
        % select scale for the whole DFA
        Q = exp(linspace(log(10),log(300),37)); 
        pts = round(Q);
        order = 2; % In the paper they use quadratic one
        F = DFA_fun(RRI,pts,order);
    
        % A: a 2x1 vector. A(1) is the scaling coefficient "alpha",
        % A(2) the intercept of the log-log regression, useful for plotting (see examples).
        A = polyfit(log(pts),log(F)',1);
        alphaAll = A(1) ;   
    
        idx = find(Q>=10 & Q<=40) ;
        A = polyfit(log(pts(idx)),log(F(idx))',1);
        alpha1 = A(1) ;
    
        idx = find(Q>=70 & Q<=300) ;
        A = polyfit(log(pts(idx)),log(F(idx))',1);
        alpha2 = A(1) ;
        features{l}(m,1) = alphaAll;
        features{l}(m,2) = alpha1;
        features{l}(m,3) = alpha2;
        
        %% PDFA
        % P = PDFA(RRI(end-63:end),64,1);
%         p =  polyfit(linspace(0,1,64), P, 1);
%         features{l}(m,54) = p(1);
        
        %% WDFA feature
        WDFAs = cell([4,1]);
        for j = 1:4
            WDFAs{j} = WDFA_list{l}{j}(m-slabel+1,:);
            features{l}(m,3+j) = max(WDFAs{j});
        end

    end
%     save('./data/ECG/features&labels_DFA.mat','features','PPG_label','PPG_label_index');
    save('./features&labels_DFA.mat','features','PPG_label','PPG_label_index');

end
clear WDFAs