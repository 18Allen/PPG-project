%% Workspace Hygiene
%clear;
%% Load data, tools
% mkdir './txt_data';
% for i = 1:30
%     writematrix(RRI_data{i}, ['sub_',num2str(i),'.txt'], 'Delimiter', 'space');
% end
%%

load('./data/allsubPPG_fs200.mat');
load('./data/allsubLabel.mat');
load('./data/processed_data.mat')
%addpath('../SST_TF_analysis/TF_anaylsis');
addpath('./lib')
%% Data parsing--------------------------------------
%% Set parameters
N_sub = size(allsubLabel,1);
PPG_data = cell(size(allsubLabel));
PPG_label = cell(size(allsubLabel));
PPG_label_index = cell(size(allsubLabel));

fs = 200;
len_epoch = 510; %in second, for Multiscale sample entropy
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
% n_features = 52; %1~15 Traditional time,16~20 Traditional freq, 21~23 DFA

% Temp function
%extract = cell([n_class,1]);
%% Parse the data to epoch interval
IHR_list = cell(size(allsubLabel));
RRI_list = cell(size(allsubLabel));
%WDFA_list = cell(size(allsubLabel));
%WDFA_curves = cell(size(allsubLabel)); %four combination n and sigma
RRI_res_list = cell(size(allsubLabel));
s = round(len_epoch/len_orig);
for i =1:N_sub
    i
    %PPG
    PPG_data{i} = buffer(allsubPPG{i},fs*len_epoch,fs*(len_epoch-len_orig))';
    PPG_data{i} = PPG_data{i}(s:end,:);
    %Label
    slabel = ceil(s/2);
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
%     WDFA_curves{i} = cell([4,1]);
%     WDFA_curves{i}{1} = WDFA_fun(RRI_res,30,30,1);
%     WDFA_curves{i}{2} = WDFA_fun(RRI_res,30,90,1);
%     WDFA_curves{i}{3} = WDFA_fun(RRI_res,90,30,1);
%     WDFA_curves{i}{4} = WDFA_fun(RRI_res,90,90,1);
%     
%     WDFA_list{i} = cell([4,1]);
%     for j = 1:4
%         WDFA_list{i}{j} = buffer(WDFA_curves{i}{j},len_epoch,(len_epoch-len_orig))';
%         WDFA_list{i}{j} = WDFA_list{i}{j}(s:end,:);
%     end
    
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
    
end

%% MSE
n_features = 20;
for l = 1:N_sub
    disp(l);
    features{l} = zeros([size(RRI_list{l},1),n_features]);
    for ep =  find(PPG_label_index{l} > 0)'
        for scale = 1:10
            tmp = MSE(RRI_list{l}{ep}, 2, 1.0, scale);
            features{l}(ep,scale) = tmp(1,1);
            features{l}(ep,10+scale) = tmp(1,2);
            clear tmp;
        end
    end
    save('features&labels_MSE.mat','features','PPG_label','PPG_label_index');
end
