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
% features = cell([n_class,1]);
% n_features = 52; %1~15 Traditional time,16~20 Traditional freq, 21~23 DFA

% Temp function
%extract = cell([n_class,1]);
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
    %PPG
    PPG_data{i} = buffer(allsubPPG{i},fs*len_epoch,fs*(len_epoch-len_orig)).';
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
    
    % Get WDFA curve from RRI_res_data
    RRI_res = RRI_res_data{i};
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

%% MMA
n_features = 1;
features = cell(N_sub,1);
for l = 1:N_sub
    disp(l);
    h_ref = MMA(RRI_data{l}.');
    features{l} = zeros([size(allsubLabel{l},1), n_features]);
    tic
    for ep = find(PPG_label_index{l} > 0)'
        hqs = MMA(RRI_list{l}{ep-slabel+1}.');
        h2qs(:,3) = hqs(:,3) + (mean(h_ref(:,3)) - mean(hqs(:,3)));
        d = sqrt(mean((hqs(:,3)-h2qs(:,3)).^2))/mean(hqs(:,3));
        % features{l}(ep,:) = hqs(:,3).';
        features{l}(ep,:) = d;
        clear hqs;
    end
    save('./features&labels_MMA.mat','features','PPG_label','PPG_label_index');
    toc
end
