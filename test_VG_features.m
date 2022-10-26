%% Workspace Hygiene
clear all;
%% Load data, tools
% mkdir './txt_data';
% for i = 1:30
%     writematrix(RRI_data{i}, ['sub_',num2str(i),'.txt'], 'Delimiter', 'space');
% end
%%
load('./data/allsubPPG_fs200_test.mat');
load('./data/allsubLabel_test.mat');
load('./processed_data.mat')
%addpath('../SST_TF_analysis/TF_anaylsis');
addpath('./lib')
%% Data parsing--------------------------------------
%% Set parameters
N_sub = size(allsubLabel,1);
PPG_data = cell(size(allsubLabel));
PPG_label = cell(size(allsubLabel));
PPG_label_index = cell(size(allsubLabel));

fs = 200;
len_epoch = 270; %in second, 4.5 min
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
n_features = 6; 
%1~44 Traditional time,45~49 Traditional freq, 50~52 DFA
%53 Higuchi fractal dimension, 54 PDFA, 55 ApEn, 56~59 WDFA

% Temp function
%extract = cell([n_class,1]);
%% Parse the data to epoch interval
IHR_list = cell(size(allsubLabel));
RRI_list = cell(size(allsubLabel));
s = round(len_epoch/len_orig/2);
for i =1:N_sub
    i
    %PPG
    PPG_data{i} = buffer(allsubPPG{i},fs*len_epoch,fs*(len_epoch-len_orig))';
    PPG_data{i} = PPG_data{i}(s:end-s+1,:);
    %Label
    PPG_label{i} = allsubLabel{i}(s:end-s+1);
    %Index (Unused = 0, Used = 1)
    index = ones(size(allsubLabel{i}));
    index([1:s-1 end-s+2:end]) = 0;
    PPG_label_index{i} = index;
    %IHR
    IHR_list{i} = IHR_data{i}(2:end);%Because the original one will have 1 extra second
    IHR_list{i} = buffer(IHR_list{i},len_epoch*4,(len_epoch-len_orig)*4)';
    IHR_list{i} = IHR_list{i}(s:end,:);
    %TBD temp. In correspond to IHR(end-12:end) =0 in data_parsing.m
    IHR_list{i}(1,1:12) = 1;IHR_list{i}(end,end-12:end) = 1;
    
    
    % RRI
    RRI_list{i} = cell(size(PPG_label{i}));
    for j = 1:length(PPG_label{i})
        idx = find(locs_data{i} > (j-1)*len_orig*upsampling_rate & locs_data{i} <= (len_epoch+(j-1)*len_orig)*upsampling_rate);
        RRI_list{i}{j} = RRI_data{i}(idx(1:end-1));
        
        % Apply change if less than 300 RRI in an epoch
%         if length(RRI_list{i}{j}) < 300
%             RRI_list{i}{j} = interp1(1:length(RRI_list{i}{j}),RRI_list{i}{j},...
%                 [linspace(1,(300-length(RRI_list{i}{j}))+1,(300-length(RRI_list{i}{j}))*2+1 ),(300-length(RRI_list{i}{j}))+2:length(RRI_list{i}{j})]);
%         end 
    end
    
end

%%  Analysis
for l= 1:N_sub
    features{l} = zeros([size(PPG_data{l},1),n_features]);
    
    for m = 1:size(PPG_data{l},1)
        m
        RRI = RRI_list{l}{m};
        min_Rslope = -inf; max_Lslope = inf;
        % build connection
        connection_table = zeros([length(RRI),length(RRI)]);
        for i = 1:length(RRI)

        end
        
        % Calculate degree
    
    end
    save('features&labels.mat','features','PPG_label','PPG_label_index');

end
clear WDFAs