%% Workspace Hygiene
clear all;
%% Load data, tools
% mkdir './txt_data';
% for i = 1:30
%     writematrix(RRI_data{i}, ['sub_',num2str(i),'.txt'], 'Delimiter', 'space');
% end
%%

load('./data/allsubPPG_fs200.mat');
load('./data/allsubLabel.mat');
load('./processed_data_all.mat')
%addpath('../SST_TF_analysis/TF_anaylsis');
addpath('./lib')
%% Data parsing--------------------------------------
%% Set parameters
N_sub = size(allsubLabel,1);
PPG_data = cell(size(allsubLabel));
PPG_label = cell(size(allsubLabel));
fs = 200;
len_epoch = 300; %in second, 5 min
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
n_features = 52; %1~15 Traditional time,16~20 Traditional freq, 21~23 DFA

% Temp function
%extract = cell([n_class,1]);
%% Parse the data to epoch interval
IHR_list = cell(size(allsubLabel));
RRI_list = cell(size(allsubLabel));
s = round(len_epoch/len_orig);
for i =1:N_sub
    i
    %PPG
    PPG_data{i} = buffer(allsubPPG{i},fs*len_epoch,fs*(len_epoch-len_orig))';
    PPG_data{i} = PPG_data{i}(s:end,:);
    %Label
    PPG_label{i} = allsubLabel{i}(s:end);
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
n_features = 45;
for l= 1:N_sub
    features{l} = zeros([size(PPG_data{l},1),n_features]);
    
    for m = 1:size(PPG_data{l},1)
        m
        RRI = RRI_list{l}{m};
 
    %% Traditional time HRV feature
    features{l}(m,1:44) = getTraditionalHRVtime(RRI,diff(RRI));
%     features{l}(m,45) = HFD(RRI,10);
    features{l}(m,45:49) = getTraditionalHRVfreq(RRI);
    % DFA features
%     if length(RRI) < 300;
%         %fprintf('ERROR!\newline')
%         RRI = [RRI,RRI(end-(300-length(RRI)-1):end)];
%     end
    
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
    features{l}(m,50) = alphaAll;
    features{l}(m,51) = alpha1;
    features{l}(m,52) = alpha2;
    
    
    
    features{l}(m,53) = HFD(RRI,20);
    end
    save('features&labels.mat','features','PPG_label');

end
