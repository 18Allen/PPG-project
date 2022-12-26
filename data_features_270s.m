%% Workspace Hygiene
clear;
%% Load data, tools
datapath = 'C:/Users/scien/Documents/TIDIS project/SST/toNengTai/toNengTai/data/';
load([datapath, '/allsubLabel.mat']);
% load([datapath, '/allsubLabel_normal.mat']);
allsubPPG = {};
% Load PPG
PPG_data_files = string(ls([datapath, '/allsubPPG_*.mat']));
% PPG_data_files = string(ls([datapath, '/allsubPPG_normal.mat']));
for i = 1:size(PPG_data_files,1)
    load([datapath, PPG_data_files{i}]);
    allsubPPG = [allsubPPG; Channels];
end

% For ECG
% load('./data/ECG/allsubECG_fs200.mat');
% load('./data/ECG/allsubLabel.mat');
% load('./data/ECG/processed_data.mat');

% load preprocessed RRI
load('./data/processed_data.mat');
load('./data/processed_data_FLOW.mat');

% addpath('../SST_TF_analysis/TF_anaylsis');
addpath('./lib')
%% Data parsing--------------------------------------
%% Set parameters
N_sub = size(allsubLabel,1);
% PPG_data = cell(size(allsubLabel));
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
features = cell([N_sub,1]);
n_features = 74;
% 1~46 Traditional time, 47~56 Traditional freq (2 HFpole features left),
% 57 ApEn, 58 Higuchi fractal dimension, 59~60 teager energy,
% 61~67 Visibility graph, 68~72 Arousal probability,
% 73 cardiopulmanory coupling, 74 PSD quality.
% DFA, PDFA and WDFA see other programs.

% Temp function
% extract = cell([n_class,1]);
%% Parse the data to epoch interval
% IHR_list = cell(size(allsubLabel));
% RRI_list = cell(size(allsubLabel));
% WDFA_list = cell(size(allsubLabel));
% WDFA_curves = cell(size(allsubLabel)); %four combination n and sigma
% RRI_res_list = cell(size(allsubLabel));
% FLOW_list = cell(size(allsubLabel));
s = round(len_epoch/len_orig);
fs2 = 4; % For frequency features
for l = 1:N_sub
    l
    slabel = ceil(s/2);
    % PPG
    PPG_data = buffer(allsubPPG{l},fs*len_epoch,fs*(len_epoch-len_orig))';
    % PPG_data{i} = buffer(allsubECG{i},fs*len_epoch,fs*(len_epoch-len_orig))';
    PPG_data = PPG_data(s:end,:);

    % Label
    PPG_label{l} = allsubLabel{l};
    %Index (Unused = 0, Used = 1)
    index = ones(size(allsubLabel{l}));
    index([1:slabel-1 end-slabel+2:end]) = 0;
    PPG_label_index{l} = index;
     
    % IHR
    IHR_list = buffer(IHR_data{l},len_epoch*4,(len_epoch-len_orig)*4)';
    IHR_list = IHR_list(s:end,:);
    %TBD temp. In correspond to IHR(end-12:end) =0 in data_parsing.m
    %IHR_list{i}(1,1:12) = 1;IHR_list{i}(end,end-12:end) = 1;

    % FLOW
    FLOW_list = buffer(FLOW{l},len_epoch*4,(len_epoch-len_orig)*4)';
    FLOW_list = FLOW_list(s:end,:);
    
    % Get WDFA curve from RRI_res_data
%     RRI_res = RRI_res_data{i};
%     WDFA_curves{i} = cell([4,1]);
%     WDFA_curves{i}{1} = WDFA_fun(RRI_res,30,30,1);
%     WDFA_curves{i}{2} = WDFA_fun(RRI_res,30,90,1);
%     WDFA_curves{i}{3} = WDFA_fun(RRI_res,90,30,1);
%     WDFA_curves{i}{4} = WDFA_fun(RRI_res,90,90,1);
%     
%     WDFA_list{i} = cell([4,1]);
%     for j = 1:4
%         WDFA_list{i}{j} = buffer(WDFA_curves{i}{j},len_epoch,(len_epoch-len_orig))';
%         WDFA_list{i}{j} = WDFA_list{i}{j}(slabel:end,:);
%     end
    
    % RRI_res_list for teager energy
    RRI_res_list = buffer(RRI_res_data{l},len_epoch,(len_epoch-len_orig))';
    RRI_res_list = RRI_res_list(s:end,:);
    
    % RRI
    RRI_list = cell(size(PPG_label{l}));
    for j = 1:length(PPG_label{l})
        idx = find(locs_data{l} > (-(len_epoch-len_orig)/2 + (j-1)*len_orig)*upsampling_rate & locs_data{l} <= ((len_epoch-len_orig)/2+(j)*len_orig)*upsampling_rate);
        RRI_list{j} = RRI_data{l}(idx(1:end-1));     
        % Apply change if less than 300 RRI in an epoch
%         if length(RRI_list{i}{j}) < 300
%             RRI_list{i}{j} = interp1(1:length(RRI_list{i}{j}),RRI_list{i}{j},...
%                 [linspace(1,(300-length(RRI_list{i}{j}))+1,(300-length(RRI_list{i}{j}))*2+1 ),(300-length(RRI_list{i}{j}))+2:length(RRI_list{i}{j})]);
%         end 
    end
    RRI_list = RRI_list(slabel:end-slabel+1); % notice: start from slabel

    %% Analysis
    features{l} = zeros([size(allsubLabel{l},1), n_features]);

    freqFeature = [];
    for f = 1:4  %f stands for HF, LF, VLF and LF2HF
        tmp = buffer(adapt_RRfreq{l}{f}, fs2*len_epoch, fs2*(len_epoch-len_orig));
        tmp = tmp(:, s:end);
        freqFeature = [median(tmp); freqFeature];
    end

    for m = find(PPG_label_index{l} > 0).'
        m
        RRI = RRI_list{m-slabel+1};
        IHR = IHR_list(m-slabel+1,:);
        RRI_res = RRI_res_list(m-slabel+1,:);
        flow_filt = FLOW_list(m-slabel+1,:);
        PPG = PPG_data(m-slabel+1,:);

        %% Quality
        PPG_cut = buffer(PPG, 10*sampling_rate).';
        qual = zeros(len_epoch/10, 1);
        for c = 1:len_epoch/10
            qual(c) = ppgSQI(PPG_cut(c,:), sampling_rate);
            % disp(qual(c));
        end
        if sum(qual>0.5)<21
            continue;
        end

        %% Traditional frequency domain feature (10)
        len = length(IHR);
        IHR = IHR - mean(IHR);
        specIHR = abs(fft(IHR)/len);
        PS = specIHR(1:len/2+1); PS(2:len/2) = 2*specIHR(2:len/2);
        VLF = 0.04; LF = 0.15; HF = 0.4;
        HFpower = sum(PS(ceil(LF*len/fs2+1):floor(HF*len/fs2+1)));
        LFpower = sum(PS(ceil(VLF*len/fs2+1):floor(LF*len/fs2+1)));
        VLFpower = sum(PS(1:floor(VLF*len/fs2+1)));
        LF2HFratio = LFpower/HFpower;

        features{l}(m,47:50) = freqFeature(:, m-slabel+1)';
        features{l}(m,51:54) = [HFpower LFpower VLFpower LF2HFratio];
        % mean repository freq and power
        [power, ff] = max(specIHR); ff = ff/len*fs2;
        features{l}(m,55) = ff;
        features{l}(m,56) = power;
                
        %% Traditional time HRV feature (46)
        features{l}(m,1:46) = getTraditionalHRVtime(RRI,diff(RRI));
    
        %% DFA features (3)
%         % select scale for the whole DFA
%         Q = exp(linspace(log(10),log(300),37)); 
%         pts = round(Q);
%         order = 2; % In the paper they use quadratic one
%         F = DFA_fun(RRI,pts,order);
%     
%         % A: a 2x1 vector. A(1) is the scaling coefficient "alpha",
%         % A(2) the intercept of the log-log regression, useful for plotting (see examples).
%         A = polyfit(log(pts),log(F)',1);
%         alphaAll = A(1) ;   
%     
%         idx = find(Q>=10 & Q<=40) ;
%         A = polyfit(log(pts(idx)),log(F(idx))',1);
%         alpha1 = A(1) ;
%     
%         idx = find(Q>=70 & Q<=300) ;
%         A = polyfit(log(pts(idx)),log(F(idx))',1);
%         alpha2 = A(1) ;
%         features{l}(m,55) = alphaAll;
%         features{l}(m,56) = alpha1;
%         features{l}(m,57) = alpha2;
    
        %% Higuchi fractal dimension (1)
        features{l}(m,57) = HFD(RRI,20);
        
        %% PDFA
%         P = PDFA(RRI(end-63:end),64,1);
%         p =  polyfit(linspace(0,1,64), P, 1);
%         features{l}(m,59) = p(1);
        
        %% ApEn of binary RRI diff in a interval of TBD (right now 5 min)
        features{l}(m,58) = ApEn(5, 0.2*std(diff(RRI)>0), diff(RRI)>0, 1);
        % rri = diff(RRI)>0; (TBD)

        %% WDFA feature
%         WDFAs = cell([4,1]);
%         for j = 1:4
%         WDFAs{j} = WDFA_list{l}{j}(m,:);
%         features{l}(m,55+j) = max(WDFAs{j});
%         end

        %% Teager energy (2)
        te = teager_energy_func(RRI_res);
        features{l}(m,59) = mean(te);
        features{l}(m,60) = std(te);

        %% Visibility graph (7)
        % [deg_mean, deg_std, ast, per_low, per_high, cc_mean, cc_std] = VG_func(RRI);
        features{l}(m,61:67) = VG_func(RRI);

        %% Arousal probability (5)
        if ~isempty(find(isnan(RRI)))
            features{l}(m,68:72) = nan;
            continue;
        end
        % Arousal probability of normalized RRI
        ap = arousal_prob_func(RRI/mean(RRI));
        if isempty(ap)
            features{l}(m,68:72) = nan;
            continue;
        end
        features{l}(m,68:72) = [max(ap), mean(ap), median(ap), min(ap), std(ap)];

        %% Cardiopulmanory coupling
        features{l}(m,73) = CPcoupling(IHR, flow_filt, 4);

        %% Quality as a feature
        features{l}(m,74) = mean(qual);
        % = SQI_eval(PPG_data{l}(m,:), length(PPG_data{l}(m,:)), length(PPG_data{l}(m,:)));

    end
    % save('./features&labels.mat','features','PPG_label','PPG_label_index');
    % save('./data/ECG/features&labels.mat','features','PPG_label','PPG_label_index');
end
clear WDFAs