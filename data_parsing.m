%% Workspace Hygiene
clear all;
%% Load data, tools
load('./data/allsubPPG_fs200_test.mat');
load('./data/allsubLabel_test.mat');
%addpath('../SST_TF_analysis/TF_anaylsis');
addpath('./lib')
%% Data parsing--------------------------------------
%% Set parameters
N_sub = size(allsubLabel,1);
PPG_data = cell(size(allsubLabel));
PPG_label = cell(size(allsubLabel));
fs = 200;
len_epoch = 30; %in second
%% Parse the data to epoch interval
for i =1:N_sub
    PPG_data{i} = buffer(allsubPPG{i},fs*len_epoch)';
    PPG_label{i} = allsubLabel{i};
end
%% Get IF-------------------------------
%% Parameter setting
% basic parameters for STFT 
sampling_rate = 200;
basicTF.win = sampling_rate*10+1;%window length(sec)(need to change!) % sample in a minute
basicTF.hop = 10;
basicTF.fs = sampling_rate ;
basicTF.fr = 0.01; % frequency resolution
basicTF.feat = 'SST11'; % two option: STFT or SST11 (one window rejection)

% advanced parameters for STFT 
advTF.num_tap = 1; % Number of tap in ConceFT
advTF.win_type = 'Gauss'; % Only 2-tap basis here
advTF.Smo = 1; % alpha smoothing; 1 = no smoothing
advTF.Rej = 0; % The bandwidth of window rejection;
advTF.ths = 1E-9; % Global threshold of STFT
advTF.HighFreq = 1.5/basicTF.fs; % highest frequency/sampling freq(need to change!) %highest I want to see
advTF.LowFreq = .5/basicTF.fs; % lowest frequency/sampling freq(need to change!)
advTF.lpc = 0;

% %% this part is used to segment input signal when it is too long(8 hours for example)
seglength = 10;%sec(need to change!)
pad = 5;%(need to change!) % 前後的pad 100s
f_component_size = round((advTF.HighFreq - advTF.LowFreq)/(basicTF.fr)*basicTF.fs);
if mod(sampling_rate,basicTF.hop*seglength) ~= 0
    error('#hops times seglength should divide the sampling rate!')
end

% This part is for curve extraction to avoid edge effect
s = 3.5;  % graph start time
t = 26.5; % graph end time
interval = basicTF.fs/ basicTF.hop;


% Loop parameters
n_class = 5;
stats = cell([n_class,3]);

% HRV
upsampling_rate = 500; % For PPG_peak_detection

% Generated feature
features = cell([n_class,1]);
n_features = 18;

% Temp function
extract = cell([n_class,1]);

%%  TF analysis
for l= 3%1:N_sub
    

    features{1} = zeros([size(PPG_data{1},1),n_features]);
    extract{l} =  zeros([size(PPG_data{1},1),upsampling_rate*len_epoch]);
    for m = 22%1:size(PPG_data{l},1)
        m
        
        PPG = PPG_data{l}(m,:);
       
    %% RR analysis----------------------------------------
        [b, a] = butter(2, [0.5, 15] / (sampling_rate / 2));
        filtered = filtfilt(b, a, PPG);
        pleth = resample(filtered,upsampling_rate,sampling_rate);
        %dPPG = derivative(1:15000,up_PPG);
        %[locs, ~, ~,~] = PPG_peakdetection2(dPPG,upsampling_rate);
        [~,locs,~,~] = PPG_peakdetection2(pleth,upsampling_rate);
        [aRRI,locs,~,~] = RRI_adjust_new(diff(locs)/upsampling_rate,0.4,2,locs(1)/upsampling_rate,upsampling_rate);
        RRI = diff(locs)./upsampling_rate;
        IF = 1./RRI;
        xx = 0:(upsampling_rate/4):length(pleth);
        IHR = spline(locs,[IF,IF(end)],xx);
        IHR(1:12) = 0;IHR(end-12:end) = 0;
    
        
     %% Plot R-peak location
%     pleth = up_PPG;
%     figure;
%     plot(pleth);
%     hold on;
%     scatter(locs,pleth(locs),'r.'); % similar to cumsum(floor(1000./RRI))+locs(1)
%     xlim([0 locs(end)])
%     %plot(dPPG)
%     %scatter(tmp_locs,pleth(tmp_locs),'g*');
%     %scatter(locs,dPPG(locs));
%     legend('','newlocs')
%     hold off

    %% Traditional time HRV feature
%     features{l}(m,1:15) = getTraditionalHRVtime(RRI,diff(RRI));
    %% DFA features
%     if length(RRI) < 300;
%         %fprintf('ERROR!\newline')
%         RRI = [RRI,RRI(end-(300-length(RRI)-1):end)];
%     end
%     
%     % select scale for the whole DFA
%     Q = exp(linspace(log(10),log(300),37)) ; 
%     %pts = round(exp(Q)) ;
%     pts = round(Q);
%     order = 2; % In the paper they use quadratic one
%     F = DFA_fun(RRI,pts,order);
% 
%     % A: a 2x1 vector. A(1) is the scaling coefficient "alpha",
%     % A(2) the intercept of the log-log regression, useful for plotting (see examples).
%     A = polyfit(log(pts),log(F)',1);
%     alphaAll = A(1) ;   
% 
%     idx = find(Q>=10 & Q<=40) ;
%     A = polyfit(log(pts(idx)),log(F(idx))',1);
%     alpha1 = A(1) ;
% 
%     idx = find(Q>=70 & Q<=300) ;
%     A = polyfit(log(pts(idx)),log(F(idx))',1);
%     alpha2 = A(1) ;

    end
    %save('features.mat','features','extract');

end
