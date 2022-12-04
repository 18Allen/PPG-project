%% Workspace Hygiene
clear;
%% Load data, tools
load('./data/ECG/allsubECG_fs200.mat');
load('./data/ECG/allsubLabel.mat');

%load('./data/ECG/allsubECG_trouble.mat');
%load('./data/allsubLabel_trouble.mat');
%%
addpath('./lib')
%% Data parsing--------------------------------------
%% Set parameters
N_sub = size(allsubLabel,1);
ECG_data = cell(size(allsubLabel));
ECG_label = cell(size(allsubLabel));
fs = 200; 
len_epoch = 30; %in second
%% Parse the data to epoch interval
for i =1:N_sub
    ECG_data{i} = allsubECG{i}; % Whole night
    ECG_label{i} = allsubLabel{i};
end

%% Get IF-------------------------------
%% Parameter setting
% basic parameters for STFT 
sampling_rate = 200;

% Loop parameters
n_class = 5;
stats = cell([n_class,3]);

% HRV
upsampling_rate = 500; % For ECG_peak_detection
% Generated feature
orig_locs_data = cell([N_sub,1]);
locs_data = cell([N_sub,1]);
IHR_data = cell([N_sub,1]);
RRI_data = cell([N_sub,1]);
RRI_res_data = cell([N_sub,1]);
adapt_RRfreq = cell([N_sub,1]);
% HIAV_data = cell([N_sub,1]);
% HIAV_locs_data = cell([N_sub,1]);

%%  TF analysis
for l= 1:N_sub
    l
    %% RR analysis----------------------------------------
    ECG = ECG_data{l};
    [b, a] = butter(2, [0.5, 25] / (sampling_rate / 2));
    filtered = filtfilt(b, a, ECG);
    pleth = resample(filtered,upsampling_rate,sampling_rate);
    %only use trough to determine cycle. 
    %pleth_p = (pleth -abs(pleth))/2;
    orig_locs = RpeakUltraLong(pleth,upsampling_rate);
    % 
    locs = orig_locs(pleth(orig_locs) < prctile(pleth(orig_locs),97));
    [aRRI,locs,~,~] = RRI_adjust_new(diff(orig_locs)/upsampling_rate,0.4,2,orig_locs(1)/upsampling_rate,upsampling_rate);
    RRI = diff(locs)./upsampling_rate;

     %% (respiratory-induced amplitude variation) RIAV as respiratory
%     [upper,~,~,~] = PPG_peakdetection2(ECG,sampling_rate);
%     [lower,~,~,~] = PPG_peakdetection2(-ECG,sampling_rate);
%     HIAV = zeros(1,length(lower));
%     for i = 1:length(HIAV)
%         a = find(upper >= lower(i));
%         if isempty(a); break; end
%         HIAV(i) = ECG(upper(a(1)))-ECG(lower(i));
%     end
%     HIAV_locs = lower(1:i-1);
%     HIAV = HIAV(1:i-1);

    %% Use entropy to classify bad signal segment of len_orig long
    len_orig = 5;
    test_ECG = buffer(allsubECG{l},len_orig*sampling_rate)';
    entropy_list = zeros([1,size(test_ECG,1)]);
    threshold_entropy = 0.0; % around 4%, but can only work for signal absent cases, other are non-significant
    
    % Mark loc with bad 
    for j = 1:size(test_ECG,1)
       [~,entropy_list(j)] = SQI_eval(test_ECG(j,:),sampling_rate*len_orig,sampling_rate*len_orig);
    end

    %%
    % Make idx_bad_RRI
    bad_list = find(entropy_list < threshold_entropy);
    idx_bad_cycles = zeros(size(locs_data{l}));
    
    % Set RRI in those segment to 0 for imputation
    test_RRI = RRI;
    for j = bad_list
        idx = locs > (j-1)*len_orig*upsampling_rate & locs <= (j)*len_orig*upsampling_rate;
        test_RRI(idx(1:end-1)|idx(2:end)) = 0;
    end
    % Plot bad segments
%     figure
%     plot(linspace(0,length(allsubECG{l})/sampling_rate,length(allsubECG{l})),allsubECG{l});
%     hold on
%     for j =1:size(test_ECG,1)
%     if ismember(j,bad_list)
%     plot([(j-1)*len_orig*sampling_rate+1:j*len_orig*sampling_rate]/sampling_rate,test_ECG(j,:),'r');
%     end
%     end
%     xlabel('time(s)')
%     ylabel('ECG')
%     title(strcat('bad segments(10s) of subject ', num2str(l),' ECG'))
%     legend('ECG','bad seg')

    %%
    % Imputation
    len_pad = 240; %In index, which is around +- 5min
    % padding
    %test_RRI = [test_RRI(1)*ones([1,len_pad]),test_RRI,test_RRI(end)*ones([1,len_pad])];
    test_RRI = [test_RRI(1:len_pad),test_RRI,test_RRI(end-len_pad+1:end)];
    Imp_RRI = RRI;
    for j = bad_list
        idx = locs > (j-1)*len_orig*upsampling_rate & locs <= (j)*len_orig*upsampling_rate;
        cal_list = find(idx(1:end-1)|idx(2:end));
        
        % If there's no RRI to impute, skip.
        if length(cal_list)<1
            % For subject 20
            continue;
        end
        s_RRI = cal_list(1); t_RRI = cal_list(end);
        % range: pad |s_RRI~t_RRI| 5 for alpha in imputation and we don't
        % take too much info from the future data
        range_imp = (len_pad+s_RRI)-len_pad:(len_pad+t_RRI)+5;
        % range_imp = (len_pad+s_RRI)-len_pad:(len_pad+t_RRI)+len_pad;
        sig = test_RRI(range_imp);
        sig([(len_pad+1):(len_pad+length(cal_list))]) = 0;
        [sig] = imputation(sig,1);
        Imp_RRI(cal_list) = sig([(len_pad+1):(len_pad+length(cal_list))]);
        test_RRI(len_pad+cal_list) = sig([(len_pad+1):(len_pad+length(cal_list))]);
    end
    % Plot 
%     figure
%     plot(RRI)
%     hold on
%     plot(Imp_RRI)
%     xlabel('index')
%     ylabel('RRI(s)')
%     legend('pre','imp')
%     hold off

    RRI = Imp_RRI;
    clear test_RRI sig cal_list idx Imp_RRI bad_list idx_bad range_imp test_ECG
    clear entropy_list idx_bad_cycles

    %%
    % IHR
    IF = 1./RRI;
    xx = 1:(upsampling_rate/4):length(pleth);
    IHR = spline(locs,[IF,IF(end)],xx);
    IHR(1:12) = 1;IHR(end-12:end) = 1;

    %% RR_res for WDFA
    % WDFA requires to take the sample 
    IF = RRI;
    xx = 1:(upsampling_rate/20):length(pleth);
    RRI_res = spline(locs,[IF,IF(end)],xx);
    RRI_res(1:80) = 1;RRI_res(end-80:end) = 1;
    RRI_res = RRI_res(1:20:end);
    
    %% STFT of whole night IHR
    IHRtmp = IHR - mean(IHR);
    fs_STFT = 4;
    LowFreq = 0.001/fs_STFT;
    HighFreq = 1.0/fs_STFT;
    fr = 0.001;% 
    [h, Dh, ~] = hermf(30*sampling_rate+1, 1, 5);
    [tfr,tfrtic,~,~] = sqSTFTbase(IHRtmp.', LowFreq, HighFreq, fr/fs_STFT, 1, h(1,:)', Dh(1,:)', 0);

    %% Get bounded-adaptive(b_) HF,LF,VLF,HF/LF
    Total = sum(abs(tfr));

    Indx = find(tfrtic*fs_STFT>= 0.15 & tfrtic*fs_STFT< 0.4);
    [~,HF_loc] = max(abs(tfr(Indx,:)));
    adapt_HF = sum(abs(tfr(Indx(1)-1+HF_loc-0.05/fr:Indx(1)-1+HF_loc+0.05/fr,:)))./Total;

    Indx = find(tfrtic*fs_STFT>= 0.04 & tfrtic*fs_STFT< 0.15);
    [~,LF_loc] = max(abs(tfr(Indx,:)));

    for y = 1:length(LF_loc)
        if Indx(1)-1+LF_loc(y) <= 0.05/fr
            LF_loc(y) = 0.05/fr-Indx(1)+2;
        end
    end
    adapt_LF = sum(abs(tfr(Indx(1)-1+LF_loc-0.05/fr:Indx(1)-1+LF_loc+0.05/fr,:)))./Total;
    adapt_VLF = sum(abs(tfr(0.003/fr:Indx(1)-1+LF_loc-0.05/fr,:)))./Total;
    adapt_LF2HF = adapt_LF./adapt_HF;
    adapt_RRfreq{l} = {adapt_HF; adapt_LF; adapt_VLF; adapt_LF2HF};
    clear tfr tfrtic;

    %% location of R peaks
    orig_locs_data{l} = orig_locs;
    locs_data{l} = locs;
    IHR_data{l} = IHR;
    RRI_data{l} = RRI;
    RRI_res_data{l} = RRI_res;
%     HIAV_data{l} = HIAV;
%     HIAV_locs_data{l} = HIAV_locs;

    %% save
    % save('processed_data.mat','orig_locs_data','locs_data','IHR_data','RRI_data','tfr_data','tfrtic_data');
    save('./data/ECG/processed_data.mat','orig_locs_data','locs_data','IHR_data','RRI_data','adapt_RRfreq','RRI_res_data');

end
