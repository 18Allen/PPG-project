function [val] = ppgSQI(cut_signal, sampling_rate)

% resamples to 20 Hz from 100 Hz
fs = 20;
cut_signal = resample(cut_signal, fs, sampling_rate);
cut_signal = detrend(cut_signal);

% signal length in seconds
n = length(cut_signal);
sec = n/fs;

% create 3rd order bandpass butterworth filter
[b,a] = butter(3, [.1 1.5]/(fs/2));

% Apply the filter after throwing away a few beginning and
% ending values to avoid delta-shift error when filtering.
% Throwing away these values should have neglible
% effect on power spectrum.
filt_signal= filtfilt(b, a, cut_signal(4:end-4));

%fast Fourier transform of signal
Fsignal = fft(filt_signal);
power = abs(Fsignal).^2/n; 
 
% f = (0:(n-8))./sec;
% f = reshape(f,n-7,1);
% subplot(2,1,2);
% plot(f(1:300),power(1:300))
% plot(f,power);

%%
% finds max power in respiratory range and its assocaited frequency, Mm
max_power = max(power(.5*sec:1.3*sec));
Mm = find(power(.5*sec:1.5*sec)==max_power);

% calculate maximum peak area (MPA) defined as the sum of the three largest continuous set of values
% around largest value M on respiratory frequency band of .1-.75Hz.
% Mm(1) to take location of first peak if there are more than one.
% MPA = power(.5*sec+Mm(1)-1) + power(.5*sec+Mm(1)-2) + power(.5*sec+Mm(1));
MPA = sum(power(.5*sec+Mm(1)-1:.5*sec+Mm(1)+1));

% calculate total respiratory range (TRA)
% sum of all values within respiratory range of .1Hz-.75Hz
resp_range = power(.1*sec:1.5*sec);
TRA = sum(resp_range);

% FFT-RQI is ratio of MPA to TRA
val = MPA/TRA;
% assign RQI value into RQIs cell array
% QIs(i,:) = QI_val;

% QIseq = QIs(:,1); figure; histogram(QIseq, 200);
% disp(length(find(QIseq<=0.4)));
% QIseq = reshape(QIs.',[],1).';
% figure;
% plot(allsubPPG{10}-mean(allsubPPG{10})); hold on; plot(QIseq.*2e4); hold off;

end
