function R = HRCFTG(ecg, Fs)
% Find the R peaks in the signal ecg sampled at Fs Hz.
% INPUT
%     ecg    : Electrocardiogram (column vector).
%     Fs     : Sampling rate of the signal in Hz.
% OUTPUT
%     R      : Locations (in samples) of the detected R peaks.
% Modified from "Fast QRS Detection with an Optimized Knowledge-Based 
% Method: Evaluation on 11 Standard ECG Databases" by M. Elgendi, 2013.
% Written by John Malik on 2019.8.1, john.malik@duke.edu.

tt = tic;
switch nargin
    case 1
        error('Specify the sampling rate Fs as the second argument.')
end

% parameters
NN = length(ecg); 
beta = .08;
QRS_length = round(.097 * Fs);
signal_length = min(5 * Fs, NN);

% bandpass filter
[b, a] = butter(3, [8, 20] / (Fs / 2));
y = filtfilt(b, a, ecg);
    
% QRS isolation 
u = y.^2;
z = movmean(u, QRS_length);

% stft parameters
hlength = 5 * Fs + 1 - rem(Fs, 2); % window length, must be odd
Lh = (hlength - 1) / 2; % window width
h = .5 * (1 - cos(2 * pi * (0:hlength-1) / (hlength-1)))'; % hann window
hop = Fs; % estimate heart rate every one second
n = 2 * Fs + 1; % resolution of frequency axis, must be odd
N = 2 * (n - 1); % number of fft points
hf = 6; % 360 bpm
lf = 0.5; % 30 bpm
t = hop:hop:NN; % samples at which to take the FT
tcol = length(t); % number of frequency estimates

% signal to take FT
x = z;

% non-negative frequency axis, cropped
fr = Fs / 2 * linspace(0, 1, n)'; 
eta = fr >= lf & fr <= hf;
fr = fr(eta);

% DFT matrix
% w = 2*pi*1i / N * (0:N-1)';
w = 2*pi*1i / N * find(eta);
D = exp(-w*(0:hlength-1));

% STFT, serial
f = zeros(tcol, 1); % rough frequency estimate
lambda = .01; % penalty weight
rSig = zeros(size(h));
for icol = 1:tcol
    ti = t(icol); 
    tau = -min(Lh, ti - 1):min(Lh, NN - ti);
    rSig(Lh + 1 + tau) = x(ti + tau);
    rSig = rSig - mean(rSig); % - mean(x(ti + tau)); % remove low-frequency content
    tfr = D * (rSig .* h); 
    if icol == 1
        [~, i] = max(abs(tfr.^2));
    else
        tfr = abs(tfr.^2);
        tfr = tfr / sum(tfr);
        [~, i] = max(tfr - lambda * (fr - fr(i)).^2); % penalty for jumping
    end
    f(icol) = fr(i);
end

% time-varying threshold
win = round(Fs * .611 ./ sqrt(f)); % Bazett's formula
[g, ~, id] = unique(win);
id = interp1(t, id, 1:NN, 'nearest', 'extrap');
V = zeros(size(u));
for i = 1:length(g)
    s = id == i;
    v = movmean(u, g(i));
    V(s) = v(s);
end

% alpha
alpha = movmean(u, signal_length) * beta;

% QRS detection
r = z > V + alpha;

% ensure detected segments are long enough
QRS = movsum(r, QRS_length) == QRS_length;

% R peak detection
[c, d] = regexp(sprintf('%i', QRS), '1+');
R = nan(length(c), 1);
for i = 1:length(c)
    [~, R(i)] = max(z(c(i):d(i)));
    R(i) = R(i) + c(i) - 1;
end

disp(['Peak detection completed in ' num2str(toc(tt)) ' seconds.'])

end


