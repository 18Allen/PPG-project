function [val] = CPcoupling(IHR, RIAV, fs)
% fs = 4;
% start = 330*30*fs; %13*30*fs; %327*30*fs;
% IHR1 = IHR(start+1:start+fs*270);
% RIAV1 = flow_filt(start+1:start+fs*270);
% figure; plot((start:start+length(IHR1)-1)./fs, IHR1);
% figure; plot((start:start+length(RIAV1)-1)./fs, RIAV1);

%%

% Assume 270 seconds
% 3 overlapping 4.5 epochs signal
IHRc = buffer(IHR, 4.5*30*fs, 4.5*30*fs/2).';
IHRc = IHRc(2:end,:);

RIAVc = buffer(RIAV, 4.5*30*fs, 4.5*30*fs/2).';
RIAVc = RIAVc(2:end,:);

cross_spec = zeros(size(IHRc));
IHRpower = zeros(size(IHRc));
RIAVpower = zeros(size(RIAVc));

for m = 1:size(cross_spec, 1)
    hr = fft(IHRc(m,:)-mean(IHRc(m,:)));
    res = fft(RIAVc(m,:)-mean(RIAVc(m,:)));
    cross_spec(m,:) = conj(hr).*res;
    IHRpower(m,:) = abs(hr).^2;
    RIAVpower(m,:) = abs(res).^2;
end

CPC = zeros(1,size(IHRc,2));
for n = 1:length(CPC)
    CPC(1,n) = abs(mean(cross_spec(:,n)))^4/(mean(IHRpower(:,n))*mean(RIAVpower(:,n)));
end

% figure; plot((0:length(CPC)-1)*fs/length(CPC), CPC);
clear IHRc RIAVc cross_spec IHRpower RIAVpower;

%% peak-to-avg ratio in power
sec = length(CPC)/fs;
% Band
freqLow = floor(0.1*sec);
freqHigh = ceil(0.5*sec);

max_power = max(CPC(freqLow:freqHigh));
Mm = find(CPC(freqLow:freqHigh)==max_power);
MPA = sum(CPC(freqLow+Mm(1)-1:freqLow+Mm(1)+1));
resp_range = CPC(freqLow:freqHigh);
TRA = sum(resp_range);
val = MPA/TRA;

end