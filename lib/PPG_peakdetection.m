function [R] = PPG(ppg, Fs)
% function [R, Q] = PPG(ppg, Fs)
%PPG Find the R peaks from the signal ppg sampled at Fs Hz.

tic
ppg = ppg(:);
ppg(isnan(ppg)) = 0;
switch nargin
    case 1
        error('Specify the sampling rate Fs as the second argument.')
end

% parameters
threshold = 0.02;
QRS_length = fix(0.11 * Fs);
beat_length = fix(0.66 * Fs);

% QRS isolation
[b, a] = butter(2, [0.5, 8] / (Fs / 2));
filtered = filtfilt(b, a, ppg);
filtered(filtered < 0) = 0;
u = filtered.^2;
V = filtfilt(ones(1, beat_length) ./ beat_length, 1, u);
indicator = filtfilt(ones(1, QRS_length) ./ QRS_length, 1, u);

mu = smooth(u, Fs*30, 'loess') ;

try threshold = filtfilt(ones(1, 5 * Fs) ./ (5 * Fs), 1, u) * threshold;
catch
    threshold = threshold .* mu ; %* mean(u);
end
t = indicator > V + threshold;

% QRS detection
[M, start] = regexp(sprintf('%i', [0 t']), '1+', 'match');
M = cellfun(@length, M);
R = nan(length(M), 1);
for i = 1:length(M)
    if M(i) > QRS_length
        [~, R(i)] = max(ppg(start(i) + 1:min(start(i) + M(i) - 1, length(ppg))));
        R(i) = R(i) + start(i);
    end
end
R0 = R(~isnan(R));

R = R0 ; 

idx = [] ;
for jj = length(R):-1:2
   if R(jj)-R(jj-1)< ceil(0.3*Fs)
       idx = [idx jj] ;
       
       if jj>2
           if R(jj) - R(jj-2)< ceil(0.3*Fs) 
                idx = [idx jj-1] ;
           end
       end
   end
end

R(unique(idx)) = [] ;



%% get Q points
% Q = zeros(size(R)) ;
% for jj = 1: length(R)
%     idx = [max(1, R(jj)-0.2*Fs): R(jj)] ;
%     [~, tmp] = min(ppg(idx)) ;
%     Q(jj) = max(1, R(jj) - 0.2*Fs) + 1 + tmp(1) ;
% end

%{
threshold = 0.02;
ppg = -ppg ;
filtered = filtfilt(b, a, ppg);
filtered(filtered < 0) = 0;
u = filtered.^2;
V = filtfilt(ones(1, beat_length) ./ beat_length, 1, u);
indicator = filtfilt(ones(1, QRS_length) ./ QRS_length, 1, u);
try threshold = filtfilt(ones(1, 5 * Fs) ./ (5 * Fs), 1, u) * threshold;
catch
    threshold = threshold * mean(u);
end
t = indicator > V + threshold;

% QRS detection
[M, start] = regexp(sprintf('%i', [0 t']), '1+', 'match');
M = cellfun(@length, M);
Q0 = nan(length(M), 1);
for i = 1:length(M)
    if M(i) > QRS_length
        [~, tmp] = max(ppg(start(i) + 1:start(i) + M(i) - 1));
        Q0(i) = tmp + start(i);
    end
end
Q0 = Q0(~isnan(Q0));





%% correct Q
Q = Q0 ;
[Qidx, Ridx] = MergePeaks(Q0, R, 0.3*100) ;
Q = zeros(size(R)) ;
for jj = 1: length(R)
    if ismember(jj, Ridx)
        tmp = find(Ridx == jj) ;
        Q(jj) = Q0(Qidx(tmp)) ;
    else
        idx = [max(1, R(jj)-0.2*100): R(jj)] ;
        [~, tmp] = min(ppg(idx)) ;
        Q(jj) = R(jj) - 0.2*100 + 1 + tmp(1) ;
    end
end
%}

%disp(['Peak detection completed in ' num2str(toc) ' seconds.'])

% figure; plot(ppg, 'black'); hold on; scatter(R, ppg(R), 'green');
end


