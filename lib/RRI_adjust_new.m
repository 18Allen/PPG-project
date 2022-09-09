function [RRI_output, R_output, added, removed] = RRI_adjust_new(RRI, low, high,firstRloc, Hz)
%RRI: RR interval (in sec) directly from R peaks
%low: low cut off thresold 
%high: high cut off thresold
%firstRloc: the first R peak (in sec)
%Hz: sampling rate
%
%output: adjusted RRI & R peaks; number of beats added/removed

if size(RRI, 1) ~= 1 ; RRI = RRI'; end

%% first remove R peaks with RRI <= low
num0 = length(RRI);
if sum(RRI(1:end-1) <= low)     
    Rbeats = [firstRloc firstRloc+cumsum(RRI)];
    bad_ind = find(RRI(1:end-1) <= low);
    
     % merge small rri to latter (what?)
    Rbeats(bad_ind+1) = []; 
    RRI = diff(Rbeats); 

end

num1 = length(RRI);
minus1 = num1-num0;


%% next adjust RRI >= high

bad_ind = find(RRI >= high);
good_ind = find(RRI > low & RRI < high);
good_rri = RRI(RRI > low & RRI < high);




if ~isempty(bad_ind)
    temp = num2cell(RRI);

    
    for ind = bad_ind
        bad_rri = RRI(ind);
        ii = 1;
        mask = good_ind > ind - 20*ii  & good_ind < ind + 20*ii;
        while sum(mask) == 0  % need larger mask
            ii = ii+1;
            mask = good_ind > ind - 20*ii & good_ind < ind + 20*ii;
        end

        temp_rri = median(good_rri(mask));
        n = round(bad_rri/temp_rri);
        if n == 1
            n = 2; % must cut if too long.
        end 
        % break large rri into smaller pieces, yet remain the total length
        % for R
        temp{ind} = (bad_rri/n)*ones(1,n);
        
    end   

    RRI = cell2mat(temp);

end

num2 = length(RRI);
plus1 = num2 - num1;


%% Second round -- adaptive adjustment

% this is the parameter chosen by the user to adjust HR discrepancy adaptively 
gamma = 0.7 ; 

%% remove small spikes from below
L = length(RRI);
bad_ind = [];
for k = 1:L-1
    startpt = max(1, k-15);
    endpt = min(L, k+15);
    temp = median([RRI(startpt:k-1) RRI(k+1:endpt)]); % havent used this
    low1 = quantile([RRI(startpt:k-1) RRI(k+1:endpt)], .25) * gamma ; 
    if RRI(k) < low1 || RRI(k) < low
        bad_ind = [bad_ind k];
    end
end

% R peak locations
Rbeats = [firstRloc firstRloc+cumsum(RRI)];
    
if length(bad_ind)
    Rbeats(bad_ind+1) = [];    
    RRI = diff(Rbeats) ; 
end
 

num3 = length(RRI);
minus2 = num3 - num2;


% small debug
Radj = [firstRloc firstRloc+cumsum(RRI)];


%% remove spikes from above
L = length(RRI);
bad_ind = [];

for k = 1:L
    startpt = max(1, k-15);
    endpt = min(L, k+15);
    high1 = quantile([RRI(startpt:k-1) RRI(k+1:endpt)], .75) * (1/gamma);
    if RRI(k) >= high1 || RRI(k) >= high
        bad_ind = [bad_ind k] ;
    end
end    

if ~isempty(bad_ind)
    good_ind = 1:L;
    good_ind(bad_ind) = [];
    good_rri = RRI(good_ind);
    temp = num2cell(RRI);


    
    for ind = bad_ind
        bad_rri = RRI(ind);
        ii = 1;
        mask = good_ind > ind - 15*ii  & good_ind < ind + 15*ii;
        while sum(mask) == 0  % need larger mask
            ii = ii+1;
            mask = good_ind > ind - 15*ii & good_ind < ind + 15*ii;
        end
        temp_rri = median(good_rri(mask));
        n = round(bad_rri/temp_rri);
        if n == 1
            n = 2;
        end 
        % break large rri into smaller pieces, yet remain the total length
        % for R
        temp{ind} = (bad_rri/n)*ones(1,n);
        
    end   
    RRI = cell2mat(temp);

end

num4 = length(RRI);
plus2 = num4 - num3;
%% summarize everythign for output
added = plus1 + plus2;
removed = minus1 + minus2;

RRI_output = RRI;

R_output = round([firstRloc firstRloc + cumsum(RRI_output)]*Hz) ;

end