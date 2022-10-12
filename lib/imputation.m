function [Imp] = imputation(sig, fs)
% find the calibration segments (assume knowledge 0)
% For convienience, sig == 111 means those that should not be taken
%%
tmp = sig==0 ; % 
[M, start] = regexp(sprintf('%i', tmp'), '1+', 'match') ;
LEN = zeros(size(M)) ;

for ppp = 1: length(LEN) ;
    LEN(ppp) = length(M{ppp}) ;
end


% find sufficiently long missing segments
% and run Taken's embedding
alpha = 5; % length of left&right window
FLen = fs*alpha ; % alpha
Imp = sig ;
avg_val = mean(sig(sig ~= 0)); % for the case where FIT does not exists
%%
if find(LEN > 0.1*fs);
    %fprintf(['\t*** Handle calibration\n']) ;

    QQ = find(LEN > 0.1*fs) ;
    
    %idx_bad_cycles = find(sig == 111);% 111 for being unselectable
    idx_bad_cycles = [find(sig == 111| sig == 0)];
%     for qi = 1:length(QQ)
%         idx_bad_cycles = [idx_bad_cycles, start(QQ(qi)): start(QQ(qi)) + LEN(QQ(qi)) - 1 ];
%     end
    for qi = 1: length(QQ)
        FITexists = 0;
        tmp = QQ(qi) ;

        % correct those zeros in the recording

        cal = start(tmp): start(tmp) + LEN(tmp) - 1 ;
        if start(tmp)-FLen < 1 ;
            llen = start(tmp)-1 ;
        else
            llen = FLen ;
        end

        if start(tmp) + LEN(tmp) + FLen > length(sig) ;
            rlen = length(sig) - (start(tmp)+LEN(tmp)-1) ;
        else
            rlen = FLen ;
        end

        % determine the pattern to fit. Pattern to the l/r of sig(tmp)
        patternL = [start(tmp)-llen: start(tmp)-1] ;
        %patternR = [start(tmp) + LEN(tmp) + 1: start(tmp) + LEN(tmp) + rlen] ;
        patternR = [start(tmp) + LEN(tmp) : start(tmp) + LEN(tmp) + rlen-1] ;
        
        X = sig([patternL patternR]) ;

        Dist = inf ;
        pivot = -1 ;
        for ppp = [1: start(tmp)-LEN(tmp)-llen-rlen start(tmp)+LEN+1: length(sig)-LEN(tmp)-rlen-llen-1]
            patternL = [ppp: ppp+llen-1] ;
            patternR = [ppp+llen + (LEN(tmp)-1) + 1: ppp+llen + (LEN(tmp)-1) + rlen] ;
            
            % Exclude those segment intersecting idx_bad_cycles in patternL,middle, and patternR
            if sum(ismember([ppp: ppp+llen + (LEN(tmp)-1) + rlen],idx_bad_cycles)) > 0
                continue
            end
            
            X0 = sig([patternL patternR]) ;
            
            
            % find the most similar pattern from the remainig signal for
            % imputation
            dd = norm(X0 - X) ;
            if dd < Dist ;            
                Dist = dd ; pivot = ppp ;
                FIT = sig(ppp+llen: ppp+llen+LEN(tmp)-1);
                FITexists = 1;
            end
        end
        if FITexists == 1
            Imp(cal) = FIT ;
        else
            Imp(cal) = ones(size(cal))*avg_val;
        end
    end
end
