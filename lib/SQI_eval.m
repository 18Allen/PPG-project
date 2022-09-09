function [TakeorNot, entropy] = SQI_eval(x, len, pad)
    % SQI_eval return the information entropy as an index of signal 
    % quuality. This method is suggested in - Optimal Signal Quality Index for
    % Photoplethysmogram Signals by Mohamed Elgendi.
    % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5597264/pdf/bioengineering-03-00021.pdf
    % ----------------------------
    % Data input:
    % x =  row signal
    % len = length of segment window to calculate entropy 
    % pad = length of overlap (with previous)
    % ----------------------------
    % Output:
    % TakeorNot = 1 to indicate the entropy > .065 (Should be adjust based on your experience)
    % entropy = 'mean' of entropies from all the window segments.
    
    
    threshold = 0.065;
    
    % cut the signal to segment of length 'len' with padding length 'pad'
    % (pad forward)
    tmp  = buffer(x,len,(len-pad))';
    tmp = tmp(len/pad:end,:);
    entropy_list = zeros([size(tmp,1),1]);
    TakeorNot = 1;
    for i = 1:size(tmp,1)
        [counts,edges] = histcounts(tmp(i,:), 'BinWidth', 1e3);
        counts = counts/size(tmp,2);
        entropy_list(i) = median(-sum((counts+ 1e-10).*log(counts+ 1e-10),1)); % prevent Inf
    end
    if min(entropy_list) < threshold
        TakeorNot = 0;
    end
    entropy = min(entropy_list);
end