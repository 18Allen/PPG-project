function output = getTraditionalHRVtime(RRI, dRRI)
% time domain features
            IHR = 1./RRI;
            dtrIHR = detrend(IHR);
            dtrRRI = detrend(RRI);
            ddtrRRI = diff(dtrRRI);
            
            % 4 features, mean & median of HR and RR (both detrended and absolute)
            % Right now, I just use RR
            mean_RRI = mean(RRI);
            mean_dtrRRI = mean(dtrRRI);
            median_RRI = median(RRI);
            median_dtrRRI = median(dtrRRI);

            % mean IHR
            mean_IHR = mean(IHR);
            mean_dtrIHR = mean(dtrIHR);
            
            % Percentile
            quantiles = [];
            sig = IHR;
            quantiles = [quantiles, quantile(sig,0.05),quantile(sig,0.1),quantile(sig,.25),...
                    median(sig),quantile(sig,.75),quantile(sig,.9),quantile(sig,.95)];
            sig = dtrIHR;
            quantiles = [quantiles, quantile(sig,0.05),quantile(sig,0.1),quantile(sig,.25),...
                    median(sig),quantile(sig,.75),quantile(sig,.9),quantile(sig,.95)];
            sig = RRI;
            quantiles = [quantiles, quantile(sig,0.05),quantile(sig,0.1),quantile(sig,.25),...
                    median(sig),quantile(sig,.75),quantile(sig,.9),quantile(sig,.95)];
            sig = dtrRRI;
            quantiles = [quantiles, quantile(sig,0.05),quantile(sig,0.1),quantile(sig,.25),...
                    median(sig),quantile(sig,.75),quantile(sig,.9),quantile(sig,.95)];
            
                
            % TBD detrend IHR. Do we need that
            %HR99 = quantile(IHR, 0.99) ;
            %HR01 = quantile(IHR, 0.01) ;
            
            
            % RRI
            SDNN = std(RRI) ;
            % RRI range. Not using max-min because I don't want glitches to
            % bother.
            rangeRRI = quantile(RRI,0.99) - quantile(RRI,.01); 
            RMSSD = sqrt( sum(dRRI.^2) ./ length(dRRI) );
            SDSD = std(dRRI);
            % Percent of RRI diff > 50 ms
            pNN50 = length(find(1000*dRRI > 50)) / length(dRRI) ;
            % mean absolute deviation
            MAD_RRI = mad(RRI);
            
            %dtrRRI
            SDNN_dtrRRI = std(dtrRRI) ;
            % RRI range. Not using max-min because I don't want glitches to
            % bother.
            range_dtrRRI = quantile(dtrRRI,0.99) - quantile(dtrRRI,.01); 
            RMSSD_dtrRRI = sqrt( sum(ddtrRRI.^2) ./ length(ddtrRRI));
            SDSD_dtrRRI = std(ddtrRRI);
            % Percent of RRI diff > 50 ms
            pNN50_dtrRRI = length(find(1000*ddtrRRI > 50)) / length(ddtrRRI) ;
            % mean absolute deviation
            MAD_dtrRRI = mad(dtrRRI);
            
%             % Additionals
%             pNN20 = length(find(1000*dRRI > 20)) / length(dRRI) ;
%             IQRNN = quantile(RRI, 0.75) - quantile(RRI, 0.25) ;
%             tmp = hist(RRI, linspace(0,3,256)) ;
%             HTI = max(tmp)/length(RRI) ;
%             skew = skewness(RRI');
%             kurt = kurtosis(RRI');
%             AHRR = quantile(RRI, 0.995) - quantile(RRI, 0.005) ;
%             
%             
%             %% nonlinear method
%             SD1 = SDSD ./ sqrt(2);
%             SD2 = sqrt( 2*SDNN.^2 - SDSD.^2 ./ 2 );
            
            %% (Modified) vectorized for better processing
            output = [mean_IHR, mean_dtrIHR, mean_RRI, mean_dtrRRI, median_RRI, median_dtrRRI,...
                        SDNN, rangeRRI, RMSSD, SDSD, pNN50, MAD_RRI,...
                     SDNN_dtrRRI, range_dtrRRI, RMSSD_dtrRRI,SDSD_dtrRRI, pNN50_dtrRRI, MAD_dtrRRI,...
                     quantiles];
%             output = [SDNN RMSSD SDSD pNN20 pNN50 IQRNN ...
%                 AHRR HTI skew kurt ...
%                 SD1 SD2 HR50 HR99 HR01];