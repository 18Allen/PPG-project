function [LF, HF, LHR, VLF, ULF ] = getTraditionalHRVfreq(RRI) 

            
            % frequency domain
            
            if length(min(RRI)) == 0
                fprintf(' ******************** BIG PROBLEMS!\n') ;
                RRI(find(RRI==0)) = [] ;
            end
            
            
            % RRI at sec
            time = cumsum(RRI) ;
            
            Fs = 2 ;            % Sampling frequency
            
            time_resample = time(1): 1/Fs : time(end) ;
            RRI_resample = interp1(time, RRI, time_resample, 'pchip', 'extrap') ;
            %QT0_resample = interp1(time, QT0, time_resample, 'pchip', 'extrap') ;
            
            if mod(length(RRI_resample),2)
                RRI_resample = RRI_resample(1:end-1);
                %QT0_resample = QT0_resample(1:end-1);
                time_resample = time_resample(1:end-1);
            end
            
            
            % remove linear trend following
            % https://www.ahajournals.org/doi/epub/10.1161/01.CIR.96.5.1557
            tmpX = [ones(1,length(time_resample)); time_resample] ;
            betahat = RRI_resample * tmpX' * inv(tmpX * tmpX') ;
            RRI_resample = RRI_resample - betahat * tmpX ;
            %betahat = QT0_resample * tmpX' * inv(tmpX * tmpX') ;
            %QT0_resample = QT0_resample - betahat * tmpX ;
            
            
            
            % get QT vs RRI cross spectrum
%             [Pxy, xi0] = cpsd(RRI_resample,QT0_resample,[],[],[],2) ;
%             [Pxx, xi0] = cpsd(RRI_resample,RRI_resample,[],[],[],2) ;
%             [Pyy, xi0] = cpsd(QT0_resample,QT0_resample,[],[],[],2) ;
%             
%             ind = find(xi0 < 0.2) ;
%             QTVI4 = sum( abs(Pxy(ind)).^2 ./ Pxx(ind) ./ Pyy(ind) ) * (xi0(2) - xi0(1)) ;
%             clear Pxx; clear Pxy; clear Pyy; clear xi0 ;
            
            
            xi = Fs * [1:length(RRI_resample)/2] / length(RRI_resample) ;
            RRIhat = fft(RRI_resample) ;
            
            P1 = abs(RRIhat(2:end/2+1)).^2 ;
            
            
            TOTAL = trapz(xi, P1);
            
            %HF power, 0.15 to 0.50 Hz
            ind = find(xi >= .15 & xi <= .5) ;
            HF = trapz(xi(ind), P1(ind)) ./ TOTAL ;
            
            %LF power, 0.04 to 0.15 Hz
            ind = find(xi >= .04 & xi <= .15) ;
            LF = trapz(xi(ind), P1(ind)) ./ TOTAL ;
            
            LHR = LF/HF ;
            
            %VLF power, 0.003 to 0.04 Hz
            ind = find(xi >= .003 & xi <= .04) ;
            if length(ind) > 2
                VLF = trapz(xi(ind), P1(ind)) ./ TOTAL ;
            else
                VLF = nan ;
            end
            
            %ULF power, 0.0001 to 0.003 Hz
            ind = find(xi >= .0001 & xi <= .003) ;
            if length(ind) > 2
                ULF = trapz(xi(ind), P1(ind)) ./ TOTAL ;
            else
                ULF = nan ;
            end
            
            clear P1 ; clear Y ; clear xi ;
            
            
            %% nonlinear method
            %SD1 = SDSD ./ sqrt(2);
            %SD2 = sqrt( 2*SDNN.^2 - SDSD.^2 ./ 2 );
end