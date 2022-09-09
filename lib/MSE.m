function [SampEn] = MSE(seq, m, r, scale)
% seq is the one-row time series
% m is the maximum length of the pattern
% r is the distance threshold
% scale is the coarse graining / partition of the time series

N = length(seq);
ylen = floor(N/scale);
%% Coarse-graining
y = zeros(1,ylen);
for k = 1:ylen
    y(1,k) = sum(seq(1,(k-1)*scale+1:k*scale))/scale;
end
%% standard deviation
s = sum(seq); s2 = sum(seq*seq');
sd = sqrt((s2-s*s/N)/(N-1));
%% Sample entropy
r_new = r*sd;
N_j = ylen-m;
cont = zeros(1,m+1);

for i = 1:N_j
    for l = i+1:N_j
        k = 0;
        while (k<m && abs(y(1,i+k)-y(1,l+k))<=r_new)
            cont(1,k+1) = cont(1,k+1)+1;
            k = k+1;
        end
        if (k==m && abs(y(1,i+m)-y(1,l+m))<=r_new)
            cont(1,m+1) = cont(1,m+1)+1;
        end
    end
end

SampEn = zeros(1,m);
for mm = 1:m
    if cont(1,mm)==0 || cont(1,mm+1)==0
        SampEn(1,mm) = -log(1/(N_j*(N_j-1)));
    else
        SampEn(1,mm) = -log(cont(1,mm+1)/cont(1,mm));
    end
end


end