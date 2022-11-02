function [hfd] = HFD(seq, kmax)

N = length(seq);
Lmk = zeros(kmax, kmax);
for k = 1:kmax
    for m = 1:k
        Lmki = 0;
        for i = 1:fix((N-m)/k)
            a = m+i*k;
            b = m+(i-1)*k;
            Lmki = Lmki + abs(seq(1,a)-seq(1,b));
        end
        Ng = (N-1)/(fix((N-m)/k)*k);
        Lmk(m,k) = (Lmki*Ng)/k;
    end
end
Lk = zeros(1,kmax);
for k = 1:kmax
    Lk(1,k) = sum(Lmk(1:k,k))/k;
end

lnLK = zeros(size(Lk));
for i = 1:kmax
    lnLK(1,i) = log(Lk(1,i));
end
lnk = zeros(size(Lk));
for i = 1:kmax
    lnk(1,i) = log(1./i);
end

b = polyfit(lnk, lnLK, 1);
hfd = b(1);

end