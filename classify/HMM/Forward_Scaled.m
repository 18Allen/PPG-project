function [Alpha, C] = Forward_Scaled(TRANS, EMIS, ini, seq)

Alpha = zeros(size(TRANS, 1), size(seq, 2));
C = zeros(1, size(seq, 2));

alpha = ini' .* EMIS(:, seq(1));
C(1) = 1/(sum(alpha));

Alpha(:, 1) = alpha * C(1);

for t = 2 : size(seq, 2)
    alpha =  diag(EMIS(:, seq(t))) * TRANS' * Alpha(:, t-1);
    C(t) = 1/(sum(alpha));
    
    Alpha(:, t) = alpha * C(t);
end





