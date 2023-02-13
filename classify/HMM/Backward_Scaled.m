function Beta = Backward_Scaled(TRANS, EMIS, seq, C)

Beta = zeros(size(TRANS, 1), size(seq, 2));

Beta(:, size(seq, 2)) = ones(size(TRANS, 1), 1) * C(size(seq, 2));

for ii = -size(seq, 2) + 1 : -1
    Beta(:, -ii) = C(-ii) * TRANS * diag(EMIS(:, seq(-ii+1))) * Beta(:, -ii+1);
end










