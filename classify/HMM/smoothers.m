function prediction = smoothers(TRANS, EMIS, ini, seq)

[Alpha, C] = Forward_Scaled(TRANS, EMIS, ini, seq);
      Beta = Backward_Scaled(TRANS, EMIS, seq, C);

prediction = zeros(1, size(Alpha, 2));

for ii = 1 : size(Alpha, 2)
    Nu = Alpha(:, ii) .* Beta(:, ii);
    [S, I] = sort(Nu, 'descend');
    prediction(ii) = I(1);
end


