function Beta = Backward(Tran, Likeli, ini, ter, O)

Beta = zeros(size(Tran, 1), size(O, 2) + 1);

Beta(:, size(O, 2)) = ter(:, :);

for t = 1 : size(O, 2) - 1
    T = size(O, 2) - t;
    for s = 1 : size(Tran, 1)
        Beta(s, T) = Tran(s, :) * (Likeli(:, O(T+1)) .* Beta(:, T+1));
    end
end

Beta(1, size(O, 2) + 1) = ini(1, :) * (Likeli(:, O(1)) .* Beta(:, 1));
