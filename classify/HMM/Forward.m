function Alpha = Forward(Tran, Likeli, ini, ter, O)

Alpha = zeros(size(Tran, 1), size(O, 2) + 1);

Alpha(:, 1) = ini' .* Likeli(:, O(1));

for t = 2 : size(O, 2)
    for s = 1 : size(Tran, 1)
        Alpha(s, t) = (Alpha(:, t-1)' * Tran(:, s)) * Likeli(s, O(t));
    end
end

Alpha(1, size(O, 2) + 1) = Alpha(:, size(O, 2))' * ter;
