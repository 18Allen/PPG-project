function [Nu, B] = Viterbi(ini_F, Tran_F, Likeli_F, O, GT)

Nu = zeros(size(Tran_F, 1), size(O, 2) + 1);
 B = zeros(size(Tran_F, 1), size(O, 2) + 1);

Nu(:, 1) = ini_F' .* Likeli_F(:, O(1));

for t = 2 : size(O, 2)
    for s = 1 : size(Tran_F, 1)
        [S, I] = sort(Nu(:, t-1) .* Tran_F(:, s) .* Likeli_F(s, O(t)), 'descend');
        Nu(s, t) = S(1);
        B(s, t) = I(1);
    end
end

[S, I] = sort(Nu(:, size(O, 2)), 'descend');

Nu(1, size(O, 2) + 1) = S(1);
B(1, size(O, 2) + 1) = I(1);
                                                %reduction
Re = zeros(1, size(O, 2));
Re(size(O, 2)) = B(1, size(O, 2) + 1);

for t = 1 : size(O, 2) - 1
    T = size(O, 2) - t;
    Re(T) = B(Re(T+1), T+1);
end
                                                %Correct_rate
A = GT - Re;
i = 0;

for j = 1 : size(GT, 2)
    if A(j) == 0
        i = i + 1;
    end
end

Corre = i/size(GT, 2);


