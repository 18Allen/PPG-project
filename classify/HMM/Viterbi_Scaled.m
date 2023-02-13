function Re = Viterbi_Scaled(TRANS, EMIS, ini, seq, STATES)

TRANS = TRANS + 10^(-10);
EMIS = EMIS + 10^(-10);
ini = ini + 10^(-10);

Nu = zeros(size(TRANS, 1), size(seq, 2) + 1);
 B = zeros(size(TRANS, 1), size(seq, 2) + 1);

TRANS_log = log(TRANS);
 EMIS_log = log(EMIS);
  ini_log = log(ini);

Nu(:, 1) = ini_log' + EMIS_log(:, seq(1));

for t = 2 : size(seq, 2)
    for s = 1 : size(TRANS, 1)
        [S, I] = sort(Nu(:, t-1) + TRANS_log(:, s) + EMIS_log(s, seq(t)), 'descend');
        Nu(s, t) = S(1);
         B(s, t) = I(1);
    end
end

[S, I] = sort(Nu(:, size(seq, 2)), 'descend');

Nu(1, size(seq, 2) + 1) = S(1);
 B(1, size(seq, 2) + 1) = I(1);
 
Re = zeros(1, size(seq, 2));
Re(size(seq, 2)) = B(1, size(seq, 2) + 1);

for t = -size(seq, 2) + 1 : -1
    Re(-t) = B(Re(-t+1), -t+1);
end
                                             %Correct_rate
A = STATES - Re;
i = 0;

for j = 1 : size(STATES, 2)
    if A(j) == 0
        i = i + 1;
    end
end

Corre = i/size(STATES, 2)

