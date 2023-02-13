function [TRANS_New, EMIS_New, ini_New] = Learning(L_TRANS, L_EMIS, L_ini, seq)

L_ter = ones(size(L_TRANS, 1), 1);

Alpha = Forward(L_TRANS, L_EMIS, L_ini, L_ter, seq);
Beta = Backward(L_TRANS, L_EMIS, L_ini, L_ter, seq);

%[Alpha, C] = Forward_Scaled(L_TRANS, L_EMIS, L_ini, seq);
%      Beta = Backward_Scaled(L_TRANS, L_EMIS, seq, C);

         T = size(seq, 2);
         N = size(L_TRANS, 2);

        Xi = zeros(T - 1, N, N);
     Gamma = zeros(T, N);

  TRANS_New = zeros(size(L_TRANS));
EMIS_New = zeros(size(L_EMIS));
   ini_New = zeros(size(L_ini));

for t = 1 : T - 1 
    for ii = 1 : N
        for jj = 1 : N
            Xi(t, ii, jj) = Alpha(ii, t) * L_TRANS(ii, jj) * L_EMIS(jj, seq(t+1)) * Beta(jj, t+1) / Alpha(1, T+1);
        end
    end
end

for t = 1 : T
    for jj = 1 : N
        Gamma(t, jj) = Alpha(jj, t) * Beta(jj, t) / Alpha(1, T + 1);
    end
end

ini_New = Gamma(1, :);

for ii = 1 : N
    for jj = 1 : N
        TRANS_New(ii, jj) = sum(Xi(:, ii, jj)) / sum(sum(Xi(:, ii, :)));
    end
end

S = 0;

for ii = 1 : size(L_EMIS, 1)
    for jj = 1 : size(L_EMIS, 2)
        for t = 1 : T
            if seq(t) == jj
                S = S + Gamma(t, ii);
            end
        end
        EMIS_New(ii, jj) = S / sum(Gamma(:, ii));
        S = 0;
    end
end

norm_1 = norm(TRANS_New - L_TRANS);
norm_2 = norm(EMIS_New - L_EMIS);


