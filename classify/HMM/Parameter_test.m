function Corre = Parameter_test(TRANS, EMIS, N)

[seq_test, states_test] = hmmgenerate(10000, TRANS, EMIS);

   seq_training = zeros(N, 10000);
states_training = zeros(N, 10000);

for ii = 1 : N
    [seq, states] = hmmgenerate(10000, TRANS, EMIS);
       seq_training(ii, :) = seq(1, :);
    states_training(ii, :) = states(1, :);
end

ini_ems = zeros(1, length(TRANS));

for ii = 1 : N
    ini_ems(states_training(ii, 1)) = ini_ems(states_training(ii, 1)) + 1;
end

ini_ems = ini_ems / sum(ini_ems);

prediction = [];

for ii = 1 : N
    [TRANS_ems, EMIS_ems] = Generate_Para_in_ems(size(TRANS), size(EMIS), seq_training(ii, :), states_training(ii, :));
    prediction_f = smoothers(TRANS_ems, EMIS_ems, ini_ems, seq_test);
    prediction = [prediction; prediction_f];
end

prediction_voting = voting(prediction, length(TRANS));

Corre = Corre_rate(prediction_voting, states_test);


