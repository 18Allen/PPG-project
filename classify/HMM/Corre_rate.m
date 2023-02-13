function Corre = Corre_rate(prediction, real_states)

N = size(real_states, 2);

Err = prediction - real_states;

n = 0;

for ii = 1 : N
    if Err(ii) == 0
        n = n + 1;
    end
end

Corre = n/N;

