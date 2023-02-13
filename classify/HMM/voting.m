function prediction_voting = voting(prediction, num_of_states)

Voting = zeros(num_of_states, size(prediction, 2));

for ii = 1 : size(prediction, 1) % Number of test_prediction
    for jj = 1 : size(prediction, 2) %prediction of one datapoint of a test
        Voting(prediction(ii, jj), jj) = Voting(prediction(ii, jj), jj) + 1;
    end
end

[~, prediction_voting] = max(Voting);

