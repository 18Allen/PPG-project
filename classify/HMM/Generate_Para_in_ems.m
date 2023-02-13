function [TRANS_ems, EMIS_ems] = Generate_Para_in_ems(Size_of_TRANS, Size_of_EMIS, seq_training, states_training)

TRANS_ems = zeros(Size_of_TRANS);
 EMIS_ems = zeros([Size_of_TRANS,Size_of_EMIS]);

for ii = 1 : length(states_training) -1
    TRANS_ems(states_training(ii), states_training(ii+1)) = TRANS_ems(states_training(ii), states_training(ii+1)) + 1;
end    

for jj = 1 : length(states_training) 
    EMIS_ems(states_training(jj), seq_training(jj)) = EMIS_ems(states_training(jj), seq_training(jj)) + 1;
end

TRANS_ems = TRANS_ems ./ sum(TRANS_ems, 2);
 EMIS_ems = EMIS_ems  ./ sum(EMIS_ems, 2);



