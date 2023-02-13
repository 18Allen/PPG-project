function [TRANS_rand, EMIS_rand] = Generate_Para_in_rand(Size_of_TRANS, Size_of_EMIS)

TRANS_rand = rand(Size_of_TRANS);
 EMIS_rand = rand(Size_of_EMIS);

TRANS_rand = TRANS_rand ./ sum(TRANS_rand, 2);
 EMIS_rand = EMIS_rand  ./ sum(EMIS_rand, 2);

 