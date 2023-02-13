function [TRANS, EMIS] = create_T_E(p_h)

TRANS = rand(5, 5);
TRANS = TRANS ./ sum(TRANS, 2);

EMIS_h = rand(5, 3);
EMIS_h = (EMIS_h * p_h) ./ sum(EMIS_h, 2);

EMIS_l = rand(5, 12);
EMIS_l = (EMIS_l * (1 - p_h)) ./ sum(EMIS_l, 2);

EMIS = zeros(5, 15);

EMIS(1, :) = [EMIS_h(1, :), EMIS_l(1, 1:12)];
EMIS(2, :) = [EMIS_l(2, 1:3) EMIS_h(2, :) EMIS_l(2, 4:12)];
EMIS(3, :) = [EMIS_l(3, 1:6) EMIS_h(3, :) EMIS_l(3, 7:12)];
EMIS(4, :) = [EMIS_l(4, 1:9) EMIS_h(4, :) EMIS_l(4, 10:12)];
EMIS(5, :) = [EMIS_l(5, 1:12) EMIS_h(5, :)];