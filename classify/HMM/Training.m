function [ini_F, Tran_F, Likeli_F, norm_1, norm_2, O, i] = Training(Tran, Likeli, O)

ini = ones(1, size(Tran, 1))/5;

i = 0;
norm_1 = 1;
norm_2 = 1;

while (norm_1 > 0.00001 || norm_2 > 0.00001) && i < 1000 
        [ini_New, Tran_New, Likeli_New, norm_1, norm_2] = Learning(Tran, Likeli, ini, O);
    
    ini = ini_New;
    Tran = Tran_New;
    Likeli = Likeli_New;
    i = i + 1;
end

ini_F = ini;
Tran_F = Tran;
Likeli_F = Likeli;






