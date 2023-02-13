function [TRANS_F, EMIS_F, ini_F, i] = Training_Scaled(TRANS_Old, EMIS_Old, ini_Old, seq)

i = 0;
L_theta_log_diff = 1;

[TRANS_New, EMIS_New, ini_New, C] = Learning_Scaled(TRANS_Old, EMIS_Old, ini_Old, seq);
 L_theta_log_Old = sum(-log(C));
 
 TRANS_Old = TRANS_New;
  EMIS_Old = EMIS_New;
   ini_Old = ini_New;

while L_theta_log_diff > 0.00001

    [TRANS_New, EMIS_New, ini_New, C] = Learning_Scaled(TRANS_Old, EMIS_Old, ini_Old, seq);
    
    L_theta_log_New = sum(-log(C));
   L_theta_log_diff = L_theta_log_New - L_theta_log_Old;
    
          TRANS_Old = TRANS_New;
           EMIS_Old = EMIS_New;
            ini_Old = ini_New;
    L_theta_log_Old = L_theta_log_New;
    
    i = i + 1;
end

  TRANS_F = TRANS_New;
   EMIS_F = EMIS_New;
    ini_F = ini_New;
    









