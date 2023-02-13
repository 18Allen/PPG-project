function [TRANS_New, EMIS_New, ini_New, C] = Learning_Scaled(TRANS_Old, EMIS_Old, ini_Old, seq)

[Alpha, C] = Forward_Scaled(TRANS_Old, EMIS_Old, ini_Old, seq);
      Beta = Backward_Scaled(TRANS_Old, EMIS_Old, seq, C);

   T = size(seq, 2);
   N = size(TRANS_Old, 2);

   Xi = zeros(N, N, T-1);
         
%L_theta_Scaled = Alpha(:, 1)' * TRANS_Old * diag(EMIS_Old(:, seq(2))) * Beta(:, 2); 
%L_theta_Scaled is the constant 1

for ii = 1 : T-1
    Xi(:, :, ii) = Alpha(:, ii) * Beta(:, ii+1)' .* EMIS_Old(:, seq(ii+1))' .* TRANS_Old;
end

Gamma = squeeze(sum(Xi, 2));

TRANS_New = sum(Xi, 3) ./ sum(Gamma, 2);

Gamma = [Gamma, Alpha(:, T) .* Beta(:, T) / C(T)];

Y = zeros(size(EMIS_Old, 2), T);
for ii = 1 : T
    Y(seq(ii), ii) = 1;
end

 EMIS_New = Gamma * Y' ./ sum(Gamma, 2);
  ini_New = Gamma(:, 1)';
