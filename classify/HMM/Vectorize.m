function Vec = Vectorize(input)

%Transform a 3-dimensional matrix into 2-dimensional matrix by vectorize
%The dimension of parameter and output are m-n-s and s-mn, respectively.

Vec_tem = [];

for ii = 1 : size(input, 1)
Vec_tem = [Vec_tem input(ii, :, :)];
end

Vec = squeeze(Vec_tem)';

