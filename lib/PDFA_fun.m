function [P] = PDFA_fun(data,pts,order)
% in column
sz = size(data);
if sz(1)< sz(2)
    data = data';
end

% PDFA
P = zeros(size(data));
N = length(data);


w = pts;
n = floor(N/w);
Nfloor = n*pts;
D = data(1:Nfloor);

y = cumsum(D-mean(D));

bin = 0:w:(Nfloor-1);
vec = 1:w;

coeff = arrayfun(@(j) polyfit(vec',y(bin(j) + vec),order),1:n,'uni',0);
y_hat = cell2mat(cellfun(@(y) polyval(y,vec),coeff,'uni',0));
P = (y - y_hat').^2;
P = cumsum(P).^0.5;

end