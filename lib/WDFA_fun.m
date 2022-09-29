function [W] = WDFA_fun(data,pts, sigma,order)
% -----------------------------------------------------
% DESCRIPTION:
% Function for the WDFA analysis.
%
% INPUTS: 
% data: a one-dimensional data vector.
% pts: sizes of the window at which to evaluate the fluctuation
% sigma:  the size for energy calcuation
%
% OUTPUTS: 
% W: A vector of size Nx1 containing the fluctuations corresponding to the
% windows specified in entries in pts. The beginning and the ending sigma
% element is not defined and thus set to 0.
% -----------------------------------------------------ß

%%
% in column
sz = size(data);
if sz(1)< sz(2)
    data = data';
end


% WDFA
N = length(data);
W = zeros([N,1]);
    
w = pts;

n = floor(N/w);
Nfloor = n*pts;
D = data(1:Nfloor);

y = cumsum(D-mean(D));
%%
bin = 0:w:(Nfloor-1);
vec = 1:w;

coeff = arrayfun(@(j) polyfit(vec',y(bin(j) + vec),order),1:n,'uni',0);
y_hat = cell2mat(cellfun(@(y) polyval(y,vec),coeff,'uni',0));
zn = y - y_hat';
mu1 = cumsum(zn);
mu2 = cumsum(zn.^2);

idx = sigma+1: length(D)-sigma;
W(idx) = (1/(2*sigma))*(mu2(idx+sigma) - mu2(idx -sigma)) ...
    - (1/(2*sigma))^2*(mu1(idx+sigma) - mu1(idx -sigma)).^2;

end