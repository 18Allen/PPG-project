function[F] = DFA_fun(data,pts,order)

% -----------------------------------------------------
% DESCRIPTION:
% Function for the DFA analysis.

% INPUTS: 
% data: a one-dimensional data vector.
% pts: sizes of the windows/bins at which to evaluate the fluctuation

% OUTPUTS: 
% F: A vector of size Nx1 containing the fluctuations corresponding to the
% windows specified in entries in pts.
% -----------------------------------------------------ß


% in column
sz = size(data);
if sz(1)< sz(2)
    data = data';
end


% DFA
npts = numel(pts);

F = zeros(npts,1);
N = length(data);


for h = 1:npts
    
    w = pts(h);
    
    n = floor(N/w);
    Nfloor = n*pts(h);
    D = data(1:Nfloor);
    
    y = cumsum(D-mean(D));
    
    bin = 0:w:(Nfloor-1);
    vec = 1:w;
    
    coeff = arrayfun(@(j) polyfit(vec',y(bin(j) + vec),order),1:n,'uni',0);
    y_hat = cell2mat(cellfun(@(y) polyval(y,vec),coeff,'uni',0));
    F(h)  = mean((y - y_hat').^2)^0.5;
end



end
