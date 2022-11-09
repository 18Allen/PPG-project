function [mean_phase, mean_module] = RLS_func(x)
if size(x,2) > size(x,1)
x = x.';
end
d = x(2:end);
p = 8;
lambda = 0.98;
delta = 1;

w = zeros([p,1]);
N = length(x)-1;

% Initialization
w(end) = 1;
P = delta*eye(p);
HFroots = zeros([N-p,1]);
%%
error_list = zeros(size(p:N));
change_list = zeros(size(p:N));
for i =p:N
     e = d(i) - w.'*x(i-p+1:i);
     g = P*x(i-p+1:i)/(lambda+x(i-p+1:i).'*P*x(i-p+1:i));
     P = (eye(p)-g*x(i-p+1:i).')*P/lambda;  
     w = w + e*g;
     
     change_list(i-p+1) = norm(e*g);
     error_list(i-p+1) = abs(e);
     rts = roots([-1;flip(w)]);
     newrts = rts(imag(rts) > 0);
     if i > 2*p
     [val, idx] = min(real(rts(imag(rts) > 0)));
        if isempty(idx)
            HFroots(i) = 0;
        else
            HFroots(i) = newrts(idx);
        end
     end
end
%%
mean_phase = mean(angle(HFroots(p+1:end)));
mean_module = mean(abs(HFroots(p+1:end)));
end