function [symbols,max_idx,nmi_list] = Feature2sym(n_label,label,predictors)
max_n = 20;
min_n = 5;
nmi_list = zeros([1,max_n-min_n+1]);
max_nmi = -100;
max_idx = 0;
for i = min_n:max_n
    K_label = kmeans(predictors,i);
    nmi_list(i) = NMI(n_label,label,i,K_label);
    if nmi_list(i) > max_nmi
        max_nmi = nmi_list(i);
        max_idx = i;
        symbols = K_label;
    end
end