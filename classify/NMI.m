function NMI_val = NMI(n_W,lamb_W, n_K, lamb_K)
%%% This program normalized mutual information between the labels 
%%% (class) of each feature vector and and the corresponding closest
%%% clusters
%%% n_W : #labels (given)
%%% lamb_W : label of each datapoint. The range of value is 1,2,...,n_W
%%% n_K : #clusters (setted)
%%% lamb_K : index of cluster calculated by k-means. The range of value is
%%% 1,2,...,n_k

n = size(lamb_W,1);
HX = 0;
for i=1:n_W
    n_i = sum(lamb_W == i);
%     if n_i == 0
%         error(strcat('No element for label '),num2str(i) ,'!!');
%     end
    tmp  =- n_i*log(n_i/n);
    if isnan(tmp)
        tmp = 0;
    end
    HX = HX +tmp;
end


HY = 0;
for i=1:n_K
    n_i = sum(lamb_K == i);
%     if n_i == 0
%         error(strcat('No element for cluster '),num2str(i) ,'!!');
%     end
    tmp = - n_i*log(n_i/n);
    if isnan(tmp)
        tmp = 0;
    end
    HY = HY + tmp;
end

IXY = 0;
for i=1:n_W
    for j=1:n_K
    n_i = sum(lamb_W == i);
    n_j = sum(lamb_K == i);
    n_ij = sum(lamb_K == i & lamb_W == i);
%     if n_i == 0
%         error(strcat('No element for intersection ('),num2str(i) ,', ',num2str(j),')!!');
%     end
    
    tmp =  n_ij*log((n*n_ij) / (n_i*n_j));
    if isnan(tmp)
        tmp = 0;
    end
    IXY = IXY + tmp;
    end
end

NMI_val = IXY/sqrt(HX*HY);
