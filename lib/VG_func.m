function [deg_mean, deg_std, ast, per_low, per_high, cc_mean, cc_std] = VG_func(RRI)

    %% 
    % build connection
    connection_table = zeros([length(RRI),length(RRI)]);
    for i = 1:length(RRI)-1
        max_Rslope = -inf; %min_Lslope = inf;
        for j = i+1:length(RRI)
            new_Rslope = (RRI(j) - RRI(i))/sum(RRI(i+1:j));
            if new_Rslope > max_Rslope
                max_Rslope = new_Rslope;
                connection_table(i,j) = 1;
                connection_table(j,i) = 1;
            else
               
            end
             % Left slope will be covered by previous ones
        end
    end
    
    % Calculate degree
    deg = sum(connection_table,2);
    deg_mean = mean(deg);
    deg_std = std(deg);
    %N = N+histcounts(sum(connection_table,2),1:60)/(length(RRI)*length(idx));
    
    % assortativity coefficient
    sum_ab = 0;
    sum_a2_plus_b2 = 0;
    sum_a_plus_b = 0;
    n_edges = 0;
    for i = 1:length(RRI)-1
        idx = find(connection_table(i,i+1:end))+i;
        % deg will always greater than 1
        if ~isempty(idx) 
        n_edges = n_edges + length(idx);
        sum_ab  = sum_ab + sum(deg(i)*deg(idx));
        sum_a2_plus_b2 = sum_a2_plus_b2 + sum(deg(i)^2 + deg(idx).^2);
        sum_a_plus_b = sum_a_plus_b + sum(deg(i) + deg(idx));
        
        end
    end
    M = sum(deg)/2;
    ast = (M*sum_ab - (sum_a_plus_b/2)^2)/ (M*sum_a2_plus_b2/2  - (sum_a_plus_b/2)^2);
    % percentage of low/high
    per_low = sum(deg <= 3)/length(RRI);
    per_high = sum(deg >= 11)/length(RRI);
    %%
    % cluster coefficient
    cc = zeros(size(RRI(2:end-1)));
    for i = 2:length(RRI)-1 % sometimes the first/last will have only 1 neighbor
        idx = find(connection_table(i,:));
        ki = length(idx);
        cc(i) = (sum(connection_table(idx,idx),'all'))/(ki*(ki-1));
    end
    cc_mean = mean(cc);
    cc_std = std(cc);

end