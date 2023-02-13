function col(V, GT)

my_colormap=[1, 0, 0;
             1, 1, 0;
             0, 1, 0;
             0, 0, 1;
             0, 0, 0];
        
labels = {'Awake', 'REM', 'N1', 'N2', 'N3'};
figure;
scatter3(V(:, 1), V(:, 2), V(:, 3), 3, GT);
colormap(my_colormap);
lcolorbar(labels, 'fontweight', 'bold');
axis tight;
title(name_title);
        
        