function [te]=teager_energy_func(RRI_res)
 te = RRI_res(2:end-1).^2 - RRI_res(1:end-2).*RRI_res(3:end);
end