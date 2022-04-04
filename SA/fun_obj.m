%Objective Function
function [misfit] = fun_obj(Us_obs, Us_cal)
    ls = length(Us_obs);
    for j = 1 : ls
        m(j) = ((Us_obs(j) - Us_cal(j)))^2;
    end
    misfit = sqrt((1/ls)*sum(m));
end