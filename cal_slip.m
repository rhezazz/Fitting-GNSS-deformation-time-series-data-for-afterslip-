%Calculated afterslip  (exponent + Logarithm)
function [Us] = cal_slip(A,B,C,D,V0,Ta,Tb,Tc,t)
    for i = 1 : length(t)
        Us(i) = A.*log(1 + (t(i)./Ta)) + B - C.*exp(-t(i)./Tb) - D.*exp(-t(i)/Tc) + V0;
    end
end