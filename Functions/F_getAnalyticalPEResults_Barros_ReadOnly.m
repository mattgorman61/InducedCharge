function [PE_sphere] = F_getAnalyticalPEResults_Barros_ReadOnly(k_med,pcharge,k_tilda,d,R)
% Returns the PE of a neutral sphere near a point charge q
%   PE is normalized, in units of 
N = 6000;
sum = 0;
epsilon_0 = 8.85*10^(-12);

for n = 0:N
    num = (1 - k_tilda)*n;
    denom = ((1+k_tilda)*n + 1)*((1+ (d))^(2*(n+1)));
    sum = sum + (num/denom);
end

% Not Normalized:
PE_sphere1 = (pcharge^2)/8/pi/epsilon_0/k_med/R*sum;

% Normalized:
PE_sphere = PE_sphere1 * epsilon_0*k_med*R/(pcharge^2);

end

