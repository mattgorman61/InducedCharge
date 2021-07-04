function [v] = F_createVector_bias(bias,numPoints,L)

N = numPoints-1;
dx = L/N;
llll = L/N * 2 / (bias+1);
dllll = linspace(((bias-1)/2*llll),0,N/2);

v = zeros(N+1,1);

v(1) = 0;
v(N+1) = L;
% v(1) = (bias+1)*llll - dllll(1);
% v(N+1) = 1/2*(bias+1)*llll + dllll(1);

for i = 2:N/2+1
    dl = dllll(i-1);
    v(i) = v(i-1) + (1/2*(bias+1)*llll - dl);
    v(N+2-i) = v(N+3-i) - (1/2*(bias+1)*llll + dl);
end


end

