function [xc,yc,zc,DeltaArea,N,NormVec,a,b,c] = F_GeometryEllip(filename)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
    
    load(filename);
    
    % Approximated area of the ellipsoid
    S = 2*pi*sqrt(a^2*b^2 + b^2*c^2 + c^2*a^2 + (a^2*b*c + a*b^2*c + a*b*c^2)/3);
    
    % The center position of each patch
    xc = p(:,1);
    yc = p(:,2);
    zc = p(:,3);
    
    N = length(xc);
    NormVec = zeros(N,3);
    
    % Patch location
    
    for i = 1:N
        % Normal vector
        Denom = sqrt(xc(i)^2/a^4+yc(i)^2/b^4+zc(i)^2/c^4);
        NormVec(i,1) = xc(i)/a^2/Denom;
        NormVec(i,2) = yc(i)/b^2/Denom;
        NormVec(i,3) = zc(i)/c^2/Denom;
    end
    
    % Approximated patch area
    DeltaArea = S/N * linspace(1,1,N)'; % Average area (should be improved)
end

