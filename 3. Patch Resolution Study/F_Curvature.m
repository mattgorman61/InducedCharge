function [Curv1,Curv2] = F_Curvature(Contact1_pf,Contact2_pf,a,b,c)
% Calculate the curvature at contact point, which is used for normal 
% contact force calculation
    x1 = Contact1_pf(1); y1 = Contact1_pf(2); z1 = Contact1_pf(3);
    x2 = Contact2_pf(1); y2 = Contact2_pf(2); z2 = Contact2_pf(3);
    
    h1 = (x1^2/a^4 + y1^2/b^4 + z1^2/c^4)^(-1/2);
    h2 = (x2^2/a^4 + y2^2/b^4 + z2^2/c^4)^(-1/2);
    
    Curv1 = 0.5*(h1^3)*((x1^2/a^2 + y1^2/b^2)/a^2/b^2 + ...
        (y1^2/b^2 + z1^2/c^2)/b^2/c^2 + (z1^2/c^2 + x1^2/a^2)/c^2/a^2);
    Curv2 = 0.5*(h2^3)*((x2^2/a^2 + y2^2/b^2)/a^2/b^2 + ...
        (y2^2/b^2 + z2^2/c^2)/b^2/c^2 + (z2^2/c^2 + x2^2/a^2)/c^2/a^2);
end

