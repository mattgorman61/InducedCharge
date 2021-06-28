function [x_supp,y_supp,z_supp] = supportfunct(s1pts,s2pts,dir)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[x1,y1,z1] = furthestPtfunct(s1pts,dir);
[x2,y2,z2] = furthestPtfunct(s2pts,-dir);

x_supp = x1-x2;
y_supp = y1-y2;
z_supp = z1-z2;

end

