function [lineBool,S,dir] = simpLinefunct(S_in,dir_in)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

S = S_in;

a = S_in(1,:);
b = S_in(2,:);

ab = b - a;
ao = -a;

lineBool = 0;

if(sameDir(ab,ao))
    dir = cross(cross(ab,ao),ab);
else
    S = a;
    dir = ao;
end

end

