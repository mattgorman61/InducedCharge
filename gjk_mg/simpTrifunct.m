function [triBool,S,dir] = simpTrifunct(S_in,dir_in)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

S = S_in;
dir = dir_in;

a = S_in(1,:);
b = S_in(2,:);
c = S_in(3,:);

ab = b - a;
ac = c - a;
ao = -a;

abc = cross(ab,ac);
% triBool = -10001;

if(sameDir(cross(abc,ac),ao))
    if(sameDir(ac,ao))
        S = [a;c];
        dir = cross(cross(ac,ao),ac);
        triBool = false;
    else
        triBool = simpLinefunct([a;b],dir);
    end
else
    if(sameDir(cross(ab,abc),ao))
        triBool = simpLinefunct([a;b],dir);
    else
        if(sameDir(abc,ao))
            dir = abc;
            triBool = false;
        else
            S = [a;c;b];
            dir = -abc;
            triBool = false;
        end
    end
end



