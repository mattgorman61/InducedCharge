function [tetraBool,S,dir] = simpTetrafunct(S_in,dir_in)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

S = S_in;
dir = dir_in;

a = S_in(1,:);
b = S_in(2,:);
c = S_in(3,:);
d = S_in(4,:);

ab = b - a;
ac = c - a;
ad = d - a;
ao = -a;

abc = cross(ab,ac);
acd = cross(ac,ad);
adb = cross(ad,ab);

% tetraBool = -10001;

if(sameDir(abc,ao))
    %S = [a;b;c;d]
    tetraBool = simpTrifunct([a;b;c],dir);
else
    if(sameDir(acd,ao))
        tetraBool = simpTrifunct([a;c;d],dir);
    else
        if(sameDir(adb,ao))
            tetraBool = simpTrifunct([a;d;b],dir);
        else
            tetraBool = true;
        end
        
    end
end

    
    
end







