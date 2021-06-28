function [sameDirFlag] = sameDir(a,b)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if(dot(a,b)>0)
    sameDirFlag = 1;
else
    sameDirFlag = 0;
end

end

