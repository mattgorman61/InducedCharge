function [xmax,ymax,zmax,maxdotprod] = furthestPtfunct(s1pts,dir)
%GJK Sub Algorithm: Support Function
%   Given points of shape1, finds farthest point in direction dir
%   
% 
% Given:
% dir ................ [x,y,z]
% s1pts .............. [x1,y1,z1] (x1,y1,z1 are column vectors, coordinates of the shape vertices)
% 
% Returns;
% xmax .............. x location of farthest point on shape1 in direction dir
% ymax .............. y location of farthest point on shape1 in direction dir
% zmax .............. z location of farthest point on shape1 in direction dir
% 

maxdotprod = -10^20;

for i =1:length(s1pts)
%     dotprod = s1pts(i,1)*dir(1) + s1pts(i,2)*dir(2) + s1pts(i,3)*dir(3);
    dotprod = dot(s1pts(i,:),dir');
    if (dotprod>maxdotprod)
        maxdotprod = dotprod;
        xmax = s1pts(i,1);
        ymax = s1pts(i,2);
        zmax = s1pts(i,3);
    end
end


end

