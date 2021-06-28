function [DetectList] = F_ContactThreatList_CONSTRUCTION(npars,rpar,a,b,c)
%F_COLLIDELIST
%{
    Returns a logical matrix of particle pairs that are about to collide.

%}
    Distance = pdist(rpar); % Interparticle distance
    DistMatrx = squareform(Distance); % Distance matrix
    
    AxisLong = max(a,max(b,c)); % Use longest axis as criterion
    DetectList = (DistMatrx <= (3*AxisLong)); % Label close pair
    
    DetectList = DetectList - eye(npars); % Remove particle itself
end
