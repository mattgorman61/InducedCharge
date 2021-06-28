function [simpFlag, S, dir] = nextsimplexfunct(S,dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Sdim = length(S(:,1) );

if(Sdim==2)
    fprintf('SIMPLEX IS CURRENTLY A LINE\n\n\n');
    [simpFlag, S, dir] = simpLinefunct(S,dir);
else
    if(Sdim==3)
        fprintf('SIMPLEX IS CURRENTLY A TRIANGLE\n\n\n');
        [simpFlag, S, dir] = simpTrifunct(S,dir);
    else
        if(Sdim==4)
            fprintf('SIMPLEX IS CURRENTLY A TETRAHEDRON\n\n\n');
            [simpFlag, S, dir] = simpTetrafunct(S,dir);
        end
    end
end

% simpFlag = -1000;




end

