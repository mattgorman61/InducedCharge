function [CollideList,ContactPoints] = F_Collision(npars,DetectList,a,b,c,...
    rpar,RMatrx,thickness)
%
%   
    % Contact point between particle pairs
    ContactPoints = zeros(npars,npars,3); % id of targer, id of another, x/y/z

    % Check each possible particle pair
    CollideList = zeros(npars);

    for i = 1: (npars-1)
        for j = (i+1):npars
            if (DetectList(i,j)==1)
                [CollideList(i,j),Contact1,Contact2]...
                    = F_CollideCheck(i,j,rpar,RMatrx,a,b,c,thickness);
                CollideList(j,i) = CollideList(i,j);
                ContactPoints(i,j,:)=Contact1;
                ContactPoints(j,i,:)=Contact2;
            end % if (CollideList(i,j)==1)
        end % for j = (i+1):npars
    end % for i = 1: (npars-1)
end

