function [CollideList,ContactPoints,vpar,rot_pf,F_par,M_par] = ...
    F_Collision(npars,Elastic,Poisson,DetectList,a,b,c,...
    rpar,RMatrx,thickness,vpar,rot_pf,frict,F_par,M_par,time,fid,erest,mass)
%
%
    LargeOverlap = 0;
    % Contact point between particle pairs
    ContactPoints = zeros(npars,npars,3); % id of targer, id of another, x/y/z

    % Check each possible particle pair
    CollideList = zeros(npars);

    for i = 1: (npars-1)
        for j = (i+1):npars
            if (DetectList(i,j)==1)
                % Contact point identification
                [CollideList(i,j),Contact1,Contact2,Ematrx_inc]...
                    = F_CollideCheck(i,j,rpar,RMatrx,a,b,c,thickness);
                CollideList(j,i) = CollideList(i,j);
                
                % Use the middle point of two contact points as the "real"
                % contact point
                if (CollideList(i,j))
                ContactPoints(i,j,:)=Contact1;
                ContactPoints(j,i,:)=Contact2;
                                
                % Linear/Angular momentum transfer in collisions
                [Fc1,Fc2,Mc1,Mc2]=F_CollMomTrans(i,j,a,b,c,Elastic,Poisson,...
                    RMatrx,Ematrx_inc,rpar,vpar,rot_pf,frict,Contact1',...
                    Contact2',time,fid,erest,mass,thickness);
                F_par(i,:) = F_par(i,:) + Fc1';
                F_par(j,:) = F_par(j,:) + Fc2';
                M_par(i,:) = M_par(i,:) + Mc1';
                M_par(j,:) = M_par(j,:) + Mc2';

            end % if (CollideList(i,j)==1)
            
        end % for j = (i+1):npars
    end % for i = 1: (npars-1)
end

