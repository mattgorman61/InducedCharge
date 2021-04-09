function [collisionList,ContactPoint_1,ContactPoint_2] = F_CollideCheckCall(EllPatchData, NEll, NpatchesEll,EulRotMatsDataIn)
% CHECKS WHETHER ELLIPSES HAVE COLLIDED
%   Alfano and Greer 2003 Journ. of Guidance, Control, and Dynamics

%a = 1; b = 2; c = 1;
%x0 = 1; y0 = 1; z0 = 1;

ContactPoint_1 = 0;

list = [];
for i = 1:NEll
    for j = i:NEll
        if(j~=i)
            a0 = EllPatchData(i*NpatchesEll,8);
            b0 = EllPatchData(i*NpatchesEll,9);
            c0 = EllPatchData(i*NpatchesEll,10);
            x0 = EllPatchData(i*NpatchesEll,11);
            y0 = EllPatchData(i*NpatchesEll,12);
            z0 = EllPatchData(i*NpatchesEll,13);
            AEul0 = EulRotMatsDataIn(:,:,i);
            ellID0 = EllPatchData(i*NpatchesEll,17);
            
            a1 = EllPatchData(j*NpatchesEll,8);
            b1 = EllPatchData(j*NpatchesEll,9);
            c1 = EllPatchData(j*NpatchesEll,10);
            x1 = EllPatchData(j*NpatchesEll,11);
            y1 = EllPatchData(j*NpatchesEll,12);
            z1 = EllPatchData(j*NpatchesEll,13);
            AEul1 = EulRotMatsDataIn(:,:,j);
            ellID1 = EllPatchData(j*NpatchesEll,17);
            
            [lColl_1,ContactPoint_1] = F_CollideCheck(a0,b0,c0,x0,y0,z0,AEul0, a1,b1,c1,x1,y1,z1,AEul1);
            [lColl_2,ContactPoint_2] = F_CollideCheck(a1,b1,c1,x1,y1,z1,AEul1, a0,b0,c0,x0,y0,z0,AEul0);
            
            if(lColl_1 == true)
                fprintf('\nCOLLISION DETECTED. ELLIDs: %g, %g\n\n',ellID0,ellID1);
                list = [list, ellID0];
                list = [list, ellID1];
            end
        end
    end
end


collisionList = list;



end

