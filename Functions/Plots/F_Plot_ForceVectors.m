function [finished] = F_Plot_ForceVectors(numSpheres,NpatchesSph,R,x,y,z,dA,nVect,x_pcs,y_pcs,z_pcs,dxs,dys,dzs,pcharge,lpCharge,lEField,sigma_b,sigma_f,k_air,k_obj,epsilon_0)
% DISPLAYS FORCE VECTORS ON EACH PATCH
%   

% Add Functions Folder to the path
currDir = pwd;
%fprintf('%s',currDir);
idcs   = strfind(currDir,'\');
newdir = currDir(1:idcs(end)-1);
addpath (newdir);


sigma = sigma_f + sigma_b;
[Fnet,Fx,Fy,Fz,F0] = F_getForces_Mult_Matrix(numSpheres,NpatchesSph,R,x,y,z,dA,nVect,x_pcs,y_pcs,z_pcs,pcharge,sigma,k_air,k_obj,epsilon_0);

% Display Force Vectors
if(lpCharge)
    scatter3(x_pcs,y_pcs,z_pcs,12,'filled','k');
end
hold on; axis equal;
scatter3(x,y,z,12,sigma_b,'filled');
% xlim([-2*R, 2*R]); ylim([-2*R, 2*R]); zlim([-2*R, 2*R]);

%%{
% Custom ColorMaps:
numLevels = 100;
cmap_cust = zeros(numLevels,3);
%{
% Red/Gray
for i = 1:numLevels
    cmap_cust(i,1) = 0.9 + 0.1*i/numLevels; 
    cmap_cust(i,2) = 0.9 - 0.9*i/numLevels;
    cmap_cust(i,3) = 0.9 - 0.9*i/numLevels;
end
%}

%{
% Purple/Gray
for i = 1:numLevels
    cmap_cust(i,1) = 0.9 + 0.1*i/numLevels; 
    cmap_cust(i,2) = 0.9 - 0.9*i/numLevels;
    cmap_cust(i,3) = 0.9 + 0.1*i/numLevels;
end
%}


% colormap(cmap_cust);

%}

colorbar;
quiver3(x,y,z,Fx,Fy,Fz,5,'r');
hold on;
q1 = quiver3(dxs(1),dys(1),dzs(1),Fnet(1,1)/F0,Fnet(1,2)/F0,Fnet(1,3)/F0,80,'k');
q1.MaxHeadSize = 2;
hold on;
if(numSpheres>1)
    q2 = quiver3(dxs(2),dys(2),dzs(2),Fnet(2,1)/F0,Fnet(2,2)/F0,Fnet(2,3)/F0,80,'k');
    q2.MaxHeadSize = 2;
end
grid on; box on;

for i = 1:numSpheres
    fprintf('\nForces on Particle %d:\n', i);
    fprintf('\tFx: %.4f \n',Fnet(i,1));
    fprintf('\tFy: %.4f \n',Fnet(i,2));
    fprintf('\tFz: %.4f \n',Fnet(i,3));
    fprintf('\n\n');
end 


if(lpCharge)
    legend('Point Charge','Surface Bound Charges','Patch Forces','Net Force','Location','southeast');
else
    legend('Surface Bound Charges','Patch Forces','Net Force','Location','southeast')
end

hold on;
tit2 = title('Force Vectors'); tit2.FontSize = 12; 
tit2.FontName = 'Times New Roman';
view(35,20);

finished = true;

end

