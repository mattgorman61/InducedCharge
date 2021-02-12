function [finished] = F_Plot_ForceVectors(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma_b,sigma_f,k_air,k_obj,epsilon_0)
% DISPLAYS FORCE VECTORS ON EACH PATCH
%   

% Add Functions Folder to the path
currDir = pwd;
%fprintf('%s',currDir);
idcs   = strfind(currDir,'\');
newdir = currDir(1:idcs(end)-1);
addpath (newdir);


sigma = sigma_f + sigma_b;
[Fnet,Fx,Fy,Fz,F0] = F_getForces_Matrix(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma,k_air,k_obj,epsilon_0);

% Display Force Vectors
scatter3(x_pc,y_pc,z_pc,12,'filled','k');
hold on;
scatter3(x,y,z,12,sigma_b,'filled');
xlim([-2*R, 2*R]); ylim([-2*R, 2*R]); zlim([-2*R, 2*R]);
colorbar;
axis square;    
tit2 = title('Force Vectors'); tit2.FontSize = 12; 
tit2.FontName = 'Times New Roman';
view(35,20);

quiver3(x,y,z,Fx/F0,Fy/F0,Fz/F0,3,'r');
quiver3(0,0,0,Fnet(1)/F0,Fnet(2)/F0,Fnet(3)/F0,0.5,'k');
axis equal; 
legend('Point Charge','Surface Bound Charges','Patch Forces','Net Force','Location','southeast')


finished = true;

end

