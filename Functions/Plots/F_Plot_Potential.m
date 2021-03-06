function [finished] = F_Plot_Potential(R,x,y,z,dA,nVect,x_pcs,y_pcs,z_pcs,pcharge,lpCharge,lEField,sigma_b,sigma_f,k_air,k_obj,epsilon_0)
% DISPLAYS FORCE VECTORS ON EACH PATCH
%   

% Add Functions Folder to the path
currDir = pwd;
%fprintf('%s',currDir);
idcs   = strfind(currDir,'\');
newdir = currDir(1:idcs(end)-1);
addpath (newdir);


sigma = sigma_f + sigma_b;
[phi,phi_0,phi_norm,theta] = F_getPotentials_Matrix(R,x,y,z,dA,nVect,x_pcs,y_pcs,z_pcs,pcharge,sigma,k_air,k_obj,epsilon_0);

% Display Potential Values
scatter3(x_pcs,y_pcs,z_pcs,12,'filled','k');
hold on;
axis equal;
scatter3(x,y,z,12,phi_norm,'filled');
%xlim([-2*R, 2*R]); ylim([-2*R, 2*R]); zlim([-2*R, 2*R]);


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

%%{
% Purple/Gray
for i = 1:numLevels
    cmap_cust(i,1) = 0.9 + 0.1*i/numLevels; 
    cmap_cust(i,2) = 0.9 - 0.9*i/numLevels;
    cmap_cust(i,3) = 0.9 + 0.1*i/numLevels;
end
%}


colormap(cmap_cust);

%}

cbar = colorbar;
cbarTitle = title(cbar, '\it{\psi}');
cbarTitle.FontName = 'Times New Roman';
cbarTitle.FontSize = 16;

% Below is for Force calculations....
%{
tit2 = title('Force Vectors'); tit2.FontSize = 12; 
tit2.FontName = 'Times New Roman';
view(35,20);

quiver3(x,y,z,Fx/F0,Fy/F0,Fz/F0,3,'b');
q = quiver3(0,0,0,Fnet(1)/F0,Fnet(2)/F0,Fnet(3)/F0,3,'k');
q.MaxHeadSize = 2;
axis equal; 
legend('Point Charge','Surface Bound Charges','Patch Forces','Net Force','Location','southeast')
%}

finished = true;

end

