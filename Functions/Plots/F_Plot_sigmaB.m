function [finished] = F_Plot_sigmaB(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma_b,k_air,k_obj,epsilon_0)
% DISPLAYS SIGMA_B, BOUND SURFACE CHARGE DENSITY, FOR EACH PATCH
%   

% Add Functions Folder to the path
currDir = pwd;
%fprintf('%s',currDir);
idcs   = strfind(currDir,'\');
newdir = currDir(1:idcs(end)-1);
addpath (newdir);

% Plot sigma_b

scatter3(x,y,z,12,sigma_b,'filled');
hold on; 
scatter3(x_pc,y_pc,z_pc,12,'filled','k');
xlim([-2*R, 2*R]); ylim([-2*R, 2*R]); zlim([-2*R, 2*R]);
colorbar;
axis square;    
tit2 = title('Surface Bound Charge'); tit2.FontSize = 12; 
tit2.FontName = 'Times New Roman';
%set(gca,'LineWidth',1.5); set(gcf,'Position',[100,100,500,400]);
view(35,20);

    
    
finished = true;



end

