function [finished] = F_Plot_sigmaB(x,y,z,x_pcs,y_pcs,z_pcs,sigma_b,lpCharge,lEField,Ext_EField_x,Ext_EField_y,Ext_EField_z)
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

if(lpCharge)
    scatter3(x_pcs,y_pcs,z_pcs,12,'filled','k');
end
%xlim([-2*R, 2*R]); ylim([-2*R, 2*R]); zlim([-2*R, 2*R]);

%%{
% Custom ColorMaps:
numLevels = 100;
cmap_cust = zeros(numLevels,3);
%%{
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


colormap(cmap_cust);

%}

colorbar;
axis equal;    
tit2 = title('Surface Bound Charge'); tit2.FontSize = 12; 
tit2.FontName = 'Times New Roman';
%set(gca,'LineWidth',1.5); set(gcf,'Position',[100,100,500,400]);
view(35,20);
xlabel('x','FontName','Times New Roman');
ylabel('y','FontName','Times New Roman');
zlabel('z','FontName','Times New Roman');

if(lEField)
    numDiv = 6;
    [xx,yy,zz] = meshgrid(linspace(min(x),max(x),numDiv),linspace(min(y),max(y),numDiv),linspace(min(z),max(z),numDiv));
    hold on;
    quiver3(xx,yy,zz, ... 
       Ext_EField_x*ones(numDiv,numDiv,numDiv),Ext_EField_y*ones(numDiv,numDiv,numDiv),Ext_EField_z*ones(numDiv,numDiv,numDiv), ...
       0.5,'g');
    %quiver3(0,0,0,1,1,1,0.25,'k','LineWidth',2);
    axis equal;
end

if(lpCharge)
    if(lEField)
        legend('Surface Charge Density','Point Charge','External E-Field','Location','south');
    else
        legend('Surface Charge Density','Point Charge','Location','south');
    end
else
    if(lEField)
        legend('Surface Charge Density','External E-Field','Location','south');
    else
        legend('Surface Charge Density','Location','south');
    end
end


    
finished = true;



end

