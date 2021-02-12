function [finished] = F_Plot_NormVectors(R,x,y,z,nVect)
% DISPLAYS NORMAL VECTORS EACH PATCH
%   

% Add Functions Folder to the path
currDir = pwd;
%fprintf('%s',currDir);
idcs   = strfind(currDir,'\');
newdir = currDir(1:idcs(end)-1);
addpath (newdir);

% Display the Normal Vectors 

scatter3(x,y,z,'filled','k');
hold on;


quiver3(x,y,z,nVect(:,1),nVect(:,2),nVect(:,3),3,'r');

xlim([-2*R, 2*R]); ylim([-2*R, 2*R]); zlim([-2*R, 2*R]);
axis square; 
tit1 = title('Normal Vectors'); tit1.FontName = 'Times New Roman';
tit1.FontSize = 12;
%set(gca,'LineWidth',1.5); set(gcf,'Position',[100,400,500,400]);

    
    
finished = true;



end

