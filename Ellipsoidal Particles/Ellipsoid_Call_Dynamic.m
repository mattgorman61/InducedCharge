clc; clear all; close all;

%% Create Multiple Ellipsoids

a = [1 1 1]; b = [2 2 2]; c = [1 1 1];
dxes = [2.5 0 1]; dyes = [0 0 1]; dzes = [0 0 0];
psiVect = [0 0 0]; % Rotation about z-axis
phiVect = [0 0 0]; % Rotation about x-axis
thetaVect = [0 0 0]; % Rotation about y-axis
NpatchesEll = 1000; 
NEll = 2;
Npatches = NEll*NpatchesEll;
EllPatchData = zeros(Npatches,17); % X,Y,Z,dA, nVx,nVy,nVz, a,b,c, x0,y0,z0, psi,phi,theta, ellID of each patch
EulRotMatsData = zeros(4,4,NEll);

frameNum = 0;
maxSteps = 50;
collflag = -1;

while (collflag < 0 && frameNum <= maxSteps)
frameNum = frameNum+1;
for n = 1:NEll
    %[x,y,z,dA,dAmat,nVect,ellID,a1,b1,c1,x0,y0,z0,psi,phi,theta,A_Eul,EllPatchData] = F_createEllipsoid(a(n),b(n),c(n),NpatchesEll,NEll,EllPatchData,dxes(n),dyes(n),dzes(n),psi(n),phi(n),theta(n),n);
                                                                                                        %a,b,c,NpatchesEll,NEll,EllPatchData,dx,dy,dz,psi,phi,theta,ellID_in
    
    [x,y,z,dA,dAmat,nVect,ellID,a1,b1,c1,x0,y0,z0,psi,phi,theta,A_Eul,EllPatchData] = F_createEllipsoid(a(n),b(n),c(n),NpatchesEll,NEll,EllPatchData,dxes(n),dyes(n),dzes(n),psiVect(n),phiVect(n),thetaVect(n),n);
    EulRotMatsData(:,:,n) = A_Eul;
    %{
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),1) = x;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),2) = y;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),3) = z;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),4) = dA;
    
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),5) = nVect(:,1);
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),6) = nVect(:,2);
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),7) = nVect(:,3);
    
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),8) = a1;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),9) = b1;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),10) = c1;
    
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),11) = x0;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),12) = y0;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),13) = z0;    
    
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),14) = ellID;
    %}
end

colors = ['r','b','k','g','m'];
figure(frameNum);
hold on;
for n = 1:NEll
    x_i = EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),1);
    y_i = EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),2); 
    z_i = EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),3);
    scatter3(x_i,y_i,z_i,'filled',colors(n));
    % quiver3(x_i,y_i,z_i,nVect(:,1),nVect(:,2),nVect(:,3),4); % Show Normal Vectors
end

axis equal;
grid on;
%view(0,90);
xlabel('x','fontsize',24,'fontname','Times New Roman','fontangle','Italic'); 
ylabel('y','fontsize',24,'fontname','Times New Roman','fontangle','Italic'); 
zlabel('z','fontsize',24,'fontname','Times New Roman','fontangle','Italic');
set(gca, 'LineWidth', 2.0 );
set(gca, 'fontsize', 24.0 );
set(gca, 'XMinorTick', 'on');
set(gca, 'Ticklength', [0.02;0.01] );
set(gca, 'YMinorTick', 'on');
set(gca, 'Ticklength', [0.02;0.01] );

%for i=1:tsteps
[collisionList,cPoint_1,cPoint_2] = F_CollideCheckCall(EllPatchData, NEll, NpatchesEll,EulRotMatsData)

if(length(collisionList)>0)
    
    scatter3(cPoint_1(1),cPoint_1(2),cPoint_1(3), 500, 'filled', 'g');
    scatter3(cPoint_2(1),cPoint_2(2),cPoint_2(3), 500, 'filled', 'm');
    
    %annotation('textbox', [0.5, 0.5, 0.25, 0.25], 'string', 'COLLISION DETECTED','FitBoxToText','on','Color','r');
    txtbox = annotation('textbox',[0.6 0.5 0 0],'string',{'COLLISION','DETECTED'},'FitBoxToText','on','Color','r','FontSize',18,'fontangle','italic','fontname','Times New Roman');
    txtbox.BackgroundColor = 'w';
    txtbox.EdgeColor = 'k';
    collflag = 1;
    
end
%end

EllMat.id = EllPatchData(:,17);
EllMat.x = EllPatchData(:,1);
EllMat.y = EllPatchData(:,2);
EllMat.z = EllPatchData(:,3);
EllMat.dA = EllPatchData(:,4);
EllMat.nvX = EllPatchData(:,5);
EllMat.nvY = EllPatchData(:,6);
EllMat.nvZ = EllPatchData(:,7);
EllMat.La = EllPatchData(:,8);
EllMat.Lb = EllPatchData(:,9);
EllMat.Lc = EllPatchData(:,10);
EllMat.x0 = EllPatchData(:,11);
EllMat.y0 = EllPatchData(:,12);
EllMat.z0 = EllPatchData(:,13);
EllMat.psi = EllPatchData(:,14);
EllMat.phi = EllPatchData(:,15);
EllMat.theta = EllPatchData(:,16);
EllMat.rotMat = EulRotMatsData(:,:,:);

 phiVect(1) = phiVect(1) + pi/2/10; % Incrementally rotate ellipse 1, 90 degrees over 20 frames
%dxes(2) = dxes(2) + 2.5/10;

if(collflag<0)
    pause;
end
end
