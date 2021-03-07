function [xc,yc,zc,NormVec,DeltaArea,NN] = F_GeometrySphere(Radius,N)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
    a = (4*pi*Radius^2)/N; % Average area per patch
    d = sqrt(a); % Length of patch
    M_Phi = round(pi*Radius/d);
    d_Phi = pi*Radius/M_Phi;
    d_Theta = a/d_Phi;
    
    % The center position of each patch
    xc_temp = linspace(0,0,N)';
    yc_temp = linspace(0,0,N)';
    zc_temp = linspace(0,0,N)';
    % Normal unit vector
    NormVec_x_temp = linspace(0,0,N)';
    NormVec_y_temp = linspace(0,0,N)';
    NormVec_z_temp = linspace(0,0,N)';
    
    count = 0;
    for i = 0:(M_Phi-1)
    
    Phi = pi*(i+0.5)/M_Phi;
    M_Theta = round(2*pi*Radius*sin(Phi)/d_Theta);
    
    % For debug
%     fprintf('==============================\n');
%     fprintf('Phi = %f, z = %f\n',Phi,Radius*cos(Phi));
    
        for j = 0:(M_Theta-1)
            count = count + 1;

            Theta = 2*pi*j/M_Theta;
            xc_temp(count)=Radius*sin(Phi)*cos(Theta);
            yc_temp(count)=Radius*sin(Phi)*sin(Theta);
            zc_temp(count)=Radius*cos(Phi);
            
            % For debug
%             fprintf('Theta = %f, x = %f,y = %f,z = %f\n',Theta,...
%                 Radius*sin(Phi)*cos(Theta),Radius*sin(Phi)*sin(Theta),...
%                 Radius*cos(Phi));
            
            % Previous definition of unit vector seems to be problematic
            %NormVec_x_temp(count) = sin(Phi)*cos(Theta);
            %NormVec_z_temp(count) = sin(Phi)*sin(Theta);
            %NormVec_y_temp(count) = cos(Phi);
            
%             NormVec_x_temp(count) = xc_temp(count)/Radius;
%             NormVec_z_temp(count) = yc_temp(count)/Radius;
%             NormVec_y_temp(count) = zc_temp(count)/Radius;            
        end % for j = 0:(M_Theta-1)
    
    end % for i = 0:(M_Phi-1)
    
    % Approximated patch area
    NN = count; % Number of generated patches
    
    xc = xc_temp(1:NN);
    yc = yc_temp(1:NN);
    zc = zc_temp(1:NN);

    NormVec(:,1) = xc_temp(1:NN)/Radius;
    NormVec(:,2) = yc_temp(1:NN)/Radius;
    NormVec(:,3) = zc_temp(1:NN)/Radius;

%     NormVec(:,1) = NormVec_x_temp(1:NN);
%     NormVec(:,2) = NormVec_y_temp(1:NN);
%     NormVec(:,3) = NormVec_z_temp(1:NN);
    
    DeltaArea = (4*pi*Radius^2)/NN * linspace(1,1,NN)'; % Average area (should be improved)
end

