function [F1,F2,M1,M2] = F_CollMomTrans(mm,nn,a,b,c,Elastic,Poisson,RMatrx,...
    Ematrx_inc,rpar,vpar,rot_pf,frict,Contact1,Contact2,time,fid,erest,mass)
%F_COLLMOMTRANS: calculate the linear/angular momentum transfer in each
%collision. Soft sphere model is used here.
    
    % Normal overlap
    delta_n = pdist([Contact1,Contact2]');
    
    % Reduced elastic modulus
    Elastic0 = Elastic/2/(1-Poisson^2);

    % Characteristic matrix of ellipsoid
    Q1 = reshape(Ematrx_inc(1,:,:),[4,4]);
    Q2 = reshape(Ematrx_inc(2,:,:),[4,4]);
    
    % Use quadratic form to find the normal vector
    A1 = Q1(1:3,1:3);
    b1 = 2*Q1(1:3,4);
    
    A2 = Q2(1:3,1:3);
    b2 = 2*Q2(1:3,4);
    
    %(1) Find the normal vector:
    % n = grad(x^T*A_1*x+b1^T*x+c1)/|grad(x^T*A_1*x+b1^T*x+c1)|
    %   = 2*A1*x + b1/|2*A1*x + b1|
    nVec1 = 2*A1*Contact1+b1;
    nVec1 = nVec1 / norm(nVec1);
    
    nVec2 = 2*A2*Contact2+b2;
    nVec2 = nVec2 / norm(nVec2);
    
    %(2) Calculate the relative velocity at contact point in inertial frame
    % 3*3 Rotation matrix
    Rot1 = reshape(RMatrx(mm,1:3,1:3),[3,3]);
    Rot2 = reshape(RMatrx(nn,1:3,1:3),[3,3]);
    
    % Rotation rate in inertial frame
    rot1 = (Rot1)'* reshape(rot_pf(mm,:),[3,1]); 
    rot2 = (Rot2)'* reshape(rot_pf(nn,:),[3,1]);
    
    % rc: particle center to contact point vector
    rc1 = Contact1 - [rpar(mm,1),rpar(mm,2),rpar(mm,3)]';
    rc2 = Contact2 - [rpar(nn,1),rpar(nn,2),rpar(nn,3)]';
    
    % Velocity at contact point
    vc1 = reshape(vpar(mm,:),[3,1]) + cross(rot1,rc1);
    vc2 = reshape(vpar(nn,:),[3,1]) + cross(rot2,rc2);
    
    vrel = vc1 - vc2;
    
    % Find the tangent direction
    vrel_tangent = vrel - dot(vrel,nVec1)*nVec1;
    if (norm(vrel_tangent) < (1E-10))
    tVec = [1;0;0];
    else
    tVec = vrel_tangent/norm(vrel_tangent);
    end
    
    % Transform Contact point: from inertial frame to particle frame
    RMatrx(mm,4,4) = 1;
    RMatrx(nn,4,4) = 1;
    
    % Tranform contact points to particle frame: x_pf = R*T*x_in
    Contact1_pf = (reshape(RMatrx(mm,:,:),[4,4]))*...
        ([Contact1;1] - [rpar(mm,1),rpar(mm,2),rpar(mm,3),0]');
    Contact2_pf = (reshape(RMatrx(nn,:,:),[4,4]))*...
        ([Contact2;1] - [rpar(nn,1),rpar(nn,2),rpar(nn,3),0]');
    
    % Curvature at contact point
    [Curv1,Curv2] = F_Curvature(Contact1_pf,Contact2_pf,a,b,c);
    
    R12 = 1/(Curv1+Curv2); % Reduced radius
    
    % magnitude of normal elastic force
    Fne_mag = 4/3*Elastic0*sqrt(R12)*(delta_n^(1.5)); 
    
    Fne1 = -Fne_mag*nVec1; %normal elastic force
    Fne2 = -Fne1;
    
    Ft1 = -frict*Fne_mag*tVec; % tangent force
    Ft2 = -Ft1;
    
    % Find the damping coefficient of sliding resistence
    alpha_t = 1.2728 - 4.2783*erest +11.087*erest^2 -22.348*erest^3 + ...
        27.467*erest^4 - 18.022*erest^5 + 4.8218*erest^6;
    k_N = 4/3*Elastic0*sqrt(R12*delta_n);
    damp = alpha_t*sqrt(mass*k_N); % damping coefficient
%     fprintf('Damping coefficient = %f\n',damp);
    Fnd1 = -damp*dot(vrel,nVec1)*nVec1;
    Fnd2 = -damp*dot(-vrel,nVec2)*nVec2;
    
    % total contact force
    F1 = Fne1 + Fnd1 + Ft1;
    F2 = Fne2 + Fnd2 + Ft2;
    
    % total contact torque in inertial frame
    M1 = cross(rc1,F1);
    M2 = cross(rc2,F2);
    
    % Output normal elastic force to file
    fprintf (fid,' %10.6f %10.6f %10.6f %10.6f %6d %6d\n',...
        time,delta_n,dot(vrel,nVec1),norm(vrel_tangent),mm,nn);%,...
        %rpar(mm,1),rpar(nn,1),Contact1(1),Contact1(2),Contact1(3),...
        %Contact2(1),Contact2(2),Contact2(3));
    
    % Print relative velocity on screen
%     fprintf(fid,'Normal Velocity = %f\n',dot(vrel,nVec1));
%     fprintf('Tangent Velocity = %f\n',norm(vrel_tangent));
    
end

