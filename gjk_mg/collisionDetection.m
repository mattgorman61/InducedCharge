% GJK Collision Detection
clc; clear all; close all;

% https://cse442-17f.github.io/Gilbert-Johnson-Keerthi-Distance-Algorithm/
% https://www.haroldserrano.com/blog/visualizing-the-gjk-collision-algorithm
% https://blog.winter.dev/2020/gjk-algorithm/

% https://blog.winter.dev/2020/epa-algorithm/

%% SHAPE 1
theta = linspace(0,2*pi,20);
phi = linspace(0,pi,20);
numpts = length(theta)*length(phi);

sx1 = zeros(numpts,1);
sy1 = zeros(numpts,1);
sz1 = zeros(numpts,1);

for t=1:length(theta)
    for p=1:length(phi)
        th = theta(t);
        ph = phi(p);
        
        sx1((t-1)*length(theta) + p) = cos(th)*sin(ph);
        sy1((t-1)*length(theta) + p) = sin(th)*sin(ph);
        sz1((t-1)*length(theta) + p) = cos(ph);
    end
end

sx1 = 2 + sx1;
sy1 = 2 + sy1;
sz1 = 3 + sz1;

s1pts = [sx1,sy1,sz1];


%% SHAPE 2
% Sphere
theta = linspace(0,2*pi,20);
phi = linspace(0,pi,20);
numpts = length(theta)*length(phi);

sx2 = zeros(numpts,1);
sy2 = zeros(numpts,1);
sz2 = zeros(numpts,1);

for t=1:length(theta)
    for p=1:length(phi)
        th = theta(t);
        ph = phi(p);
        
        sx2((t-1)*length(theta) + p) = cos(th)*sin(ph);
        sy2((t-1)*length(theta) + p) = sin(th)*sin(ph);
        sz2((t-1)*length(theta) + p) = cos(ph);
    end
end

sx2 = 2.5 + sx2;
sy2 = 5 + sy2;
sz2 = 3.5 + sz2;

s2pts = [sx2,sy2,sz2];


%% DISPLAY SHAPES
figure;
hold on; grid on; box on;
axis equal;
view(35,20);

S = [1.4215    0.3086   -3.4910
   -1.2200    0.5675   -3.1652
    1.9732    2.3236   15.9458
   -0.2015   -1.2589   -4.9727];

scatter3(sx1,sy1,sz1,'b','filled');
scatter3(sx2,sy2,sz2,'r','filled');

% Shape Center points
cs1 = [mean(sx1),mean(sy1),mean(sz1)];
cs2 = [mean(sx2),mean(sy2),mean(sz2)];

scatter3(cs1(1), cs1(2), cs1(3),'b','filled' );
scatter3(cs2(1), cs2(2), cs2(3),'r','filled' );
% 
% figure(1000);
% hold on; box on; grid on;
% scatter3(1.4215, 0.3086, -3.4910, 'k', 'filled');
% scatter3(-1.2200, 0.5675, -3.1652, 'g', 'filled');
% scatter3(1.9732, 2.3236, 15.9458, 'm', 'filled');
% scatter3(-0.2015, -1.2589, -4.9727, 'c', 'filled');
% scatter3(0,0,0,'r','filled');

%% GJK ALGORITHM

% Original direction is from center of one shape to center of the other shape
initdir_arb = [cs2(1)-cs1(1), cs2(2)-cs1(2), cs2(3)-cs1(3)];

% Create Simplex. Max dimension of Simplex S is 4 (tetrahedron case)
[x_simp1,y_simp1,z_simp1] = furthestPtfunct(s1pts, initdir_arb);
S = [x_simp1,y_simp1,z_simp1];

dir = -S(1,:);

while(true)

    [x_supp,y_supp,z_supp] = supportfunct(s1pts,s2pts,dir);
    supp = [x_supp,y_supp,z_supp];
    
    if(dot(supp,dir)<=0)
        fprintf('\n\n NO COLLISION \n\n');
        break;
    end
    
    S = [supp; S]; %#ok
    
    [bool, S, dir] = nextsimplexfunct(S,dir);
    if( bool == 1)
        fprintf('\n\n COLLISION \n\n');
        break;
    end
     
end



 



