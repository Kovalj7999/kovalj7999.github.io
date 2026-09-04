% =========================================================================
% Title:       EEG Source Localization via Regularization and Truncated SVD
% Course:      Numerical Methods for Data Mining
% Author:      Chiara Giavalisco
% Description: Simulation of a 3D spherical head volume conductor with 32 
%              internal current dipoles (sources) and 5 surface electrodes.
%              Solves the underdetermined inverse problem (G*x = V) using:
%                1. Minimum-norm Tikhonov-like pseudo-inversion: G' * inv(G*G')
%                2. Moore-Penrose pseudo-inverse (pinv) with SVD tolerance thresholding.
% =========================================================================

clear; clc; close all

%% 1. Geometric Head Model & Source Generation
% Approximate the head volume conductor as a sphere of radius R = 20
[X,Y,Z] = sphere; 
x=20*X;
y=20*Y;
z=20*Z;

surf(x,y,z)
axis equal

% Define 32 internal neural sources randomly distributed inside the volume
% Vectorize mesh coordinates into single columns
x1=X(:);
y1=Y(:);
z1=Z(:);

% Sample 32 random points and scale coordinates to keep them strictly inside (radius < 19)
r=randi([1,441],32,1);
s=rand(32,1);
x1=x1(r)*19.*s;
y1=y1(r)*19.*s;
z1=z1(r)*19.*s;

% Visualize the 3D volume conductor alongside the internal sources
surf(x,y,z)
alpha(0.2)
hold on 
plot3(x1,y1,z1,'ob')
title('32 sources within the sphere')

%% 2. Electrode Sensor Placement
% Define coordinates for 5 sensors positioned on the sphere surface (radius = 20)
figure
surf(x,y,z)
alpha(0.2)
x2=[0 0 0 17.3 -17.3];
y2=[0 20 -20 0 0];
z2=[20 0 0 -10 -10];

hold on 
plot3(x2(1),y2(1),z2(1),'o','Linewidth',3)
plot3(x2(2),y2(2),z2(2),'o','Linewidth',3)
plot3(x2(3),y2(3),z2(3),'o','Linewidth',3)
plot3(x2(4),y2(4),z2(4),'o','Linewidth',3)
plot3(x2(5),y2(5),z2(5),'o','Linewidth',3)
hold off
title('5 electrodes at the surface of the sphere')

%% 3. Synthetic EEG Time Series Generation (Forward Dynamics)
% Simulate multichannel autoregressive signals: 5 channels x 100 time points
numit=100;
t=4;
V=zeros(5,numit+3); % Pre-allocation: [5 channels x (numit + transient buffer)]

% Initial conditions
V(1,t-3)=1.1;
V(1,t-2)=1.6;
V(1,t-1)=0.95*sqrt(2)*V(1,t-1)-0.9025*V(1,t-2)+rand(1);
V(4,t-1)=1.45; 
V(5,t-1)=0.9;

% Autoregressive signal propagation across channels with additive Gaussian noise
for i=4:numit+3
    w=randn(5,1);
    V(1,i)=0.95*sqrt(2)*V(1,i-1)-0.9025*V(1,i-2);
    V(2,i)=0.5*V(1,i-2);
    V(3,i)=-0.4*V(1,i-3);
    V(4,i)=-0.5*V(1,i-2)+0.25*sqrt(2)*V(4,i)+0.25*sqrt(2)*V(5,i-1);
    V(5,i)=0.25*sqrt(2)*V(4,i-1)+0.25*sqrt(2)*V(5,i-1);
    V(:,i)=V(:,i)+w;
end

% Remove the transient warm-up window to yield 100 actual time samples
V=V(:,4:end);

%% 4. Lead Field (Gain) Matrix Construction
% Build forward operator G (5 sensors x 32 sources) via inverse squared Euclidean distance
for i=1:5
    X2(i,1:32)=x2(i);
    Y2(i,1:32)=y2(i);
    Z2(i,1:32)=z2(i);
    for j=1:32
        X1(1:5,j)=x1(j);
        Y1(1:5,j)=y1(j);
        Z1(1:5,j)=z1(j);
    end
end
R=(X1-X2).^2+(Y1-Y2).^2+(Z1-Z2).^2;
c=1;
G=c./R;

%% 5. Inverse Problem Formulation & Regularization
% Underdetermined system (5 sensors << 32 sources): G*x = V
% Solution 1: Right pseudo-inverse (minimum-norm solution x = G' * (G*G')^-1 * V)
x = G'*(inv(G*G'))*V;

% Solution 2: Truncated SVD via Moore-Penrose pseudo-inverse with varying singular-value tolerances
% Singular values <= tol are set to zero to filter noise sensitivity
xpi1=pinv(G,10e-1)*V;
xpi2=pinv(G,10e-2)*V;
xpi3=pinv(G,10e-3)*V;
xpi4=pinv(G,10e-4)*V;

% Frobenius norm error comparison against the analytic minimum-norm estimate
err3=norm(x-xpi3)
err4=norm(x-xpi4)

%% 6. Reconstruction Quality & Sensor Signal Tracking
t=linspace(1,100);

% --- Channel 1 Reconstruction Comparison ---
figure
tiledlayout(2,2);
nexttile
plot(t,V(1,:),'g') % Ground truth recorded potential
title('Signal 1 in time')
nexttile
plot(t,G(1,:)*x,'b') % Forward re-projection via closed-form minimum-norm
title('Signal 1 in time with regularization')
nexttile
plot(t,G(1,:)*xpi4,'m') % Re-projection with low singular value thresholding (tol = 10e-4)
title('Signal 1 in time with pinv and tol=10e-4')
nexttile
plot(t,G(1,:)*xpi3,'r') % Re-projection with over-regularized thresholding (tol = 10e-3)
title('Signal 1 in time with pinv and tol=10e-3')

% --- Channel 2 Reconstruction Comparison ---
figure
tiledlayout(2,2);
nexttile
plot(t,V(2,:),'g')
title('Signal 2 in time')
nexttile
plot(t,G(2,:)*x,'b')
title('Signal 2 in time with regularization')
nexttile
plot(t,G(2,:)*xpi4,'m')
title('Signal 2 in time with pinv and tol=10e-4')
nexttile
plot(t,G(2,:)*xpi3,'r')
title('Signal 2 in time with pinv and tol=10e-3')

% --- Channel 3 Reconstruction Comparison ---
figure
tiledlayout(2,2);
nexttile
plot(t,V(3,:),'g')
title('Signal 3 in time')
nexttile
plot(t,G(3,:)*x,'b')
title('Signal 3 in time with regularization')
nexttile
plot(t,G(3,:)*xpi4,'m')
title('Signal 3 in time with pinv and tol=10e-4')
nexttile
plot(t,G(3,:)*xpi3,'r')
title('Signal 3 in time with pinv and tol=10e-3')

% --- Channel 4 Reconstruction Comparison ---
figure
tiledlayout(2,2);
nexttile
plot(t,V(4,:),'g')
title('Signal 4 in time')
nexttile
plot(t,G(4,:)*x,'b')
title('Signal 4 in time with regularization')
nexttile
plot(t,G(4,:)*xpi4,'m')
title('Signal 4 in time with pinv and tol=10e-4')
nexttile
plot(t,G(4,:)*xpi3,'r')
title('Signal 4 in time with pinv and tol=10e-3')

% --- Channel 5 Reconstruction Comparison ---
figure
tiledlayout(2,2);
nexttile
plot(t,V(5,:),'g')
title('Signal 5 in time')
nexttile
plot(t,G(5,:)*x,'b')
title('Signal 5 in time with regularization')
nexttile
plot(t,G(5,:)*xpi4,'m')
title('Signal 5 in time with pinv and tol=10e-4')
nexttile
plot(t,G(5,:)*xpi3,'r')
title('Signal 5 in time with pinv and tol=10e-3')

%% 7. Residual Error Surface Plots
% Absolute tracking error over time: |V - G*x| for each sensor channel
err=abs(V-G*x);
figure
contourf(err)
colorbar
title('Error in time with regularization')

errpi4=abs(V-G*xpi4);
figure
contourf(errpi4)
colorbar
title('Error in time with pinv and tol=10e-4')

errpi3=abs(V-G*xpi3);
figure
contourf(errpi3)
colorbar
title('Error in time with pinv and tol=10e-3')
