clc
close all
clear all

% writerObj = VideoWriter('Animation.avi'); % File for saving the animation into avi file
% open(writerObj);

dt = 10; %Simulation time step
T_end = 10*60*60; % Simulation time
T = 0 : dt : T_end; 


X_1 = zeros(6,length(T)); % First Satellite position in Earth-centered inertial reference frame
X_2 = zeros(6,length(T)); % Second Satellite position in Earth-centered inertial reference frame
xsat_diff_12 = zeros(6,length(T)); % Relative state vector
mu = 3.986*10^14; % Gravitational parameter
R_earth = 6.371e6; % Earth radius
q1 = 1e-6; % First Satellite Magnetic Charge

%% Chaser satellite parameters 
Hight = 500e3; % Hight of the satellite orbit
Rad = R_earth + Hight;
incl_1 = 51.7*pi/180; % Inclination
epsilon = 0; %excentricity
phi = 0; % longitude of the ascending node
omega = 0; % argument of the pericenter
t_pi = 0; % pericenter time
t_cur = 0; % current time
approx = 0; % first step for Newton Method
[r, v, D] = orbitalMotionKeplerian(mu, Rad, epsilon, phi, omega, incl_1, t_pi, t_cur, approx);


X_1(1:3,1) = r; 
X_1(4:6,1) = v; % Absolute initial state vector of the satellite

Omega = cross(X_1(1:3,1),X_1(4:6,1))/norm(X_1(1:3,1))^2; % Orbital angular velocity
A = orbital_dcm(X_1(:,1)); % Transition matrix to orbital reference frame

%% Second satellite parameters
C_rel = [0; 1000; 0; 0; 500; 0]; % Linear motion equation parameters

dX = trajectory(norm(Omega),C_rel,0); % Relative vector
X_2(:,1) = [X_1(1:3,1)+A'*dX(1:3)';X_1(4:6,1)+A'*dX(4:6)'+cross(Omega,A'*dX(1:3)')]; % Absolute initial state vector of the second satellite

i=1;
% Relative state vector calculation
xsat_diff_12(:,i) = [A*(X_2(1:3,i) - X_1(1:3,i));A*((X_2(4:6,i) - X_1(4:6,i))-cross(Omega,(X_2(1:3,i) - X_1(1:3,i))))];

Time_length_plot = 0.9*90*60; % time for the length of the plot in animation
scrsz = get(0,'ScreenSize');
fig = figure('Color', [1 1 1],'Position',[1 1 scrsz(3) scrsz(4)])
load('topo.mat','topo','topomap1'); % For Earth plotting
whos topo topomap1 

%% Main cycle
for i = 2:1:length(T) 
    
    C_constants = coord2const(xsat_diff_12(:,i-1), norm(Omega)); % Relative motion constants calculation 
    
    force(:,i) = zeros(3,1); % acting on the satellite force (now it zeros)
    [~,X_new] = ode45(@(t,X) Right_part(t,X,A'*force(:,i)),[0:1:dt],X_1(:,i-1)); % Integration of the satellite motion equations
    X_1(:,i) = X_new(end,:)';
    
    force(:,i)=zeros(3,1); % acting on the satellite force (now it zeros)
    [~,X_new] = ode45(@(t,X) Right_part(t,X,force(:,i)),[0:1:dt],X_2(:,i-1)); % Integration of the second satellite motion equations
    X_2(:,i) = X_new(end,:)';
  
    Omega = cross(X_1(1:3,i),X_1(4:6,i))/norm(X_1(1:3,i))^2; % Orbital angular velocity
    
    B = DipoleMagModel(X_1(1:3,i) , T(i)); % Magnetic field model 
    L_f = LorentzForce(q1,X_1(1:3,i),X_1(4:6,i),B); % Lorentz force input
    
    A = orbital_dcm(X_1(:,i)); % Transition matrix to orbital reference frame
   
    xsat_diff_12(:,i) = [A*(X_2(1:3,i) - X_1(1:3,i));A*((X_2(4:6,i) - X_1(4:6,i))-cross(Omega,(X_2(1:3,i) - X_1(1:3,i))))]; % Relative state vector in orbital reference frame
    
    %% Draw animation
    if rem(i,int16(length(T)/500))==0 %% Annimation is plotted each length(T)/500 steps
        i/length(T)
        
        if T(i)<Time_length_plot
            interval = [1:i];
        else
            i_start = int16(Time_length_plot/dt);
            interval = [i-i_start:i];
        end
        subplot(1,2,1)
        plot3(xsat_diff_12(1,interval),xsat_diff_12(2,interval),xsat_diff_12(3,interval),'b','LineWidth',2)
        hold on        
        plot3(xsat_diff_12(1,i),xsat_diff_12(2,i),xsat_diff_12(3,i),'ob','LineWidth',3)
        hold on 
        plot3(0,0,0,'or','LineWidth',3)
        axis equal
        grid on
        hold off 
        xlabel('x, m')
        ylabel('y, m')
        zlabel('z, m')
        title('Relative trajectory','FontSize',20)
        
        subplot(1,2,2)
        plot3(X_1(1,interval),X_1(2,interval),X_1(3,interval),'b','LineWidth',2)
        hold on
        plot3(X_2(1,interval),X_2(2,interval),X_2(3,interval),'r','LineWidth',2)

        plot3(X_1(1,i),X_1(2,i),X_1(3,i),'ob','LineWidth',3)
        hold on
        plot3(X_2(1,i),X_2(2,i),X_2(3,i),'or','LineWidth',3)
        
        % Earth plotting
        [x,y,z] = sphere(50);
        [x, y, z] = EarthRotation(x, y, z, i*dt);
        props.AmbientStrength = 0.1;
        props.DiffuseStrength = 1;
        props.SpecularColorReflectance = .5;
        props.SpecularExponent = 20;
        props.SpecularStrength = 1;
        props.FaceColor= 'texture';
        props.EdgeColor = 'none';
        props.FaceLighting = 'phong';
        props.Cdata = topo;
        surface(x*R_earth,y*R_earth,z*R_earth,props);
        light('position',[-1 0 1]);
        light('position',[-1.5 0.5 -0.5], 'color', [.6 .2 .2]);
        axis equal
        view(45,45)
        xlabel('X, m')
        ylabel('Y, m')
        zlabel('Z, m')
        legend('1st satellite','2nd satellite','Location','NorthEast')
        grid on
        hold off
        
        pause(0.1)   
%         F = getframe(fig);  % for saving the animation into the avifile
%         writeVideo(writerObj,F);
    end

end
% close(writerObj);