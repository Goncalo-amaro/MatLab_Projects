%% Two charged Satellites system - Draft 2
clc
clear all
close all

% writerObj = VideoWriter('Animation.avi'); % File for saving the animation into avi file
% open(writerObj);

dt    = 10; %Simulation time step
T_end = 1.5*24*60*60; % Simulation time
T     = 0 : dt : T_end;

const_range = [0.25, 0.5, 0.75, 1];
number_simu = 2;

for inter = 1:length(const_range) - 1 
   
    n1 = random('unif',-0.1,0.1);
    m1 = n1 / abs(n1);
    val_C1 = random('unif', const_range(inter), const_range(inter + 1));
    inic_C11(inter) = m1 * val_C1;
    
    n2 = random('unif',-0.1,0.1);
    m2 = n2 / abs(n2);
    val_C2 = random('unif', const_range(inter), const_range(inter + 1));
    inic_C12(inter) = m2 * val_C2;
    
end

conv_time = zeros(number_simu, length(const_range) - 1);

for num_par = 1:length(const_range) - 1
    for s = 1:number_simu
        
        X_des   = zeros(6,length(T)); % Imaginary desired trajectory for the system (x_sat_12) to achieve - in ECRF
        X_1     = zeros(6,length(T)); % First Satellite position in Earth-centered inertial reference frame (ECRF)
        X_2     = zeros(6,length(T)); % Second Satellite position in Earth-centered inertial reference frame (ECRF)
        mu      = 3.986*10^14; % Gravitational parameter
        R_earth = 6.371e6; % Earth radius
        
        %% Chaser satellite parameters
        Hight   = 500e3; % Hight of the satellite orbit
        Rad     = R_earth + Hight;
        incl_1  = 51.7*pi/180; % Inclination
        epsilon = 0; %excentricity
        phi     = 0; % longitude of the ascending node
        omega   = 0; % argument of the pericenter
        t_pi    = 0; % pericenter time
        t_cur   = 0; % current time
        approx  = 0; % first step for Newton Method
        
        [r, v, D] = orbitalMotionKeplerian(mu, Rad, epsilon, phi, omega, incl_1, t_pi, t_cur, approx);
        
        X_des(1:3,1) = r;
        X_des(4:6,1) = v; % Absolute initial state vector of the satellite
        
        Omega = cross(X_des(1:3,1),X_des(4:6,1))/norm(X_des(1:3,1))^2; % Orbital angular velocity
        A     = orbital_dcm(X_des(:,1)); % Transition matrix to orbital reference frame
        
        %% First satellite parameters
        C1_g1 = inic_C11(num_par);
        C3_g1 = random('unif',-5, 5);
        C_rel_1 = [C1_g1; random('unif', -5, 5); random('unif', -5, 5); C3_g1; random('unif', -5, 5); random('unif', -5, 5)];
        
        dX_1     = trajectory(norm(Omega),C_rel_1,0); % Relative vector
        X_1(:,1) = [ X_des(1:3,1) + A'*dX_1(1:3)'; X_des(4:6,1) + A'*dX_1(4:6)' + cross(Omega,A'*dX_1(1:3)') ]; % Absolute initial state vector of the second satellite
        
        
        %% Second satellite parameters
        C1_g2 = inic_C12(num_par);
        C3_g2 = random('unif',-5, 5);
        C_rel_2 = [C1_g2; random('unif', -5, 5); random('unif', -5, 5); C3_g2; random('unif', -5, 5); random('unif', -5, 5)];
        
        dX_2     = trajectory(norm(Omega),C_rel_2,0); % Relative vector
        X_2(:,1) = [ X_des(1:3,1) + A'*dX_2(1:3)'; X_des(4:6,1) + A'*dX_2(4:6)' + cross(Omega,A'*dX_2(1:3)') ]; % Absolute initial state vector of the second satellite
        
        
        %% Relative state vector calculation
        xsat_diff_1(:,1) = [A*(X_1(1:3,1) -X_des(1:3,1)); A*((X_1(4:6,1) -X_des(4:6,1)) - cross(Omega,(X_1(1:3,1) -X_des(1:3,1))))];
        xsat_diff_2(:,1) = [A*(X_2(1:3,1) -X_des(1:3,1)); A*((X_2(4:6,1) -X_des(4:6,1)) - cross(Omega,(X_2(1:3,1) -X_des(1:3,1))))];
        xsat_diff(:,1)   = [A*(X_2(1:3,1) -X_1(1:3,1)); A*((X_2(4:6,1) -X_1(4:6,1)) - cross(Omega,(X_2(1:3,1) -X_1(1:3,1))))];
        
        
        %% constants for the intial conditions
        i = 1;
        C_des  = [0; 5; 5; 5; 5; 5];
        B_2d = sqrt(C_des(2)^2 + C_des(3)^2);
        B_3d = - 3 * omega * C_des(1) * T(i) + C_des(4);
        B_4d = sqrt(C_des(5)^2 + C_des(6)^2);

        Omega    = cross(X_des(1:3,1),X_des(4:6,1))/norm(X_des(1:3,1))^2; % Orbital angular velocity
        m        = 1;
        q_max    = 10e-6;
        q_mean_1 = 0;
        q_mean_2 = 0;
        pkk      = 10e-5;
        
        %% Main cycle
        for i = 2:1:length(T)
            
            i/length(T)
            s
            num_par
            
            omega_1 = norm(Omega);
            
            k_a = 1e-5;
            k_b = 1e-5;
            k_x = 1e-8;
            k_z = 1e-8;
            
            k_y = 1e-5;
            
            Bm_1     = DipoleMagModel(X_1(1:3,i-1) , T(i-1)); % Magnetic field model for sat 1
            Bm_2     = DipoleMagModel(X_2(1:3,i-1) , T(i-1)); % Magnetic field model for sat 2
            
            %% constants sat 1
            C_constants1 = coord2const(xsat_diff_1(:,i-1), norm(Omega)); % Relative motion constants calculation
            
            B_11(i) = C_constants1(1);
            B_21(i) = sqrt(C_constants1(2)^2 + C_constants1(3)^2);
            B_31(i) = -3 * C_constants1(1) * omega_1 * T(i) + C_constants1(4);
            B_41(i) = sqrt(C_constants1(5)^2 + C_constants1(6)^2);
            
            delta_B31(i) = (B_31(i) - B_3d);
            
            a1_sinb1(i) = C_constants1(3) / (sqrt(C_constants1(2)^2 + C_constants1(3)^2));
            a1_cosb1(i) = C_constants1(2) / (sqrt(C_constants1(2)^2 + C_constants1(3)^2));
            
            a2_sinb1(i) = - C_constants1(6) / (sqrt(C_constants1(5)^2 + C_constants1(6)^2));
            a2_cosb1(i) =   C_constants1(5) / (sqrt(C_constants1(5)^2 + C_constants1(6)^2));
            
            req_y1 = -k_y * ((B_41(i) -B_4d) * sin(a2_sinb1(i))) * (1/omega_1);
            
            if abs(B_11(i)) <= 0.05 && abs(delta_B31(i)) <= 2.5
                
                req_x1 = -k_x * (omega_1) * ( -B_11(i) -2 * ((B_21(i) -B_2d) * sin(omega_1 * T(i) + a1_sinb1(i))));
                req_z1 = -k_z * ((B_21(i) -B_2d)*sin(omega_1 * T(i) + a1_cosb1(i))  - 2 * delta_B31(i)); % ERROR CHECKING!!!!!!
                
            else
                req_x1 = -k_a * (1/omega_1) * ( -B_11(i));
                req_z1 = 0.5 * ( -3 * B_11(i) * omega_1^2 - k_b * omega_1 * delta_B31(i)); % ERROR CHECKING!!!!!!
  
            end
            
            %% constants sat 2
            C_constants2 = coord2const(xsat_diff_2(:,i-1), norm(Omega)); % Relative motion constants calculation
            
            B_12(i) = C_constants2(1);
            B_22(i) = sqrt(C_constants2(2)^2 + C_constants2(3)^2);
            B_32(i) = -3 * C_constants2(1) * omega_1 * T(i) + C_constants2(4);
            B_42(i) = sqrt(C_constants2(5)^2 + C_constants2(6)^2);
            
            delta_B32(i) = (B_3d - B_32(i));
            
            a1_sinb2(i) = C_constants2(3) / (sqrt(C_constants2(2)^2 + C_constants2(3)^2));
            a1_cosb2(i) = C_constants2(2) / (sqrt(C_constants2(2)^2 + C_constants2(3)^2));
            
            a2_sinb2(i) = - C_constants2(6) / (sqrt(C_constants2(5)^2 + C_constants2(6)^2));
            a2_cosb2(i) =   C_constants2(5) / (sqrt(C_constants2(5)^2 + C_constants2(6)^2));
            
            req_y2 = -k_y * ((B_42(i) - B_4d) * sin(a2_sinb2(i))) * (1/omega_1);
            
            if abs(B_12(i)) <= 0.05 && abs(delta_B32(i)) <= 2.5
                req_x2 = -k_a * (1/omega_1) * ( -B_12(i) -2 * ((B_22(i) -B_2d) * sin(omega_1 * T(i) + a1_sinb2(i))));
                req_z2 = 0.5 * ( -3 * B_12(i) * omega_1^2 + k_b * omega_1 * B_32(i));
            else
                req_x2 = -k_x * (1/omega_1) * ( -B_12(i));
                req_z2 = -k_z * (1/omega_1) * ( -3 * B_12(i) * delta_B32(i) * omega_1 -(2/omega_1) * delta_B32(i));
            end
            
            B_1dif(i) = B_12(i) -B_11(i);
            B_2dif(i) = B_22(i) -B_21(i);
            B_3dif(i) = B_32(i) -B_31(i);
            B_4dif(i) = B_42(i) -B_41(i);
            
            
            %% charge calculation for sat 2 -main cycle-
            L_f_2             = LorentzForce(1,m,X_2(1:3,i-1),X_2(4:6,i-1),Bm_2); % try previous q_mean instead of q=1
            L_f_orbital_req_2 = A * L_f_2;
            % L_f_val_2(:,i)    = L_f_2;
            
            q1_2 = - req_x2 / L_f_orbital_req_2(1);
            q2_2 = - req_y2 / L_f_orbital_req_2(2);
            q3_2 = - req_z2 / L_f_orbital_req_2(3); % take another look to the values - later
            
            % req_force_2(i,:) = [req_x,req_y,req_z]';
            
            if abs(q1_2) > q_max
                q1_2 = q_max * sign(q1_2);
            end
            
            if abs(q3_2) > q_max
                q3_2 = q_max * sign(q3_2);
            end
            
            if abs(q2_2) > q_max
                q2_2 = q_max * sign(q2_2);
            end
            
            q_1_2(i) = q1_2;
            q_2_2(i) = q2_2;
            q_3_2(i) = q3_2;
            
            % q_val_2(i)  = q1_2 + q3_2;
            q_sign_2(i) = sign(q1_2 + q3_2 + q2_2);
            q_general_2 = [q_1_2(i); q_2_2(i); q_3_2(i)];
            q_mean_2(i) = norm(q_general_2) * q_sign_2(i);
            
            if abs(q_mean_2(i)) > q_max
                q_mean_2(i) = q_max * sign(q_mean_2(i));
            end
            
            
            %% charge calculation for sat 1 -main cycle-
            L_f_1             = LorentzForce(1,m,X_1(1:3,i-1),X_1(4:6,i-1),Bm_1); % try previous q_mean instead of q=1
            L_f_orbital_req_1 = A * L_f_1;
            % L_f_val(:,i)    = L_f_1;
            
            q1_1 = - req_x1 / L_f_orbital_req_1(1);
            q2_1 = - req_y1 / L_f_orbital_req_1(2);
            q3_1 = - req_z1 / L_f_orbital_req_1(3); % take another look to the values - later
            
            % req_force_1(i,:) = [req_x,req_y,req_z]';
            
            if abs(q1_1) > q_max
                q1_1 = q_max * sign(q1_1);
            end
            
            if abs(q3_1) > q_max
                q3_1 = q_max * sign(q3_1);
            end
            
            if abs(q2_1) > q_max
                q2_1 = q_max * sign(q2_1);
            end
            
            q_1_1(i) = q1_1;
            q_2_1(i) = q2_1;
            q_3_1(i) = q3_1;
            
            %q_val_1(i)  = q1_1 + q3_1;
            q_sign_1(i) = sign(q1_1 + q3_1 + q2_1);
            q_general_1 = [q_1_1(i); q_2_1(i); q_3_1(i)];
            q_mean_1(i) = norm(q_general_1) * q_sign_1(i);
            
            if abs(q_mean_1(i)) > q_max
                q_mean_1(i) = q_max * sign(q_mean_1(i));
            end
            
            L_f_1              = LorentzForce(q_mean_1(i),m,X_1(1:3,i-1),X_1(4:6,i-1),Bm_1);
            L_f_orbital_1(i,:) = (A * L_f_1);
            
            L_f_2              = LorentzForce(q_mean_2(i),m,X_2(1:3,i-1),X_2(4:6,i-1),Bm_2);
            L_f_orbital_2(i,:) = (A * L_f_2);
            
            
            %% forces implementation -main cycle-
            force(:,i) = zeros(3,1);
            [~,X_new]  = ode45(@(t,X) Right_part(t,X,force(:,i)),[0:1:dt],X_des(:,i-1)); % Integration of the second satellite motion equations
            X_des(:,i)     = X_new(end,:)';
            
            force_2(:,i) = L_f_orbital_2(i,:); % force acting on the satellite 2
            [~,X_new_2]  = ode45(@(t,X) Right_part(t,X,A'*L_f_orbital_2(i,:)'/m),[0:1:dt],X_2(:,i-1)); % Integration of the second satellite motion equations
            X_2(:,i)     = X_new_2(end,:)';
            
            force_1(:,i) = L_f_orbital_1(i,:); % force acting on the satellite 1
            [~,X_new_1]  = ode45(@(t,X) Right_part(t,X,A'*L_f_orbital_1(i,:)'/m),[0:1:dt],X_1(:,i-1)); % Integration of the first satellite motion equations
            X_1(:,i)     = X_new_1(end,:)';
            
            Omega = cross(X_des(1:3,i),X_des(4:6,i))/norm(X_des(1:3,i))^2; % Orbital angular velocity
            
            B     = DipoleMagModel(X_1(1:3,i) , T(i)); % Magnetic field model
            A     = orbital_dcm(X_des(:,i)); % Transition matrix to orbital reference frame
            
            xsat_diff_1(:,i) = [ A*(X_1(1:3,i) - X_des(1:3,i)); A*((X_1(4:6,i) - X_des(4:6,i)) - cross(Omega,(X_1(1:3,i) - X_des(1:3,i)))) ]; % Relative state vector in orbital reference frame
            xsat_diff_2(:,i) = [ A*(X_2(1:3,i) - X_des(1:3,i)); A*((X_2(4:6,i) - X_des(4:6,i)) - cross(Omega,(X_2(1:3,i) - X_des(1:3,i)))) ]; % Relative state vector in orbital reference frame
            xsat_diff(:,i) = [ A*(X_2(1:3,i) - X_1(1:3,i)); A*((X_2(4:6,i) - X_1(4:6,i)) - cross(Omega,(X_2(1:3,i) - X_1(1:3,i)))) ]; % Relative state vector in orbital reference frame
            
            %% time record of the convergence
            if abs(B_11(i)) < 0.25 && abs(B_12(i)) < 0.25
                P(1,i) = 1;
            else
                P(1,i) = 0;
            end
            
            if abs(B_21(i) - B_2d) < 30 && abs(B_22(i) - B_2d) < 30
                P(2,i) = 1;
            else
                P(2,i) = 0;
            end
            
            if abs(delta_B31(i)) < 30 && abs(delta_B32(i)) < 30
                P(3,i) = 1;
            else
                P(3,i) = 0;
            end
            
            if abs(B_41(i) - B_4d) < 30 && abs(B_42(i) - B_4d) < 30
                P(4,i) = 1;
            else
                P(4,i) = 0;
            end
            
            if sum( P(:,i)) == 4 && conv_time(s, num_par) == 0
                conv_time(s, num_par) = i;
            end
        end
        
        figure
        plot3(xsat_diff_1(1,:),xsat_diff_1(3,:),xsat_diff_1(2,:),'b',xsat_diff_2(1,:),xsat_diff_2(3,:),xsat_diff_2(2,:),'r'); % xsat_diff(1,:),xsat_diff(3,:),xsat_diff(2,:),'k'
        grid on
        xlabel('x')
        ylabel('z')
        zlabel('y')
        title('sat1, sat2')
        legend('sat1','sat2')
        
        figure
        plot(T, B_11(:),T, B_12(:),'r')
        hold on
        title('B1')
        legend('sat1','sat2')

        figure
        plot(T, B_21(:),T, B_22(:),'r')
        hold on
        title('B2')
        legend('sat1','sat2')
        
        figure
        plot(T, B_31(:),T, B_32(:),'r')
        hold on
        title('B3')
        legend('sat1','sat2')
        
        figure
        plot(T, B_41(:),T, B_42(:),'r')
        hold on
        title('B4')
        legend('sat1','sat2')
        
        if conv_time(s, num_par) ~= 0
            
            size_v1 = length(B_11(conv_time(s, num_par):end));
            size_v2 = length(B_12(conv_time(s, num_par):end));
            
            constants_error_1{num_par}(:,s) = [abs(sum(B_11(conv_time(s, num_par):end))) / size_v1;
                                               abs(sum(B_21(conv_time(s, num_par):end) - B_2d)) / size_v1;
                                               abs(sum(delta_B31(conv_time(s, num_par):end))) / size_v1;
                                               abs(sum(B_41(conv_time(s, num_par):end) - B_4d)) / size_v1];
            
            constants_var_1{num_par}(:,s) = [sum(B_11(conv_time(s, num_par):end).^2) / size_v1;
                                             sum((B_21(conv_time(s, num_par):end) - B_2d).^2) / size_v1;
                                             sum(delta_B31(conv_time(s, num_par):end).^2) / size_v1;
                                             sum((B_41(conv_time(s, num_par):end) - B_4d).^2) / size_v1];
            
            constants_error_2{num_par}(:,s) = [abs(sum(B_12(conv_time(s, num_par):end))) / size_v2;
                                               abs(sum(B_22(conv_time(s, num_par):end) - B_2d)) / size_v2;
                                               abs(sum(delta_B32(conv_time(s, num_par):end))) / size_v2;
                                               abs(sum(B_42(conv_time(s, num_par):end) - B_4d)) / size_v2];
            
            constants_var_2{num_par}(:,s) = [sum(B_12(conv_time(s, num_par):end).^2) / size_v2;
                                             sum((B_22(conv_time(s, num_par):end) - B_2d).^2) / size_v2;
                                             sum(delta_B32(conv_time(s, num_par):end).^2) / size_v2;
                                             sum((B_42(conv_time(s, num_par):end) - B_4d).^2) / size_v2];
        else
            
            constants_error_1{num_par}(:,s) = zeros(4,1) * NaN;
            
            constants_var_1{num_par}(:,s) = zeros(4,1) * NaN;
            
            constants_error_2{num_par}(:,s) = zeros(4,1) * NaN;
            
            constants_var_2{num_par}(:,s) = zeros(4,1) * NaN;
        end
    end
end


%% figures

for i = 1:num_par
    for j = 1:s
        
        B1_err1(j,i) = constants_error_1{i}(1,j);
        B2_err1(j,i) = constants_error_1{i}(2,j);
        B3_err1(j,i) = constants_error_1{i}(3,j);
        B4_err1(j,i) = constants_error_1{i}(4,j);
        
        B1_var1(j,i) = constants_var_1{i}(1,j);
        B2_var1(j,i) = constants_var_1{i}(2,j);
        B3_var1(j,i) = constants_var_1{i}(3,j);
        B4_var1(j,i) = constants_var_1{i}(4,j);
        
        B1_err2(j,i) = constants_error_2{i}(1,j);
        B2_err2(j,i) = constants_error_2{i}(2,j);
        B3_err2(j,i) = constants_error_2{i}(3,j);
        B4_err2(j,i) = constants_error_2{i}(4,j);
        
        B1_var2(j,i) = constants_var_2{i}(1,j);
        B2_var2(j,i) = constants_var_2{i}(2,j);
        B3_var2(j,i) = constants_var_2{i}(3,j);
        B4_var2(j,i) = constants_var_2{i}(4,j);
        
    end
end

%% code efectiveness check - Sat 1 -
% figure
% boxplot(B1_err1(:,:), const_range(2:end))
% title('Error deviation B_1  value for sat 1')
% ylabel('Meters (m)')
% xlabel('Initial Constant 1 deviation (m)')
% 
% figure
% boxplot(B2_err1(:,:), const_range(2:end))
% title('Error deviation B_2  value for sat 1')
% ylabel('Meters (m)')
% xlabel('Initial constant 1 deviation (m)')
% 
% figure
% boxplot(B3_err1(:,:), const_range(2:end))
% title('Error deviation B_3  value for sat 1')
% ylabel('Meters (m)')
% xlabel('Initial Constant 1 deviation (m)')
% 
% figure
% boxplot(B4_err1(:,:), const_range(2:end))
% title('Error deviation B_4  value for sat 1')
% ylabel('Meters (m)')
% xlabel('Initial Constant 1 deviation (m)')
% 
% figure
% boxplot(B1_var1(:,:), const_range(2:end))
% title('Variancy deviation B_1  value for sat 1')
% ylabel('?')
% xlabel('Initial Constant 1 deviation (m)')
% 
% figure
% boxplot(B2_var1(:,:), const_range(2:end))
% title('Variancy deviation B_2  value for sat 1')
% ylabel('?')
% xlabel('Initial constant 1 deviation (m)')
% 
% figure
% boxplot(B3_var1(:,:), const_range(2:end))
% title('Variancy deviation B_3  value for sat 1')
% ylabel('?')
% xlabel('Initial Constant 1 deviation (m)')
% 
% figure
% boxplot(B4_var1(:,:), const_range(2:end))
% title('Variancy deviation B_4  value for sat 1')
% ylabel('?')
% xlabel('Initial Constant 1 deviation (m)')

%% code efectiveness check - Sat 2 -
% figure
% boxplot(B1_err2(:,:), const_range(2:end))
% title('Error deviation B_1 value for sat 2')
% ylabel('Meters (m)')
% xlabel('Initial Constant 1 deviation (m)')
% 
% figure
% boxplot(B2_err2(:,:), const_range(2:end))
% title('Error deviation B_2 value for sat 2')
% ylabel('Meters (m)')
% xlabel('Initial constant 1 deviation (m)')
% 
% figure
% boxplot(B3_err2(:,:), const_range(2:end))
% title('Error deviation B_3 value for sat 2')
% ylabel('Meters (m)')
% xlabel('Initial Constant 1 deviation (m)')
% 
% figure
% boxplot(B4_err2(:,:), const_range(2:end))
% title('Error deviation B_4 value for sat 2')
% ylabel('Meters (m)')
% xlabel('Initial Constant 1 deviation (m)')
% 
% figure
% boxplot(B1_var2(:,:), const_range(2:end))
% title('Variancy deviation B_1 value for sat 2')
% ylabel('?')
% xlabel('Initial Constant 1 deviation (m)')
% 
% figure
% boxplot(B2_var2(:,:), const_range(2:end))
% title('Variancy deviation B_2 value for sat 2')
% ylabel('?')
% xlabel('Initial constant 1 deviation (m)')
% 
% figure
% boxplot(B3_var2(:,:), const_range(2:end))
% title('Variancy deviation B_3 value for sat 2')
% ylabel('?')
% xlabel('Initial Constant 1 deviation (m)')
% 
% figure
% boxplot(B4_var2(:,:), const_range(2:end))
% title('Variancy deviation B_4 value for sat 2')
% ylabel('?')
% xlabel('Initial Constant 1 deviation (m)')


