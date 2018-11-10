% Karthik Mysore Srinivasa
% December 8th 2014
%function Torque_s = Controller(w_s,q_i_s)
clear all
clc
I = [27 17 25];% In Satellite Frame[kg-m^2]; Default values are [27 17 25]
I_Mat = diag(I);% Principal Inertia Matrix
eps = 0.001; % Epsilon
kp = 50;% Propotional Constant
kv = 50;% Velocity Constant
Tf = 100000;% Enter Simulation Time seconds
Ts = 0.1;% Step Time in seconds
K = zeros(3,Tf);% Dummy matrix to store Angular Velocity in rad/sec at every instant
G = zeros(4,Tf);% Dummy Matrix to store Quaternions at every instant
U = zeros(3, Tf);% Dummy Matrix to store Control values
m_coils = zeros(3,Tf);% Dummy Matrix to store Residual Dipole Moments in A-m^2 at every instant
Torque_s = zeros(3,Tf);% Dummy Matrix to store Torque in N-m at every instant
Y = zeros(3,Tf);% dummy Matrix to store Euler Angles in Degrees at every instant
q_s_c = [1;0;0;0];% Quaternion Rotation from Satellite to Controller Frame
q_init_s = [1;0;0;0];...[-0.4873; 0.0345; 0.0542; 0.8709];% Initial Quaternions
w_init_s = [0.01; 0.01; 0.01];% Initial Angular Velocity in rad/sec
T = zeros(1, Tf); % Dummy Vector to store Time in seconds as it increments
% %% Compute Magnetic Field vector using World Magnetic Model
for i = 1:Tf
    T(i) = i-1;
    
    if T(i) == 0
        
        
        % References: SGP4_Setup by Brandon Jackson
        % Orbit (TLE)
        longstr1 = '1 25544U 14067A 14342.51579142 .00005418 00000-0 10235-3 0 3046';
        longstr2 = '2 25544 51.6494 243.7352 0003674 255.5105 239.9716 15.50141065841738';% Initialize and Start the SGP4 propagator
        SGP4_Setup(longstr1, longstr2) % Orbit Propagator Setup
        ECEF_Init = sgp4(T(i))*1000; % Position Vector for each Time Instant in ECEF Frame w.r.t ECI Frame
        LLA_Init = ecef2lla(ECEF_Init'); % Transform to Geodetic Frame w.r.t ECI Frame
        Latitude = LLA_Init(1); % Compute Real-Time Latitude
        Longitude = LLA_Init(2); % Compute Real-Time Longitude
        Altitude = LLA_Init(3)/1000; % Compute Real-Time Altitude
        %% Compute Magnetic Field vector using World Magnetic Model
        [b_0_t] = igrfmagm(Altitude, Latitude, Longitude, decyear(2014,12,8))*1E-9; % Convert nT to T;
        NED_Init = dcmecef2ned(Latitude, Longitude)*ECEF_Init; % Direction Cosine Matrix for transformation from ECEF to NED Frame
        %% First time when Attitude and Rate are NOT available
        u = -((eps^2*kp*q_init_s(2:4))+(eps*kv*w_init_s)); % Control Law
        U(:,i) = u;
        A_q = quat2dcm(q_init_s'); % Attitude Matrix
        b_t = (A_q)*b_0_t'; % Magnetic Field in Satellite Body Frame
        S_b_t = [0, b_t(3), -b_t(2);...
            -b_t(3), 0, b_t(1);...
            b_t(2), -b_t(1), 0]; % Skew-Symmetric Matrix
        Transpose = S_b_t';
        
        Abs = (norm(b_0_t)^2)^(-1);
        m = Abs*Transpose*u;
        m_coils(:, i) = m; % Residual Dipole Moment in A-m^2
        T_coils = cross(m_coils(:, i), b_t); % Control Torque in N-m
        %T_coils = S_b_t*m_coils(:, i);
        Torque_s(:, i) = T_coils; % Control Torque in N-m
        [w_s, q_i_s]= Satellite_Kinematics(I, w_init_s...
            , Torque_s(:, i), q_init_s...
            , q_s_c...
            , Ts); % Angular Velocity and Quaternions at the next time instant
        G(:, 1) = q_init_s; % Quaternions at the first instant
        G(:, i+1) = q_i_s(:);
        K(:, i+1) = w_s; % Angular Velocity at the first instant
        K(:, 1) = w_init_s;
        [Y(1,i), Y(2,i), Y(3,i)] = quat2angle(G(:,i)', 'ZXZ'); % Euler Angles in Radians
        Y(1,i) = Y(1,i)*(180/pi);Y(2,i) = Y(2,i)*(180/pi);Y(3,i) = Y(3,i)*(180/pi); % Euler Angles in Degrees
        
    else
        % References: SGP4_Setup by Brandon Jackson
        % Orbit (TLE)
        longstr1 = '1 25544U 14067A 14342.51579142 .00005418 00000-0 10235-3 0 3046';
        longstr2 = '2 25544 51.6494 243.7352 0003674 255.5105 239.9716 15.50141065841738';
        % Initialize and Start the SGP4 propagator
        % References: SGP4_Setup by Brandon Jackson
        SGP4_Setup(longstr1, longstr2) % Orbit Propagator Setup
        ECEF_Init = sgp4(T(i))*1000; % Position Vector for each Time Instant in ECEF Frame w.r.t ECI Frame
        LLA_Init = ecef2lla(ECEF_Init'); % Transform to Geodetic Frame w.r.t ECI Frame
        Latitude = LLA_Init(1); % Compute Real-Time Latitude
        Longitude = LLA_Init(2); % Compute Real-Time Longitude
        if LLA_Init(3)/1000 < 600000
            Altitude = LLA_Init(3)/1000; % Compute Real-Time Altitude
        else
            Altitude = 600000;
        end
        %% Compute Magnetic Field vector using World Magnetic Model
        [b_0_t] = igrfmagm(Altitude, Latitude, Longitude, decyear(2014,12,8))*1E-9; % Convert nT to T;
        NED_Init = dcmecef2ned(Latitude, Longitude)*ECEF_Init;
        %% If Attitude and Rate Feedback available, compute Control Input'u'
        u = -((eps^2*kp*((G((2:4), (i)))))+(eps*kv*((K(:, (i)))))); %Control Law
        A_q = quat2dcm(G(:,i)'); % Attitude Matrix
        b_t = (A_q)*b_0_t'; % Magnetic Field in Satellite Body Frame
        S_b_t = [0, b_t(3), -b_t(2);...
            -b_t(3), 0, b_t(1);...
            b_t(2), -b_t(1), 0]; % Skew-Symmetric Matrix
        Transpose = S_b_t';
        U(:,i) = u;
        Abs = (norm(b_0_t)^2)^(-1);
        m = Abs*Transpose*u;
        m_coils(:, i) = m; % Residual Dipole Moment in A-m^2
        T_coils = cross(m_coils(:, i), b_t); % Control Torque in N-m
        %T_coils = S_b_t*m_coils(:, i);
        Torque_s(:, i) = T_coils; % Control Torque in N-m
        [w_s, q_i_s]= Satellite_Kinematics(I, w_init_s...
            , Torque_s(:, i), q_init_s...
            , q_s_c...
            , Ts); % Angular Velocity and Quaternions at the next instant
        G(:, i+1) = q_i_s(:);
        K(:, i+1) = w_s;
        [Y(1,i), Y(2,i), Y(3,i)] = quat2angle(G(:,i)', 'ZXZ'); % Euler Angles in Radians
        Y(1,i) = Y(1,i)*(180/pi);Y(2,i) = Y(2,i)*(180/pi);Y(3,i) = Y(3,i)*(180/pi); % Euler Angles in Degrees
    end
end
p = 0:(Tf);
figure(1);
hold on
plot(p, K(1,:),'r');
plot(p, K(2,:),'b');
plot(p, K(3,:),'g');
xlabel('Time in seconds');
ylabel('Angular Velocity in radians per second');
title('Angular Velocity in radians per second');
legend('w_1','w_2','w_3');
set(gca,'XLim',[0 Tf]);
hold off
figure(2);
hold on
plot(p, G(1,:),'r');
plot(p, G(2,:),'b');
plot(p, G(3,:),'g');
plot(p, G(4,:),'m');
xlabel('Time in seconds');
ylabel('Quaternions');
title('Quaternions');
legend('q1','q2','q3','q4');
set(gca,'XLim',[0 Tf]);
hold off
figure(3);
hold on
plot(T, Torque_s(1,:),'r');
plot(T, Torque_s(2,:),'b');
plot(T, Torque_s(3,:),'g');
xlabel('Time in seconds');
ylabel('Torque in Newton-meter');
title('Torque in Newton-meter');
legend('Torque1','Torque2','Torque3');
hold off
figure(4);
hold on
plot(T, m_coils(1,:),'r');
plot(T, m_coils(2,:),'b');
plot(T, m_coils(3,:),'g');
xlabel('Time in seconds');
ylabel('Dipole Moment in Ampere - meters^2');
title('Dipole Moment in Ampere - meters^2');
legend('m1','m2','m3');
hold off
figure(5);
hold on
plot(T, Y(1,:),'r');
plot(T, Y(2,:),'b');
plot(T, Y(3,:),'g');
xlabel('Time in seconds');
ylabel('Euler Angles in Degrees');
title('Euler Angles in Degrees');
legend('theta1','theta2','theta3');
hold off