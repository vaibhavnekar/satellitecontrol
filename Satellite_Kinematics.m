function [w_s, q_i_s]= Satellite_Kinematics(I_c, w_init_s, Torque_s, q0_i_s, q_s_c, Ts)
%#codegen
persistent w_c_km1 q_i_c_km1
% Convert torque from satellite to controller frame
Torque_c_temp = quatmultiply(quatmultiply(quatinv(q_s_c), [0; Torque_s]), q_s_c);
Torque_c = Torque_c_temp(2:4,1);
if isempty(w_c_km1) || isempty(q_i_c_km1)
    % First time this function is run, pass through initial conditions
    w_s = w_init_s;
    q_i_s = q0_i_s;
    % Store for next time
    w_c_temp = quatmultiply(quatmultiply(quatinv(q_s_c), [0; w_s]), q_s_c);
    w_c_km1 = w_c_temp(2:4);
    q_i_c_km1 = quatmultiply(q_i_s, q_s_c);
else
    % Second time this function is run, time to integrate
    x_km1 = [q_i_c_km1; w_c_km1];
    % Integrate
    k1 = Kinematics(x_km1, I_c, Torque_c);
    k2 = Kinematics(x_km1 + 0.5*k1*Ts, I_c, Torque_c);
    k3 = Kinematics(x_km1 + 0.5*k2*Ts, I_c, Torque_c);
    k4 = Kinematics(x_km1 + k3*Ts, I_c, Torque_c);
    x_k = x_km1 + (k1 + 2*k2 + 2*k3 + k4)*Ts/6;
    % Store for next time
    q_i_c_km1 = x_k(1:4);
    w_c_km1 = x_k(5:7);
    % Rotate and return
    q_i_s = quatmultiply(q_i_c_km1, quatinv(q_s_c));
    w_s_temp = quatmultiply(quatmultiply(q_s_c, [0; w_c_km1]), quatinv(q_s_c));
    w_s = w_s_temp(2:4);
end
end

function output = skew_matrix(x)
% Returns the skew symmetric matrix of the input vector
output = [0 -x(3) x(2);...
    x(3) 0 -x(1);...
    -x(2) x(1) 0];
end

function results = Kinematics(x, I, Torque_c)
q = x(1:4);
w = x(5:7);
I_Mat = diag(I);
q_dot = 0.5.*[0 -w'; w -skew_matrix(w)]*q;
w_dot = I_Mat^(-1)*(-skew_matrix(w)*(I_Mat*w) + Torque_c);
% w_dot = [0; 0; 0];
% w_dot(1,1) = ( Torque_c(1) - w(2)*w(3)*(I(3,3) - I(2,2)) )/I(1,1);
% w_dot(2,1) = ( Torque_c(2) - w(1)*w(3)*(I(1,1) - I(3,3)) )/I(2,2);
% w_dot(3,1) = ( Torque_c(3) - w(1)*w(2)*(I(2,2) - I(1,1)) )/I(3,3);
results = [q_dot; w_dot];
end

function qres = quatmultiply(q, r)
q = q';
r = r';
% Calculate vector portion of quaternion product
% vec = s1*v2 + s2*v1 + cross(v1,v2)
vec = [q(:,1).*r(:,2) q(:,1).*r(:,3) q(:,1).*r(:,4)] + ...
    [r(:,1).*q(:,2) r(:,1).*q(:,3) r(:,1).*q(:,4)]+...
    [ q(:,3).*r(:,4)-q(:,4).*r(:,3) ...
    q(:,4).*r(:,2)-q(:,2).*r(:,4) ...
    q(:,2).*r(:,3)-q(:,3).*r(:,2)];
% Calculate scalar portion of quaternion product
% scalar = s1*s2 - dot(v1,v2)
scalar = q(:,1).*r(:,1) - q(:,2).*r(:,2) - ...
    q(:,3).*r(:,3) - q(:,4).*r(:,4);
qres = [scalar vec]';
end

function qinv = quatinv(qin)
q_conj = [qin(1); -qin(2:4)];
qinv = q_conj./sqrt(sum(qin.^2));
end
