function q_delta = omega2quatdelta(omega , delta_t)
%convert from rotaional velocity to a delta quaternion
alpha_delta = norm(omega)*delta_t
e_delta1 = omega/norm(omega);
q_delta = [cos(alpha_delta/2);e_delta1*sin(alpha_delta/2)];


end

