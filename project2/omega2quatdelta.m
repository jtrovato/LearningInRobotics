function q_delta = omega2quatdelta(omega , delta_t)
if norm(omega) == 0
    q_delta = [1 0 0 0]';
    return
end
%convert from rotaional velocity to a delta quaternion
alpha_delta = norm(omega)*delta_t;
e_delta = omega/norm(omega);
q_delta = [cos(alpha_delta/2);e_delta*sin(alpha_delta/2)];


end

