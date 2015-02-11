function q= rot2quat(rot)
if norm(rot) == 0
    q = [1 0 0 0]';
    return
end
alpha = norm(rot);
e = rot/norm(rot);
q = [cos(alpha/2);e*sin(alpha/2)];

end

