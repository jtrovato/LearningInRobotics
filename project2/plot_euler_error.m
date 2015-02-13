[x,y,z] = rot2euler(kf_R);
[xr,yr,zr] = rot2euler(rots);
plot(imu_ts, squeeze(x),'r', imu_ts, squeeze(y), 'g', imu_ts, squeeze(z),'b')
hold on
plot(vicon_ts, squeeze(xr), 'r',  vicon_ts, squeeze(yr), 'g', vicon_ts, squeeze(zr),'b', 'LineWidth', 2)
title('Euler Angles in ukf and vicon')
xlabel('time (s)');
ylabel('rad')
legend('theta_x (ukf)', 'theta_y (ukf)', 'theta_z (ukf)', 'theta_x (vicon)', 'theta_y (vicon)', 'theta_z (vicon)');
grid on