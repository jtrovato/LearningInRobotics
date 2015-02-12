%new plot
function hand = newrotplot2(R_kf,R_vicon, hand)
    x= [2 0 0]';
    y= [0 2 0]';
    z= [0 0 2]';
    xr = R_kf*x;
    yr = R_kf*y;
    zr = R_kf*z;
    xr2 = R_vicon*x;
    yr2 = R_vicon*y;
    zr2 = R_vicon*z;
    
    if(isempty(hand))
        subplot(1,2,1)
        hand(1) = plot3([0, xr(1)],[0,xr(2)],[0,xr(3)], 'r', 'LineWidth', 2);
        hold on;
        hand(2) = plot3([0, yr(1)],[0,yr(2)],[0,yr(3)], 'g', 'LineWidth', 2);
        hand(3) = plot3([0, zr(1)],[0,zr(2)],[0,zr(3)], 'b', 'LineWidth', 2);
        axis equal;
        grid on;
        axis([-2 2 -2 2 -2 2]);
        title('kalman');
        subplot(1,2,2)
        hand(4) = plot3([0, xr2(1)],[0,xr2(2)],[0,xr2(3)], 'r', 'LineWidth', 2);
        hold on;
        hand(5) = plot3([0, yr2(1)],[0,yr2(2)],[0,yr2(3)], 'g', 'LineWidth', 2);
        hand(6) = plot3([0, zr2(1)],[0,zr2(2)],[0,zr2(3)], 'b', 'LineWidth', 2);
        hold off
        axis equal;
        grid on;
        title('vicon');
        axis([-2 2 -2 2 -2 2]);
        disp('first plot');
    else
        set(hand(1), 'Xdata',[0, xr(1)], 'YData',[0,xr(2)],'ZData',[0,xr(3)]);
        set(hand(2), 'Xdata',[0, yr(1)], 'YData',[0, yr(2)],'ZData',[0, yr(3)]);
        set(hand(3), 'Xdata',[0, zr(1)], 'YData',[0, zr(2)],'ZData',[0, zr(3)]);
        set(hand(4), 'Xdata',[0, xr2(1)], 'YData',[0,xr2(2)],'ZData',[0,xr2(3)]);
        set(hand(5), 'Xdata',[0, yr2(1)], 'YData',[0, yr2(2)],'ZData',[0, yr2(3)]);
        set(hand(6), 'Xdata',[0, zr2(1)], 'YData',[0, zr2(2)],'ZData',[0, zr2(3)]);
    
    end
    drawnow;
end