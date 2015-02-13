%new plot
function hand = newrotplot(R, hand)
    x= [2 0 0]';
    y= [0 2 0]';
    z= [0 0 2]';
    xr = R*x;
    yr = R*y;
    zr = R*z;
    
    if(isempty(hand))
        hand(1) = plot3([0, xr(1)],[0,xr(2)],[0,xr(3)], 'r', 'LineWidth', 2);
        hold on;
        hand(2) = plot3([0, yr(1)],[0,yr(2)],[0,yr(3)], 'g', 'LineWidth', 2);
        hand(3) = plot3([0, zr(1)],[0,zr(2)],[0,zr(3)], 'b', 'LineWidth', 2);
        hold off
        axis equal;
        grid on;
        axis([-2 2 -2 2 -2 2]);
        disp('first plot');
    else
        set(hand(1), 'Xdata',[0, xr(1)], 'YData',[0,xr(2)],'ZData',[0,xr(3)]);
        set(hand(2), 'Xdata',[0, yr(1)], 'YData',[0, yr(2)],'ZData',[0, yr(3)]);
        set(hand(3), 'Xdata',[0, zr(1)], 'YData',[0, zr(2)],'ZData',[0, zr(3)]);
    end
    drawnow;
end