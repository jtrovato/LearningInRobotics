%new plot
function hand = newrotplot(R, hand)
    x= [2 0 0]';
    y= [0 2 0]';
    z= [0 0 2]';
    xr = R*x;
    yr = R*y;
    zr = R*z;
    
    if(isempty(hand))
        xhand = plot3([0, xr(1)],[0,xr(2)],[0,xr(3)], 'r', 'LineWidth', 3);
        hold on;
        yhand = plot3([0, yr(1)],[0,yr(2)],[0,yr(3)], 'g', 'LineWidth', 3);
        zhand = plot3([0, zr(1)],[0,zr(2)],[0,zr(3)], 'b', 'LineWidth', 3);
        hold off
        axis equal;
        grid on;
        axis([-2 2 -2 2 -2 2]);
    else
        xhand = hand(1); yhand = hand(2); zhand = hand(3);
        set(xhand, 'Xdata',[0, xr(1)], 'YData',[0,xr(2)],'ZData',[0,xr(3)]);
        set(yhand, 'Xdata',[0, yr(1)], 'YData',[0, yr(2)],'ZData',[0, yr(3)]);
        set(zhand, 'Xdata',[0, zr(1)], 'YData',[0, zr(2)],'ZData',[0, zr(3)]);
    end
    hand = [xhand yhand zhand];
    drawnow
end