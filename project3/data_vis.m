% visualize data
close all;

label = 1;
trial = 1;

d = data(label).trials;
ts = d(:,1);
s = 1:length(ts);


% figure();
% plot(s,d(:,2),s,d(:,3),s,d(:,4),s,d(:,5),s,d(:,6),s,d(:,7));
% legend('a_x','a_y','a_z','\omega_x','\omega_y','\omega_z');
figure();
plot(s,d(:,2),s,d(:,3),s,d(:,4));
legend('a_x','a_y','a_z');
% figure();
% plot(s,d(:,5),s,d(:,6),s,d(:,7));
% legend('\omega_x','\omega_y','\omega_z');
%[x,~,button] = ginput;

% 
% figure();
% plot(ts,d(:,2),ts,d(:,3),ts,d(:,4),ts,d(:,5),ts,d(:,6),ts,d(:,7));
% legend('a_x','a_y','a_z','\omega_x','\omega_y','\omega_z');
% figure();
% plot(ts,d(:,2),ts,d(:,3),ts,d(:,4));
% legend('a_x','a_y','a_z');
% figure();
% plot(ts,d(:,5),ts,d(:,6),ts,d(:,7));
% legend('\omega_x','\omega_y','\omega_z');

%%
K = 7;
[labels, centroids] = kmeans(d(:,2:4), K); %kmeans on the acceleration

dis = d(2:4);
for i=1:K
    el = sum(labels==i);
    dis(labels==i, :) = repmat(centroids(i,:), el, 1);
end
figure();
plot(s,dis(:,1),s,dis(:,2),s,dis(:,3));
legend('a_x','a_y','a_z');


% figure();
% grid on;
% plot3(d(:,2),d(:,3),d(:,4),'b*');
% k=4;
% figure();
% grid on;
% for i=1:k
%     inds = (labels == k);
%     plot3(d(inds,2), d(inds,3), d(inds,4), '*');
%     hold on;
% end