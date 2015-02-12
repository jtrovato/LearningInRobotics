clear all;
imu = load(['./ESE650 P2/imu/', 'imuRaw1.mat']);
imu_ts = imu.ts;
vals = imu.vals;

wscale =  0.016;
wbias = [373.6;375.5;370];
w = [vals(5,:);vals(6,:);vals(4,:)];
w = bsxfun(@minus, w, wbias)*wscale;

sig = zeros(1,length(imu_ts));
for j = 2:length(imu_ts)
    dt = imu_ts(j)-imu_ts(j-1);
    sig(j) = sig(j-1)+w(2,j)*dt;
end
figure()
plot(imu_ts,sig)

