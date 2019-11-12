clear; clc; close all;
load task_real.mat;
%dt = mean(diff(timeIMU));
%steps = size(zAcc,2);


%% Noise power calculation GNSS
plot_periodogram(zGNSS(1,:),1,1); hold on;
plot_periodogram(zGNSS(2,:),1,1); hold on;
sigma_x_y = 0.05;%exp(6.87/20)^2;
plot_periodogram(sigma_x_y*randn(1,3683),1,1);

sigma_z = exp(-20/20)^2;
plot_periodogram(zGNSS(3,:),1,2); hold on;
plot_periodogram(sigma_z*randn(1,3683),1,2);


%% Power squared 
gnssPower(1) = sum(gnssErr(1,:).^2/(length(gnssErr(1,:))-1));
gnssPower(2) = sum(gnssErr(2,:).^2/(length(gnssErr(2,:))-1));
gnssPower(3) = sum(gnssErr(3,:).^2/(length(gnssErr(3,:))-1));

figure(17);
clf;
title('Measurement noise: calculated vs. ideal')
subplot(311);
plot(gnssErr(1,:)); hold on; plot(sqrt(gnssPower(1))*randn(900,1)); legend('calculated','ideal')

subplot(312);
plot(gnssErr(2,:)); hold on; plot(sqrt(gnssPower(2))*randn(900,1));

subplot(313);
plot(gnssErr(3,:)); hold on; plot(sqrt(gnssPower(3))*randn(900,1));


%% Noise power calculation IMU Accelaration

N = 90000;
for k = 3500:8500
    display(k);
    
    vel1 = xtrue(4:6,k);
    vel2 = xtrue(4:6,k+1);
    R = quat2rotmat(xtrue(7:10,k));
    bias(:,k) = xtrue(11:13,k);
    
    accTrue(:,k) = inv(R)*((vel2-vel1)/(dt));
    imuAccErr(:,k) = zAcc(:,k)-accTrue(:,k)-bias(:,k);
end
%%
imuAccErr1 = [dtrend(zAcc(1,3500:8500)) dtrend(zAcc(1,35000:39000))];
imuAccErr2 = [dtrend(zAcc(2,3500:8500)) dtrend(zAcc(2,33000:39000))];
imuAccErr3 = [dtrend(zAcc(3,6000:9000)) dtrend(zAcc(3,40000:39000))];

imuAccPower(1) = sum(imuAccErr1.^2/(length(imuAccErr1)-1));
imuAccPower(2) = sum(imuAccErr2.^2/(length(imuAccErr2)-1));
imuAccPower(3) = sum(imuAccErr3.^2/(length(imuAccErr3)-1));


figure(18);
clf;
title('Measurement noise: calculated vs. ideal')
subplot(311);
plot(imuAccErr1); hold on; plot(sqrt(imuAccPower(1))*randn(length(imuAccErr1),1)); legend('calculated','ideal')

subplot(312);
plot(imuAccErr2); hold on; plot(sqrt(imuAccPower(2))*randn(length(imuAccErr2),1));

subplot(313);
plot(imuAccErr3); hold on; plot(sqrt(imuAccPower(3))*randn(length(imuAccErr3),1));


%%
imuGyroErr1 = [dtrend(zGyro(1,5000:9000)) dtrend(zGyro(1,65000:69000))];
imuGyroErr2 = [dtrend(zGyro(2,6000:9000)) dtrend(zGyro(2,67000:69000))];
imuGyroErr3 = [dtrend(zGyro(3,6000:9000)) dtrend(zGyro(3,66000:69000))];

imuGyroPower(1) = sum(imuGyroErr1.^2/(length(imuGyroErr1)-1));
imuGyroPower(2) = sum(imuGyroErr2.^2/(length(imuGyroErr2)-1));
imuGyroPower(3) = sum(imuGyroErr3.^2/(length(imuGyroErr3)-1));


figure(18);
clf;
title('Measurement noise: calculated vs. ideal')
subplot(311);
plot(imuGyroErr1); hold on; plot(sqrt(imuGyroPower(1))*randn(length(imuGyroErr1),1)); legend('calculated','ideal')

subplot(312);
plot(imuGyroErr2); hold on; plot(sqrt(imuGyroPower(2))*randn(length(imuGyroErr2),1));

subplot(313);
plot(imuGyroErr3); hold on; plot(sqrt(imuGyroPower(3))*randn(length(imuGyroErr3),1));





