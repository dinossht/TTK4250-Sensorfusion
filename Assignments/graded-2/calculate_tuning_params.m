clear; clc; close all;
load task_simulation.mat;
dt = mean(diff(timeIMU));
steps = size(zAcc,2);

%% Measurement noise
% GNSS Position  measurement
p_std =  [0.300    0.300    0.508]'; % Measurement noise
RGNSS = diag(p_std.^2);

% accelerometer
qA = 10*1.4101e-04; % accelerometer measurement noise covariance
qAb = 1000*1.4101e-05; % accelerometer bias driving noise covariance
pAcc = 1e-9;  % accelerometer bias reciprocal time constant

qG = 1000*2.0552e-07; % gyro measurement noise covariance
qGb = 5000*2.0552e-08;  % gyro bias driving noise covariance
pGyro = 1e-9; % gyrp bias reciprocal time constant

%% Estimator
eskf = ESKF(qA, qG, qAb, qGb, pAcc, pGyro);

%% Noise power calculation GNSS

N = 90000;
GNSSk = 1;
for k = 1:N
    display(k);
    if  timeIMU(k) >= timeGNSS(GNSSk)
        [gnssErr(:,GNSSk),~] = eskf.innovationGNSS(xtrue(:,k), ones(15,15), zGNSS(:,GNSSk), RGNSS, leverarm);
        GNSSk = GNSSk  + 1;        
    end
end

%% Power squared 
gnssPower(1) = sum(gnssErr(1,:).^2/(length(gnssErr(1,:))-1));
gnssPower(2) = sum(gnssErr(2,:).^2/(length(gnssErr(2,:))-1));
gnssPower(3) = sum(gnssErr(3,:).^2/(length(gnssErr(3,:))-1));

figure(17);
clf;
title('Measurement noise: calculated vs. ideal')
subplot(311);
plot(gnssErr(1,:)); hold on; plot(sqrt(gnssPower(1))*randn(900,1)); legend('calculated','ideal'); ylim([-1 1]);

subplot(312);
plot(gnssErr(2,:)); hold on; plot(sqrt(gnssPower(2))*randn(900,1));ylim([-1 1]);

subplot(313);
plot(gnssErr(3,:)); hold on; plot(sqrt(gnssPower(3))*randn(900,1));ylim([-1 1]);


%% Noise power calculation IMU Accelaration
N = 90000;
for k = 1:N-2
    display(k);
    
    vel1 = xtrue(4:6,k);
    vel2 = xtrue(4:6,k+2);
    R = quat2rotmat(xtrue(7:10,k));
    bias_g(:,k) = xtrue(11:13,k);
    
    bias_a(:,k) = xtrue(14:16,k);
    
    accTrue(:,k) = inv(R)*((vel2-vel1)/(2*dt)-[0; 0; 9.82]);    
    imuAccErr(:,k) = zAcc(:,k)-accTrue(:,k)-bias_g(:,k);
end

%%  Noise power calculation IMU Accelaration
imuAccErr = accTrue - zAcc(:,2:end-1);

imuAccErr(1,:) = imuAccErr(1,:) - movmean(imuAccErr(1,:),5);
imuAccErr(2,:) = imuAccErr(2,:) - movmean(imuAccErr(2,:),5);
imuAccErr(3,:) = imuAccErr(3,:) - movmean(imuAccErr(3,:),5);

imuAccErr(1,:) = min(0.05, max(-0.05, imuAccErr(1,:)));
imuAccErr(2,:) = min(0.05, max(-0.05, imuAccErr(2,:)));
imuAccErr(3,:) = min(0.05, max(-0.05, imuAccErr(3,:)));

imuAccPower(1) = sum(imuAccErr(1,:).^2/(length(imuAccErr(1,:))-1));
imuAccPower(2) = sum(imuAccErr(2,:).^2/(length(imuAccErr(2,:))-1));
imuAccPower(3) = sum(imuAccErr(3,:).^2/(length(imuAccErr(3,:))-1));

display(imuAccPower);
m = size(imuAccErr,2);
figure(18);
clf;
title('Measurement noise: calculated vs. ideal')
subplot(311);
plot(imuAccErr(1,:)); hold on; plot(sqrt(imuAccPower(1))*randn(m,1)); legend('calculated','ideal'); ylim([-0.06 0.06]);

subplot(312);
plot(imuAccErr(2,:)); hold on; plot(sqrt(imuAccPower(2))*randn(m,1));legend('calculated','ideal'); ylim([-0.06 0.06]);

subplot(313);
plot(imuAccErr(3,:)); hold on; plot(sqrt(imuAccPower(3))*randn(m,1));legend('calculated','ideal'); ylim([-0.06 0.06]);


%%  Noise power calculation IMU Gyros
imuGyroErr(1,:) = zGyro(1,:)-movmean(zGyro(1,:),5);
imuGyroErr(2,:) = zGyro(2,:)-movmean(zGyro(2,:),5);
imuGyroErr(3,:) = zGyro(3,:)-movmean(zGyro(3,:),5);

imuGyroErr(1,:) = min(2e-3, max(-2e-3, imuGyroErr(1,:)));
imuGyroErr(2,:) = min(2e-3, max(-2e-3, imuGyroErr(2,:)));
imuGyroErr(3,:) = min(2e-3, max(-2e-3, imuGyroErr(3,:)));

imuGyroPower(1) = sum(imuGyroErr(1,:).^2/(length(imuGyroErr(1,:))-1));
imuGyroPower(2) = sum(imuGyroErr(2,:).^2/(length(imuGyroErr(2,:))-1));
imuGyroPower(3) = sum(imuGyroErr(3,:).^2/(length(imuGyroErr(3,:))-1));

display(imuGyroPower);
m = size(imuGyroErr,2);
figure(19);
clf;
title('Measurement noise: calculated vs. ideal')
subplot(311);
plot(imuGyroErr(1,:)); hold on; plot(sqrt(imuGyroPower(1))*randn(m,1)); legend('calculated','ideal'); ylim([-2e-3 2e-3]);

subplot(312);
plot(imuGyroErr(2,:)); hold on; plot(sqrt(imuGyroPower(2))*randn(m,1));legend('calculated','ideal'); ylim([-2e-3 2e-3]);

subplot(313);
plot(imuGyroErr(3,:)); hold on; plot(sqrt(imuGyroPower(3))*randn(m,1));legend('calculated','ideal'); ylim([-2e-3 2e-3]);
