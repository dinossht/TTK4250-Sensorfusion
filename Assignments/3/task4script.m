clc; clear; close all;
% load data
usePregen = true; % you can generate your own data if set to false
if usePregen
    load task4data.mat;
else
    K = 100;
    Ts = 2.5;
    r = 5;
    q = [0.005, 1e-6*pi]; % q(1): CV noise (effective all the time), q(2): CT noise (effective only in turns)
    init.x = [0; 0; 2; 0; 0];
    init.P = diag([25, 25, 3, 3, 0.0005].^2)
    [Xgt, Z] = simulate_atc(q, r, K, init, false);
end


figure(1); clf; hold on; grid on;
plot(Xgt(1,:), Xgt(2,:));
scatter(Z(1,:), Z(2, :));


%%
% tune single filters
r = 5;                      % pos. measurement noise covariance            
qCV = 0.05;                 % acceleration covariance                      
qCT = [0.005 , 0.000025];   % acceleration, turn rate covariance           

% choose model to tune
s = 1;                                                      

% make models
models =  cell(2,1);
models{1} = EKF(discreteCVmodel(qCV, r));
models{2} = EKF(discreteCTmodel(qCT, r));

% % % % allocate
xbar = zeros(5, 100);
Pbar = zeros(5, 5, 100);
xhat = zeros(5, 100);
Phat = zeros(5, 5, 100);
NIS = zeros(100, 1);

% initialize filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xbar(:, 1) = [0 0 1 1 1]';
Pbar(:,:,1) = diag([50 50 10 10 10].^2);

% filter
for k = 1:K
    NISes(k) = models{s}.NIS(Z(:,k),xbar(:,k),Pbar(:,:,k));
    
    % estimate
    [xhat(:,k),Phat(:,:,k)] = ...
        models{s}.update(Z(:,k),xbar(:,k),Pbar(:,:,k));

    if k < K
        [xbar(:,k+1),Pbar(:,:,k+1)] = ...
            models{s}.predict(xhat(:,k),Phat(:,:,k),Ts);
    end    
end 

% errors
poserr = sqrt(sum((xhat(1:2,:) - Xgt(1:2,:)).^2, 1));
% posRMSE = 
velerr = sqrt(sum((xhat(3:4, :) - Xgt(3:4, :)).^2, 1));
% velRMSE = 
posRMSE = sqrt(mean(sum((xhat(1:2,:)-Xgt(1:2,:)).^2,1))); % position RMSE
velRMSE = sqrt(mean(sum((xhat(3:4,:)-Xgt(3:4 ,:)).^2,1))); % velocity RMSE


%%%%%%%%%%%%
% show results
figure(3); clf; grid on; hold on;
plot(Xgt(1,:), Xgt(2,:));
plot(xhat(1,:), xhat(2, :));
title(sprintf('q = %f, r = %f, posRMSE = %f, velRMSE= %f',qCV, r, posRMSE, velRMSE));
%%%%%%%%%%%%

%%
% consistency
confidenceInterval = chi2inv([0.025,0.975],K*2)/K;
ANIS = mean(NISes);

% plot
figure(2); clf; hold on; grid on;
plot(xhat(1,:), xhat(2,:));
scatter(Z(1,:), Z(2, :));
title(sprintf('posRMSE = %.3f, velRMSE = %.3f',posRMSE, velRMSE))

figure(3); clf; hold on; grid on;
plot(xhat(5,:))
plot(Xgt(5,:))
ylabel('omega')

figure(4); clf;
subplot(3,1,1)
plot(NISes); grid on;
ylabel('NIS')
subplot(3,1,2);
plot(poserr); grid on;
subplot(3,1,3)
plot(velerr); grid on;
%%

% tune IMM by only looking at the measurements
r = 5;                      % pos. measurement noise covariance            
qCV = 0.0025;                 % acceleration covariance                      
qCT = [0.005 , 0.00005];   % acceleration, turn rate covariance 

PI11 = 0.95;
PI22 = 0.95;
PI = [PI11, (1 - PI22); (1 - PI11), PI22];
assert(all(sum(PI, 1) == [1, 1]),'columns of PI must sum to 1')

% make model
models =  cell(2,1);
models{1} = EKF(discreteCVmodel(qCV, r));
models{2} = EKF(discreteCTmodel(qCT, r));
imm = IMM(models, PI);

% allocate
xbar = zeros(5, 2, K); % dims: state, models, time
Pbar = zeros(5, 5, 2, K); % dims: state, state, models, time
probbar = zeros(2, K);
xhat = zeros(5, 2, K);
xest = zeros(5, K);
Pest = zeros(5, 5, K);
Phat = zeros(5, 5, 2, K);
probhat = zeros(2, K);
NIS = zeros(K, 1);
NISes = zeros(2, K);

% initialize
x0 = [0; 0; 2; 0; 0];
xbar(:, :, 1) = repmat(x0, [1, 2]);

P0 = diag ([25 , 25 , 3, 3, 0.0005].^2);
Pbar(:, : ,:, 1) = repmat(P0,[1,1,2]);

p10 = 0.9;
sprobs0 = [p10; (1 - p10)]
probbar(:, 1) = sprobs0;

% filter
for k=1:100
    [NIS(k), NISes(:, k)] = imm.NIS(Z(:, k), probbar(:, k),xbar(:, :, k), Pbar(:,:,:,k));

    [probhat(:, k), xhat(:, :, k), Phat(:, :, :, k)] = ...
        imm.update(Z(:, k), probbar(:, k), xbar(:, :, k), Pbar(: ,:, :, k));

    [xest(:, k), Pest(:, :, k)] = imm.estimate(probhat(:, k), xhat(: ,:, k), Phat(:, :, :, k));
    NEES(k) = (xest(1:4, k) - Xgt(1:4, k))' * (Pest(1:4, 1:4, k) \ (xest(1:4, k) - Xgt(1:4, k)));
    if k < 100
        [probbar(:, k+1), xbar(:, :, k+1), Pbar(:, :, :, k+1)] =...
            imm.predict(probhat(:, k), xhat(:, :, k), Phat(:, :, :, k), Ts);
    end
end

% consistency
confidenceInterval = chi2inv([0.025, 0.975], K * 2) / K;
ANIS = mean(NIS)

% plot
figure(5); clf; hold on; grid on;
plot(xest(1,:), xest(2,:));
scatter(Z(1,:), Z(2, :));

figure(6); clf; hold on; grid on;
plot(xest(5,:))
ylabel('omega')

figure(7); clf;
plot(probhat');
grid on;
ylabel('Pr(s)')

figure(8); clf; hold on; grid on;
plot(NIS)
plot(NISes')
ylabel('NIS')
%%
% tune IMM by looking at ground truth
r = 5;                      % pos. measurement noise covariance            
qCV = 0.0025;                 % acceleration covariance                      
qCT = [0.005 , 0.00005];   % acceleration, turn rate covariance 

PI11 = 0.95;
PI22 = 0.95;
PI = [PI11, (1 - PI22); (1 - PI11), PI22];
assert(all(sum(PI, 1) == [1, 1]),'columns of PI must sum to 1')

% make model
models =  cell(2,1);
models{1} = EKF(discreteCVmodel(qCV, r));
models{2} = EKF(discreteCTmodel(qCT, r));
imm = IMM(models, PI);

% allocate
xbar = zeros(5, 2, K);
Pbar = zeros(5, 5, 2, K);
probbar = zeros(2, K);
xhat = zeros(5, 2, K);
xest = zeros(5, K);
Pest = zeros(5, 5, K);
Phat = zeros(5, 5, 2, K);
probhat = zeros(2, K);
NIS = zeros(K, 1);
NISes = zeros(2, K);
NEES = zeros(K, 1);

% initialize
x0 = [0; 0; 2; 0; 0];
xbar(:, :, 1) = repmat(x0, [1, 2]);

P0 = diag ([25 , 25 , 3, 3, 0.0005].^2);
Pbar(:, : ,:, 1) = repmat(P0,[1,1,2]);

p10 = 0.9;
sprobs0 = [p10; (1 - p10)]
probbar(:, 1) = sprobs0;

% filter
for k=1:100
    [NIS(k), NISes(:, k)] = imm.NIS(Z(:, k), probbar(:, k),xbar(:, :, k), Pbar(:,:,:,k));

    [probhat(:, k), xhat(:, :, k), Phat(:, :, :, k)] = ...
        imm.update(Z(:, k), probbar(:, k), xbar(:, :, k), Pbar(: ,:, :, k));

    [xest(:, k), Pest(:, :, k)] = imm.estimate(probhat(:, k), xhat(: ,:, k), Phat(:, :, :, k)) ;
    NEES(k) = (xest(1:4, k) - Xgt(1:4, k))' * (Pest(1:4, 1:4, k) \ (xest(1:4, k) - Xgt(1:4, k)));    
    if k < 100
        [probbar(:, k+1), xbar(:, :, k+1), Pbar(:, :, :, k+1)] =...
            imm.predict(probhat(:, k), xhat(:, :, k), Phat(:, :, :, k), Ts);
    end
end

% errors
poserr = sqrt(sum((xest(1:2,:) - Xgt(1:2,:)).^2, 1));
posRMSE = sqrt(mean(poserr.^2));  % not true RMSE (which is over monte carlo simulations)
velerr = sqrt(sum((xest(3:4, :) - Xgt(3:4, :)).^2, 1));
velRMSE = sqrt(mean(velerr.^2)); % not true RMSE (which is over monte carlo simulations)

peakPosDeviation = max(poserr) ;
peakVelDeviation = max(velerr) ;

% consistency
confidenceIntervalNIS = chi2inv([0.025 , 0.975], K * 2) / K;
ANIS = mean(NIS)
confidenceIntervalNEES = chi2inv([0.025 , 0.975], K * 4) / K;
ANEES = mean(NEES)

% plot
figure(9); clf; hold on; grid on;
plot(xest(1,:), xest(2,:));
plot(Xgt(1,:), Xgt(2, :));
title(sprintf('posRMSE = %.3f, velRMSE = %.3f',posRMSE, velRMSE))

figure(10); clf; hold on; grid on;
plot(xest(5,:))
plot(Xgt(5,:))

figure(11); clf;
plot(probhat');
grid on;

figure(12); clf;
subplot(4,1,1); 
plot(poserr); grid on;
ylabel('position error')
subplot(4,1,2);
plot(velerr); grid on;
ylabel('velocity error')
subplot(4,1,3); hold on; grid on;
plot(NIS)
plot(NISes')
ciNIS = chi2inv([0.05, 0.95], 2);
inCI = sum((NIS >= ciNIS(1)) .* (NIS <= ciNIS(2)))/K;
plot([1,K], repmat(ciNIS',[1,2])','r--')
text(104, -2, sprintf('%.2f%% inside CI', inCI),'Rotation',90);
ylabel('NIS');
subplot(4,1,4);
plot(NEES); grid on; hold on;
ylabel('NEES');
ciNEES = chi2inv([0.05, 0.95], 4);
inCI = sum((NIS >= ciNEES(1)) .* (NIS <= ciNEES(2)))/K;
plot([1,K], repmat(ciNEES',[1,2])','r--')
text(104, -5, sprintf('%.2f%% inside CI', inCI),'Rotation',90);
