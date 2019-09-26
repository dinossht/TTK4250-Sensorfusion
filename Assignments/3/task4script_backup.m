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
r = % ... ;
qCV = % ... ;
qCT = [% ..., % ...];

% choose model to tune
s = %...;

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

% initialize filter
xbar(:, 1) = %...
Pbar(:, : ,1) = %...

% filter
for k = 1:K
    ...
end

% errors
% poserr = sqrt(sum((xhat(1:2,:) - Xgt(1:2,:)).^2, 1));
% posRMSE = 
% velerr = sqrt(sum((xhat(3:4, :) - Xgt(3:4, :)).^2, 1));
% velRMSE = 

% consistency
confidenceInterval = % ...
ANIS = mean(NIS)

% plot
figure(2); clf; hold on; grid on;
plot(xhat(1,:), xhat(2,:));
scatter(Z(1,:), Z(2, :));
% title(sprintf('posRMSE = %.3f, velRMSE = %.3f',posRMSE, velRMSE))

figure(3); clf; hold on; grid on;
plot(xhat(5,:))
% plot(Xgt(5,:))
ylabel('omega')

figure(4); clf;
%subplot(3,1,1)
plot(NISes); grid on;
ylabel('NIS')
% subplot(3,1,2);
% plot(poserr); grid on;
% subplot(3,1,3)
% plot(velerr); grid on;
%%
% tune IMM by only looking at the measurements
r = %...;
qCV = %...;
qCT = [%..., %...];
PI = [%..., %...; %..., %...];
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
xbar(:, :, 1) = repmat(%..., [1, 2]);
Pbar(:, : ,:, 1) = repmat(%...,[1,1,2]);
probbar(:, 1) = [%...; %...];

% filter
for k=1:100
    [NIS(k), NISes(:, k)] = ...
    [probhat(:, k), xhat(:, :, k), Phat(:, :, :, k)] = ...
        %...
    [xest(:, k), Pest(:, :, k)] = %...
    if k < 100
        [probbar(:, k+1), xbar(:, :, k+1), Pbar(:, :, :, k+1)] = ...
            %...
    end
end

% consistency
confidenceInterval = % ...
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
r = %...;
qCV = %...;
qCT = [%..., %...];
PI = [%...,%...; %..., %...];
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
xbar(:, :, 1) = repmat(%..., [1, 2]);
Pbar(:, : ,:, 1) = repmat(%...,[1,1,2]);
probbar(:, 1) = [%...; %...];

% filter
for k=1:100
    [NIS(k), NISes(:, k)] = %...
    
    [probhat(:, k), xhat(:, :, k), Phat(:, :, :, k)] = %...
    
    [xest(:, k), Pest(:, :, k)] = %...
    
    NEES(k) = %...
    if k < 100
        [probbar(:, k+1), xbar(:, :, k+1), Pbar(:, :, :, k+1)] = %...
    end
end

% errors
poserr = sqrt(sum((xest(1:2,:) - Xgt(1:2,:)).^2, 1));
posRMSE = %... % not true RMSE (which is over monte carlo simulations)
velerr = sqrt(sum((xest(3:4, :) - Xgt(3:4, :)).^2, 1));
velRMSE = %... % not true RMSE (which is over monte carlo simulations)
% peakPosDeviation = 
% peakVelDeviation = 

% consistency
confidenceIntervalNIS = % ... 
ANIS = mean(NIS)
confidenceIntervalNEES = % ... 
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