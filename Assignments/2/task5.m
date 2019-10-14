clc; clear; close all;

%% get and plot the data
usePregen = true % choose between own generated data and pregenerated
if usePregen
    load task5data.mat
    fprintf('K = %i time steps with sampling intervall Ts = %f sec', K, Ts)
    figure(1); clf; grid on; hold on;
    % show ground truth and measurements
    plot(Xgt(1,:), Xgt(2,:));
    scatter(Z(1, :), Z(2, :));
    title('Data')
    % show turnrate
    figure(2); clf; grid on;
    plot(Xgt(5, :));
    xlabel('time step')
    ylabel('turn rate')
%%
else
    % rng(...) % random seed can be set for repeatability
    % inital state distribution
    x0 = [0, 0, 1, 1, 0]';
    P0 = diag([50, 50, 10, 10, pi/4].^2);
    % model parameters
%   % commented out to be able to run pregen without filling this in
%     qtrue = [...; ...];
%     rtrue = ...;
%     % sampling interval a lenght
%     K = ...;
%     Ts = ...;
    % get data
    [Xgt, Z] = sampleCTtraj(K, Ts, x0, P0, qtrue, rtrue);
    % show ground truth and measurements
    figure(1); clf; grid on; hold on;
    plot(Xgt(1,:), Xgt(2,:));
    scatter(Z(1, :), Z(2, :));
    title('Data')
    % show turnrate
    figure(2); clf; grid on;
    plot(Xgt(5, :));
    xlabel('time step')
    ylabel('turn rate')
end
%%
% 5 a: tune by hand and comment -- FILL IN THE DOTS

% allocate
xbar = zeros(4, K); 
xhat = zeros(4, K);
Pbar = zeros(4, 4, K);
Phat = zeros(4, 4, K);

% set parameters
q = 4;
r = 10;

% create the model and estimator object
model = discreteCVmodel(q, r);
ekf = EKF(model);

% initialize
xbar(:, 1) = [0 0 1 1]';
Pbar(:,:,1) = diag([50 50 10 10].^2);

for k = 1:K
    % estimate for this round
    [xhat(:,k),Phat(:,:,k)] = ...
        ekf.update(Z(:,k),xbar(:,k),Pbar(:,:,k));

    % Prediction for next round
    if k < K
        [xbar(:,k+1),Pbar(:,:,k+1)] = ...
            ekf.predict(xhat(:,k),Phat(:,:,k),Ts);
    end    
end

% calculate a performance metric
%RMSE = @(x, x_hat) (sqrt(mean((x' - x_hat').^2)));
posRMSE = sqrt(mean(sum((xhat(1:2,:)-Xgt(1:2,:)).^2,1))); % position RMSE
velRMSE = sqrt(mean(sum((xhat(3:4,:)-Xgt(3:4 ,:)).^2,1))); % velocity RMSE

% show results
figure(3); clf; grid on; hold on;
plot(Xgt(1,:), Xgt(2,:));
plot(xhat(1,:), xhat(2, :));
title(sprintf('q = %f, r = %f, posRMSE = %f, velRMSE= %f',q, r, posRMSE, velRMSE));
%%

% Task 5 b and c -- FILL IN THE DOTS

% parameters for the parameter grid
Nvals = 20;
qlow = 0.05;
qhigh = 10;
rlow = 0.1;
rhigh = 100;

% set the grid on logscale (not mandatory)
qs = logspace(log10(qlow), log10(qhigh), Nvals);
rs = logspace(log10(rlow), log10(rhigh), Nvals);

% allocate estimates
xbar = zeros(4, K);
Pbar = zeros(4, 4, K);
xhat = zeros(4, K);
Phat = zeros(4, 4, K);

% allocate for metrics over the grid
NIS = zeros(Nvals, Nvals, K);
NEES = zeros(Nvals, Nvals, K);


% other values of interest that can be stored
% only if you want to investigate something, like bias etc.
innovs = zeros ( Nvals , Nvals , 2 , K );
innovCorr = zeros ( Nvals , Nvals , 11 , 4) ;
NormInnov = zeros ( Nvals , Nvals , 2, K );
Ss = zeros ( Nvals , Nvals , 2 ,2, K) ;


% loop through the grid and estimate for each pair
% initialize (the same for all parameters can be used )
xbar(:,1) = [0 0 1 1]';
Pbar(:,:,1) = diag([50 50 10 10].^2) ;
% run through the grid and estimate
for qi = 1: Nvals
    display(qi/Nvals);
    for ri = 1: Nvals
        model = discreteCVmodel(qs(qi),rs(ri)) ;
        ekf = EKF(model);
        for k = 1: K
            NIS(qi,ri,k) = ekf.NIS(Z(:,k),xbar(:,k),Pbar(:,:,k));
            
            [xhat(:,k),Phat(:,:,k)] = ...
                ekf.update(Z(:,k),xbar(:,k),Pbar(:,:,k));
            
            NEES(qi,ri,k) = (Xgt(1:4,k)-xhat(:,k))'*...
                (Phat(:,:,k)\(Xgt(1:4,k)-xhat(:,k)));
            
            if k < K
                [xbar(:,k+1),Pbar(:,:,k+1)] = ...
                    ekf.predict(xhat(:,k),Phat(:,:,k),Ts);
            end
        end
        % some other quantities of interest
        [innovs(qi,ri,:,k),Ss(qi,ri,:,:,k)] = ...
            ekf.innovation(Z(:,k),xbar(:,k),Pbar(:,:,k));
        
        NormInnov(qi,ri,:,k) = ...
            chol(squeeze(Ss(qi,ri,:,:,k)))'\ ...
            squeeze(innovs(qi,ri,:,k));
        
        innovCorr(qi,ri,:,:) = ...
            xcorr(squeeze(NormInnov(qi,ri,:,:))',5) ;
    end
end

% calculate averages r x q x time
ANEES = mean(NEES,3);
ANIS = mean(NIS,3);
Ainnov = mean(innovs,4) ;

%%
% Task 5 b: ANIS plot -- FILL IN THE DOTS

% specify the  probabilities for confidence regions and calculate
alphas = [0.05 , 0.95];
CINIS = chi2inv ( alphas , 2* K )/ K;
disp(CINIS);

% plot
[qq ,rr] = meshgrid(qs, rs); % creates the needed grid for plotting
figure(4); clf; grid on;
surfc (qs , rs , ANIS' ) ; hold on ; % note transpose

caxis([0, 10])
[C, H] = contour3( qs , rs , ANIS' , CINIS ) ; % note transpose
i = 1;
while i <= size(C, 2)
    istart = i + 1;
    iend = i + C(2, i);
    % plots the countours on the surface
    plot3(C(1,istart:iend), C(2,istart:iend),ones(1,C(2, i))*C(1,i), 'r');
    i = i + C(2, i) + 1;
end
xlabel('q')
ylabel('r')
zlabel('ANIS')
zlim([0, 10])
%%
% Task 5 c: ANEES plot

% specify the  probabilities for confidence regions and calculate
alphas = [0.05 , 0.95];
CINEES = chi2inv ( alphas , 4* K )/ K;
% [...]
surfc ( qs ,rs , ANEES' ) ; hold on ; % note transpose
contour3 ( qs , rs , ANEES' , CINEES ); % note transpose

% plot
[qq ,rr] = meshgrid(qs, rs); % creates the needed grid for plotting
figure(8); clf; grid on;
surfc ( qs ,rs , ANEES' ); hold on ; % note transpose

caxis([0, 50])
[C, H] = contour3( qs , rs , ANEES', CINEES ); % note transpose

i = 1;
while i <= size(C, 2)
    istart = i + 1;
    iend = i + C(2, i);
    % plots the countours on the surface
    plot3(C(1,istart:iend), C(2,istart:iend),ones(1,C(2, i))*C(1,i), 'r');
    i = i + C(2, i) + 1;
end
xlabel('q')
ylabel('r')
zlabel('ANEES')
zlim([0, 50])
