% get and plot the data
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
q = 5 * eye(4);
r = 2 * eye(2);

% create the model and estimator object
model = discreteCVmodel(q, r);
ekf = EKF(model);

% initialize
K_gain = [eye(2),          zeros(2, 2);
     (1/Ts) * eye(2), - (1/Ts) * eye(2)];
xhat(:, 1) = K_gain * [Z(:, 1); Z(:, 2)];
Phat(:, :, 1) = [r,         1/Ts * r;
                 1/Ts * r,  2/(Ts^2) * r + Ts/3 * eye(2)];

for k = 3:(K-1)
    % estimate
    [xp, Pp] = ekf.predict(xhat(:, k), Phat(:, :, k), Ts);
    xbar(:, k) = xp;
    Pbar(:, :, k) = Pp;
    % innovate
    [vk, Sk] = ekf.innovation(Z(:, k + 1), xp, Pp);
    % update
    [xupd, Pupd] = ekf.update(Z(:, k + 1), xp, Pp);
    xhat(:, k + 1) = xupd;
    Phat(:, :, k + 1) = Pupd; 
end

% calculate a performance metric
RMSE = @(x, x_hat) (sqrt(mean((x' - x_hat').^2)));
posRMSE = RMSE(Xgt(1:2, :), xhat(1:2, :)); % position RMSE
velRMSE = RMSE(Xgt(3:4, :), xhat(3:4, :)); % velocity RMSE

% show results
figure(3); clf; grid on; hold on;
plot(Xgt(1,:), Xgt(2,:));
plot(xhat(1,:), xhat(2, :));
title(sprintf('q = %f, r = %f, posRMSE = %f, velRMSE= %f',q, r, posRMSE, velRMSE));
%%
% Task 5 b and c -- FILL IN THE DOTS

% parameters for the parameter grid
Nvals = 100;
qlow = 0.1;
qhigh = 100;
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

% initialize (the same for all parameters can be used)
xbar(:, 1) = K_gain * [Z(:, 1); Z(:, 2)];
Pbar(:, : , 1) = [r,         1/Ts * r;
                  1/Ts * r,  2/(Ts^2) * r + Ts/3 * eye(2)];


% loop through the grid and estimate for each pair
for i = 1:Nvals % q = qs
    for j = 1:Nvals % r = rs        
        % create the model and estimator object
        model = discreteCVmodel(qs(i) * eye(4), rs(j) * eye(2));
        ekf = EKF(model);
        for k = 3:(K-1)
            % estimate
            [xp, Pp] = ekf.predict(xhat(:, k), Phat(:, :, k), Ts);
            xbar(:, k) = xp;
            Pbar(:, :, k) = Pp;
            % innovate
            % [vk, Sk] = ekf.innovation(Z(:, k + 1), xp, Pp);
            % update
            [xupd, Pupd] = ekf.update(Z(:, k + 1), xp, Pp);
            xhat(:, k + 1) = xupd;
            Phat(:, :, k + 1) = Pupd; 
            NIS(i, j, k) = ekf.NIS(Z(:, k + 1), xhat(:, k + 1), Phat(:, k + 1));
        end
        
    end
end

% calculate averages
%ANEES = ...
ANIS = sum(NIS);
%%
% Task 5 b: ANIS plot -- FILL IN THE DOTS

% specify the  probabilities for confidence regions and calculate
%alphas = ...
%CINIS = ...; % the confidence bounds, Hint: inverse CDF.
%disp(CINIS);

% plot
[qq ,rr] = meshgrid(qs, rs); % creates the needed grid for plotting
figure(4); clf; grid on;
%surf(...); hold on;
caxis([0, 10])
%[C, H] = contour(...);
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
%alphas = ...
%CINEES = ...; % the confidence bounds, Hint inverse CDF
disp(CINEES);

% plot
[qq ,rr] = meshgrid(qs, rs); % creates the needed grid for plotting
figure(8); clf; grid on;
%surf(...); hold on;
caxis([0, 50])
%[C, H] = contour(...);
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
%%
% anything extra: