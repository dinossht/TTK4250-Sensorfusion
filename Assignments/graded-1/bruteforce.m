%% 
clc; clear; close all;

%% load data
usePregen = true; % you can generate your own data if set to false
if usePregen
    load task2data.mat;
else
    K = 100;
    Ts = 2.5;
    r = 5;
    q = [0.005, 1e-6*pi]; % q(1): CV noise (effective all the time), q(2): CT noise (effective only in turns)
    init.x = [0; 0; 2; 0; 0];
    init.P = diag([25, 25, 3, 3, 0.0005].^2);
    % detection and false alarm
    PDtrue = 0.9;
    lambdatrue = 1e-4;
    [Xgt, Z, a] = simulate_atc_track(Ts, K, q, r, init, PDtrue, lambdatrue, false);
end


%% Iterating over q & r

% parameters for the parameter grid
Nvals = 10;

qlow = 0.01;
qhigh = 0.1;
rlow = 2.5;
rhigh = 7.5;

%qlow = 0.01;
%qhigh = 0.09;
%rlow = 1;
%rhigh = 9;

% set the grid on logscale (not mandatory)
qs = logspace(log10(qlow), log10(qhigh), Nvals);
rs = logspace(log10(rlow), log10(rhigh), Nvals);

% sensor 
%r = 6;
lambda = 1e-3;
PD = 0.8;
gateSize = 10^2;

% dynamic models
%qCV = 0.1;
qCT = [0.005, 0.000025];
x0 = [0; 0; 2; 0; 0];
P0 = diag([25, 25, 3, 3, 0.0005].^2);

% markov chain (you are free to parametrize this in another way)
PI11 = 0.95;
PI22 = 0.95;
p10 = 0.5;  % initial mode probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PI = [PI11, (1 - PI22); (1 - PI11), PI22]; assert(all(sum(PI, 1) == [1, 1]),'columns of PI must sum to 1')
sprobs0 = [p10; (1 - p10)]; assert(sum(sprobs0) == 1, 'initial mode probabilities must sum to 1');

% allocate
xbar = zeros(5, 2, K);
Pbar = zeros(5, 5, 2, K);
probbar = zeros(2, K);
xhat = zeros(5, 2, K);
xest = zeros(5, K);
Pest = zeros(5, 5, K);
Phat = zeros(5, 5, 2, K);
probhat = zeros(2, K);
NEES = zeros(Nvals, Nvals, K);
NEESpos = zeros(Nvals, Nvals, K);
NEESvel = zeros(Nvals, Nvals, K);

% initialize
xbar(:, :, 1) = repmat(x0, [1, 2]);
Pbar(:, : ,:, 1) = repmat(P0,[1,1,2]);
probbar(:, 1) = sprobs0;

% error init
posRMSE = zeros(Nvals,Nvals);
velRMSE = zeros(Nvals,Nvals); 
peakPosDeviation = zeros(Nvals,Nvals);
peakVelDeviation = zeros(Nvals,Nvals); 

tic;
for qi = 1:Nvals
    display(qi/Nvals);
    for ri = 1:Nvals
        % make model
        models =  cell(2,1);
        models{1} = EKF(discreteCVmodel(qs(qi), rs(ri)));
        models{2} = EKF(discreteCTmodel(qCT, rs(ri)));
        imm = IMM(models, PI);
        tracker = IMMPDAF(imm, lambda, PD, gateSize);
        

        % filter
        for k=1:K
            % Update
            [probhat(:, k), xhat(:, :, k), Phat(:, :, :, k)] = tracker.update(Z{k}, probbar(:, k), xbar(:, :, k), Pbar(:, :, :, k));

            % Total state mean and cov
            [xest(:, k), Pest(:, :, k)] =  tracker.imm.estimate(probhat(: , k), xhat(:, :, k), Phat(:, :, :, k));

            NEES(qi,ri,k) = (xest(:, k) - Xgt(:, k))' * (Pest(:, :, k) \ (xest(:, k) - Xgt(:, k)));
            NEESpos(qi,ri,k) = (xest(1:2, k) - Xgt(1:2, k))' * (Pest(1:2, 1:2, k ) \ (xest(1:2, k) - Xgt(1:2, k)));
            NEESvel(qi,ri,k) = (xest(3:4, k) - Xgt(3:4, k))' * (Pest(3:4, 3:4, k ) \ (xest(3:4, k) - Xgt(3:4, k)));
            if k < 100
                [probbar(:, k+1), xbar(:, :, k+1), Pbar(:, :, :, k+1)] = ...
                    tracker.predict(probhat(:, k), xhat(:, :, k), Phat(:, :, :, k), Ts);
            end
        end
        
        %% Calculate RMSE
        poserr = sqrt(sum((xest(1:2,:) - Xgt(1:2,:)).^2, 1));
        posRMSE(qi,ri) = sqrt(mean(poserr.^2)); % not true RMSE (which is over monte carlo simulations)
        velerr = sqrt(sum((xest(3:4, :) - Xgt(3:4, :)).^2, 1));
        velRMSE(qi,ri) = sqrt(mean(velerr.^2)); % not true RMSE (which is over monte carlo simulations)
        peakPosDeviation(qi,ri) = max(poserr);
        peakVelDeviation(qi,ri) = max(velerr);        
        
    end
end
toc;
ANEES = mean(NEES,3);
ANEESpos = mean(NEESpos,3);
ANEESvel = mean(NEESvel,3);

%% Plotting

% specify the  probabilities for confidence regions and calculate
alphas = [0.05 , 0.95];
CINEES = chi2inv ( alphas , 2* K )/ K;
disp(CINEES);

% plot
[qq ,rr] = meshgrid(qs, rs); % creates the needed grid for plotting
figure(4); clf; grid on;
surf(qs , rs , ANEES' ) ; hold on ; % note transpose


caxis([0, 10])
[C, H] = contour3( qs , rs , ANEES' , CINEES ) ; % note transpose
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
figure(4); clf; grid on;

%qs_l = linspace((qlow), (qhigh), Nvals);
%rs_l = linspace((rlow), (rhigh), Nvals);

surf(qs, rs, posRMSE' ) ; hold on ; % note transpose


xlabel('q')
ylabel('r')
zlabel('POS_RMSE')
%zlim([0, 10])