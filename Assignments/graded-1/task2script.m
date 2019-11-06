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
%%
off = 0;
K = 100;


% plot measurements close to the trajectory
figure(1); clf; hold on; grid on;
Zplotdata = [];
plotMeasDist = 45;
for k = 1:K
   toPlot = false(size(Z{off+k},2),1);
   for j = 1:size(Z{off+k}, 2)
        v = Z{off+k}(:, j) - Xgt(1:2, off+k);
        toPlot(j) = v' * v <= plotMeasDist^2;
   end
   Zplotdata = [Zplotdata, Z{off+k}(:, toPlot)];
end
set(gca, 'ColorOrderIndex', 2)
scatter(Zplotdata(1,:), Zplotdata(2,:));
set(gca, 'ColorOrderIndex', 1)
plot(Xgt(1,1+off:K+off),Xgt(2,1+off:K+off), 'LineWidth',1.5);
title('True trajectory and the nearby measurements')
%%
%{
% play measurement movie. Remember that you can click the stop button.
figure(2); clf; grid on;
set(gcf,'Visible','on')
plotpause = 0.1; % sets a pause in between time steps if it goes to fast
minZ = min(Xgt(1:2, :), [], 2);
maxZ = max(Xgt(1:2, :), [], 2);
minZ = minZ - 0.05*(maxZ - minZ);
maxZ = maxZ + 0.05*(maxZ - minZ);
for k = 1:K
    scatter(Z{k}(1,:), Z{k}(2,:))
    title(sprintf('measurements at step %d', k))
    axis([minZ(1), maxZ(1), minZ(2), maxZ(2)])
    grid on;
    drawnow;
    pause(plotpause);
end
%}
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look at individual EKF-PDAs
%{
r = 5; 
qCV = 0.05;                % acceleration covariance                      
qCT = [0.5 , 0.0025];      % acceleration, turn rate covariance   

lambda = 1e-3;
PD = 0.8;
gateSize = 5^2;
% choose model to tune
s = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make models
models =  cell(2,1);
models{1} = EKF(discreteCVmodel(qCV, r));
models{2} = EKF(discreteCTmodel(qCT, r));

tracker =  cell(2,1);
tracker{1} = PDAF(models{1}, lambda, PD, gateSize);
tracker{2} = PDAF(models{2}, lambda, PD, gateSize);

% % % % allocate
xbar = zeros(5, K);
Pbar = zeros(5, 5, K);
xhat = zeros(5, K);
Phat = zeros(5, 5, K);
NEES = zeros(K, 1);
NEESpos = zeros(K, 1);
NEESvel = zeros(K, 1);

% initialize filter
xbar(:, 1) = Xgt(:,off+1);
Pbar(:, : ,1) = diag([25, 25, 3, 3, 0.0005].^2);

% filter
for k = 1:K
    [xhat(:, k) , Phat(:, :, k)] = tracker{s}.update(Z{off+k}, xbar(:, k), Pbar(:,: ,k));
    NEES(k) = ((xhat(1:3+s, k) - Xgt(1:3+s, off+k))' / squeeze(Phat(1:3+s, 1:3+s, k))) * (xhat(1:3+s, k) - Xgt(1:3+s, off+k));
    NEESpos(k) = ((xhat(1:2, k) - Xgt(1:2, off+k))' / squeeze(Phat(1:2, 1:2, k))) * (xhat(1:2, k) - Xgt(1:2, off+k));
    NEESvel(k) = ((xhat(3:4, k) - Xgt(3:4, off+k))' / squeeze(Phat(3:4, 3:4, k))) * (xhat(3:4 ,k) - Xgt(3:4 ,off+k));
    if k < K
        [xbar(:, k+1), Pbar(:, :, k+1)] = tracker{s}.predict(xhat(:, k), Phat(:, :, k), Ts);
    end
end

% errors
poserr = sqrt(sum((xhat(1:2,:) - Xgt(1:2,off+1:off+K)).^2, 1));
posRMSE = sqrt(mean(poserr.^2));
velerr = sqrt(sum((xhat(3:4, :) - Xgt(3:4, off+1:off+K)).^2, 1));
velRMSE = sqrt(mean(velerr.^2));

% consistency
chi2inv([0.025, 0.975], K*2)/K
ANEESpos = mean(NEESpos)
ANEESvel = mean(NEESvel)

chi2inv([0.025, 0.975], K*4)/K  %% 4*k = antall frihetsgrader
ANEES = mean(NEES)

% plot
figure(3); clf; hold on; grid on;
plot(xhat(1,:), xhat(2,:));
plot(Xgt(1, :), Xgt(2, :))
axis('equal')
%scatter(Z(1,:), Z(2, :));
title(sprintf('posRMSE = %.3f, velRMSE = %.3f',posRMSE, velRMSE))

figure(4); clf; hold on; grid on;
plot(xhat(5,:))
plot(Xgt(5,off+1:off+K))
legend('est','true');
ylabel('omega') 



figure(5); clf;
subplot(3,1,1);
plot(NEES); grid on; hold on;
ylabel('NEES');
ciNEES = chi2inv([0.05, 0.95], 4);
inCI = sum((NEES >= ciNEES(1)) .* (NEES <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
text(104, -5, sprintf('%.2f%% inside CI', inCI),'Rotation',90);

subplot(3,1,2);
plot(NEESpos); grid on; hold on;
ylabel('NEESpos');
ciNEES = chi2inv([0.05, 0.95], 2);
inCI = sum((NEESpos >= ciNEES(1)) .* (NEESpos <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
text(104, -5, sprintf('%.2f%% inside CI', inCI),'Rotation',90);

subplot(3,1,3);
plot(NEESvel); grid on; hold on;
ylabel('NEESvel');
ciNEES = chi2inv([0.05, 0.95], 2);
inCI = sum((NEESvel >= ciNEES(1)) .* (NEESvel <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
text(104, -5, sprintf('%.2f%% inside CI', inCI),'Rotation',90);
%}

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMM-PDA

% sensor 
r = 6;
lambda = 0.0000982474;
PD = 0.95;
gateSize = 5^2;

% dynamic models
qCV = 0.05;
qCT = [0.5, 0.001];
x0 = Xgt(:,off+1);%[0; 0; 2; 0; 0];
P0 = diag([25, 25, 3, 3, 0.0005].^2);

% markov chain (you are free to parametrize this in another way)
PI11 = 0.95;
PI22 = 0.9;
p10 = 0.9;  % initial mode probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PI = [PI11, (1 - PI22); (1 - PI11), PI22]; assert(all(sum(PI, 1) == [1, 1]),'columns of PI must sum to 1')
sprobs0 = [p10; (1 - p10)]; assert(sum(sprobs0) == 1, 'initial mode probabilities must sum to 1');


% make model
models =  cell(2,1);
models{1} = EKF(discreteCVmodel(qCV, r));
models{2} = EKF(discreteCTmodel(qCT, r));
imm = IMM(models, PI);
tracker = IMMPDAF(imm, lambda, PD, gateSize);

% allocate
xbar = zeros(5, 2, K);
Pbar = zeros(5, 5, 2, K);
probbar = zeros(2, K);
xhat = zeros(5, 2, K);
xest = zeros(5, K);
Pest = zeros(5, 5, K);
Phat = zeros(5, 5, 2, K);
probhat = zeros(2, K);
NEES = zeros(K, 1);
NEESpos = zeros(K, 1);
NEESvel = zeros(K, 1);

% initialize
xbar(:, :, 1) = repmat(x0, [1, 2]);
Pbar(:, : ,:, 1) = repmat(P0,[1,1,2]);
probbar(:, 1) = sprobs0;

% filter
for k=1:K
    % Update
    [probhat(:, k), xhat(:, :, k), Phat(:, :, :, k)] = tracker.update(Z{k}, probbar(:, k), xbar(:, :, k), Pbar(:, :, :, k));
    
    % Total state mean and cov
    [xest(:, k), Pest(:, :, k)] =  tracker.imm.estimate(probhat(: , k), xhat(:, :, k), Phat(:, :, :, k));
    
    NEES(k) = (xest(:, k) - Xgt(:, off+k))' * (Pest(:, :, k) \ (xest(:, k) - Xgt(:, off+k)));
    NEESpos(k) = (xest(1:2, k) - Xgt(1:2, off+k))' * (Pest(1:2, 1:2, k ) \ (xest(1:2, k) - Xgt(1:2, k)));
    NEESvel(k) = (xest(3:4, k) - Xgt(3:4, off+k))' * (Pest(3:4, 3:4, k ) \ (xest(3:4, k) - Xgt(3:4, k)));
    if k+off < K+off
        [probbar(:, k+1), xbar(:, :, k+1), Pbar(:, :, :, k+1)] = ...
            tracker.predict(probhat(:, k), xhat(:, :, k), Phat(:, :, :, k), Ts);
    end
end

% errors
poserr = sqrt(sum((xest(1:2,:) - Xgt(1:2,1+off:K+off)).^2, 1));
posRMSE = sqrt(mean(poserr.^2)); % not true RMSE (which is over monte carlo simulations)
velerr = sqrt(sum((xest(3:4, :) - Xgt(3:4, 1+off:K+off)).^2, 1));
velRMSE = sqrt(mean(velerr.^2)); % not true RMSE (which is over monte carlo simulations)
peakPosDeviation = max(poserr);
peakVelDeviation = max(velerr);

% consistency
CI2K = chi2inv([0.025, 0.975], K*2)/K
ANEESpos = mean(NEESpos)
ANEESvel = mean(NEESvel)

CI4K = chi2inv([0.025, 0.975], K*4)/K
ANEES = mean(NEES)

% plot
figure(6); clf; hold on; grid on;
plot(xest(1,:), xest(2,:));
plot(Xgt(1,1+off:K+off), Xgt(2, 1+off:K+off));
axis('equal')
title(sprintf('posRMSE = %.3f, velRMSE = %.3f, peakPosDev = %.3f, peakVelDev = %.3f',posRMSE, velRMSE, peakPosDeviation, peakVelDeviation))

figure(7); clf; hold on; grid on;
plot(xest(5,:))
plot(Xgt(5,1+off:K+off))

figure(8); clf;
plot(probhat');
grid on;

figure(9); clf;
subplot(2,1,1); 
plot(poserr); grid on;
ylabel('position error')
subplot(2,1,2);
plot(velerr); grid on;
ylabel('velocity error')

figure(10); clf;
subplot(3,1,1);
plot(NEES); grid on; hold on;
ylabel('NEES');
ciNEES = chi2inv([0.05, 0.95], 4);
inCI = sum((NEES >= ciNEES(1)) .* (NEES <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
text(104, -5, sprintf('%.2f%% inside CI', inCI),'Rotation',90);

subplot(3,1,2);
plot(NEESpos); grid on; hold on;
ylabel('NEESpos');
ciNEES = chi2inv([0.05, 0.95], 2);
inCI = sum((NEESpos >= ciNEES(1)) .* (NEESpos <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
text(104, -5, sprintf('%.2f%% inside CI', inCI),'Rotation',90);

subplot(3,1,3);
plot(NEESvel); grid on; hold on;
ylabel('NEESvel');
ciNEES = chi2inv([0.05, 0.95], 2);
inCI = sum((NEESvel >= ciNEES(1)) .* (NEESvel <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
text(104, -5, sprintf('%.2f%% inside CI', inCI),'Rotation',90);
%}
%{
%%
%estimation "movie"
mTL = 0.2; % maximum transparancy (between 0 and 1);
plotpause = 0;%1; % lenght to pause between time steps;
plotRange = 2:K % the range to go through
N = 50; % number of points to use for ellipse;

%k = 31; assert(all([k > 1, k <= K]), 'K must be in proper range')
figure(11); clf; grid on; hold on; axis equal;
set(gcf,'Visible','on')
co = get(gca, 'colorOrder');
minAx = min(cell2mat(Z'));
maxAx = max(cell2mat(Z'));
thetas = (0:N)*2*pi/N;
for k = plotRange
    % tracker
    subplot(1,2,1);
    cla; grid on; hold on; axis equal;
    title(sprintf('time step %d', k));
    gated = tracker.gate(Z{k}, probbar(:, k), xbar(:, :, k), Pbar(:, :, :, k));
    minG = 1e20 * ones(2,1);
    maxG = zeros(2,1);
    [skupd, xkupd, Pkupd] = tracker.conditionalUpdate(Z{k}(:, gated), probbar(:, k), xbar(:, :, k), Pbar(:, :, :, k));
    beta = tracker.associationProbabilities(Z{k}(:, gated), probbar(:, k), xbar(:, :, k), Pbar(:, :, :, k));
    for s = 1:2
        plot(squeeze(xhat(1, s, 1:(k-1))), squeeze(xhat(2, s, 1:(k-1)))','Color', co(s,:));
        axis([minAx(1), maxAx(1), minAx(2), maxAx(2)])
        for j = 1:size(xkupd, 3)
            csj = mTL * co(s, :) + (1 - mTL) * (beta(j)*skupd(s, j)*co(s, :) + (1 - beta(j)*skupd(s, j)) * ones(1, 3)); % transparancy
            pl = plot([squeeze(xhat(1, s, k-1)); xkupd(1, s, j)], [squeeze(xhat(2, s, k-1)); xkupd(2, s, j)], '--', 'Color', csj);
            axis([minAx(1), maxAx(1), minAx(2), maxAx(2)])
            %alpha(pl, beta(j)*skupd(s, j));
            drawnow;
            covEllData = chol(Pkupd(1:2,1:2, s, j))' * [cos(thetas); sin(thetas)] + xkupd(1:2, s, j);
            pl = plot(covEllData(1, :), covEllData(2, :), '--', 'Color', csj);
        end
        
        [vk, Sk] = imm.modeFilters{s}.innovation([0,0], squeeze(xbar(:, s, k)), squeeze(Pbar(:, :, s, k)));
        gateData = chol(Sk)' * [cos(thetas); sin(thetas)] * sqrt(tracker.gateSize) + squeeze(xbar(1:2, s, k));
        plot(gateData(1, :),gateData(2, :), '.--', 'Color', co(s,:))
        scatter(Z{k}(1, :), Z{k}(2, :), 'rx')
        if a(k) > 0
            scatter(Z{k}(1, a(k)), Z{k}(2, a(k)), 'gx')
        end
        minGs = min(gateData, [], 2);
        minG = minGs .* (minGs < minG) + minG .* (minG < minGs);
        maxGs = max(gateData, [], 2);
        maxG = maxGs .* (maxGs > maxG) + maxG .* (maxG > maxGs);
        
    end
    scale = 1;
    minAx = minG - scale * (maxG - minG);
    maxAx = maxG + scale * (maxG - minG);
    axis([minAx(1), maxAx(1), minAx(2), maxAx(2)])
    %legend()
    
    % mode probabilities
    subplot(1,2,2)
    cla; grid on; hold on;
    for s = 1:2
        plot(probhat(s, 1:(k - 1)), 'Color', co(s, :));
        for j = 1:size(xkupd, 3)
            csj = mTL * co(s, :) + (1 - mTL) * (beta(j)*skupd(s, j)*co(s, :) + (1 - beta(j)*skupd(s, j)) * ones(1, 3)); % transparancy
            plot([k-1, k], [probhat(s, k-1), skupd(s, j)], '--', 'color', csj)
        end
    end
    axis([1, plotRange(end), 0, 1])
    drawnow;
    pause(plotpause)
end
%}