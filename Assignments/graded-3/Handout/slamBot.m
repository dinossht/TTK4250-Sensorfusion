clear; close all; clc;
%rosinit;
sub = rossubscriber('slam_data');
%%
doAsso = true;

Q = diag([0.01 0.01 3*pi/180].^2);  
R = diag([0.03 10*pi/180].^2);

JCBBalphas = [0.5, 0.05]; % first is for joint compatibility, second is individual 
slam = EKFSLAM(Q, R, doAsso, JCBBalphas, [0.05, 0]);

K = 100;

%% individual

% allocate
xpred = cell(1, K);
Ppred = cell(1, K);
xhat = cell(1, K);
Phat = cell(1, K);
a = cell(1, K);

% init
xpred{1} = [0;0;0]; % we start at the correct position for reference
Ppred{1} = 0.0001*eye(3);%(3, 3); % we also say that we are 100% sure about that

%axAsso = gca;
N = K;
doAssoPlot = true; % set to true to se the associations that are done

fig = figure(4);
ax = gca;

cla(ax);
for k = 1:N
    msg = receive(sub);
    raw = split(msg.Data);
    
    dist_u = 0.01 * str2num(raw{1});
    new_ang = pi * str2num(raw{2}) / 180;
    zOdo = [dist_u; 0; new_ang];
    
    sonar = str2num(raw{3});
    z{k} = [0.01 * sonar; 0];
    
    %display(zOdo)
    
    
	display(N-k);
    
    display(sonar)
    
    if sonar == -1
        [xhat{k}, Phat{k}] = slam.predict(xpred{k}, Ppred{k}, [0;0;0]);
    else
        [xhat{k}, Phat{k}, NIS(k), a{k}] =  slam.update(xpred{k}, Ppred{k}, z{k});
    end    
    
    if k < K
        [xpred{k + 1}, Ppred{k + 1}] = slam.predict(xhat{k}, Phat{k}, zOdo);
    end
    
    
    % checks
    %if size(xhat{k},1) ~= size(Phat{k},1)
    %    error('dimensions of mean and covariance do not match')
    %end
    %{
    if doAssoPlot && k > 1 %&& any(a{k} == 0) % uncoment last part to only see new creations
        cla(axAsso); hold on;grid  on;
        zpred = reshape(slam.h(xpred{k}), 2, []);
        scatter(axAsso, z{k}(1, :), z{k}(2, :));
        scatter(axAsso, zpred(1, :), zpred(2, :));
        plot(axAsso, [z{k}(1, a{k}>0); zpred(1, a{k}(a{k}>0))], [z{k}(2, a{k}>0); zpred(2, a{k}(a{k}>0))], 'r', 'linewidth', 2)
        
        legend(axAsso, 'z', 'zbar', 'a')
        title(axAsso, sprintf('k = %d: %s', k, sprintf('%d, ',a{k})));
        %pause();
    end
    %}
    cla(ax); hold on;
    scatter(ax, 0.3, 0, 'r^')
    scatter(ax, 0.3, 0.3, 'r^')
    scatter(ax, 0.6, 0.3, 'r^')
    scatter(ax, 0.9, 0.0, 'r^')
    
    scatter(ax, xhat{k}(4:2:end), xhat{k}(5:2:end), 'b*')
    scatter(ax, xhat{k}(1), xhat{k}(2), 'rx')
    
%     if k > 1 % singular cov at k = 1
%         el = ellipse(xhat{k}(1:2),Phat{k}(1:2,1:2),5,200);
%         plot(ax,el(1,:),el(2,:),'b');
%     end
%     
%     for ii=1:((size(Phat{k}, 1)-3)/2)
%        rI = squeeze(Phat{k}(3+[1,2]+(ii-1)*2,3+[1,2]+(ii-1)*2)); 
%        el = ellipse(xhat{k}(3 + (1:2) + (ii-1)*2),rI,5,200);
%        plot(ax, el(1,:),el(2,:),'b');
%     end
%     windx = 0.5;
%     windy = 0.5;
%     xlim([xhat{k}(1)-windx, xhat{k}(1)+windx])
%     ylim([xhat{k}(2)-windy, xhat{k}(2)+windy])
    title(ax, sprintf('k = %d',k))
    grid(ax, 'on');
    

    %display(xhat{k})
end

%% run a movie
pauseTime = 2.26;
fig = figure(4);
ax = gca;
windx = 10;
windy = 2;
pause();
for k = 1:N
    cla(ax); hold on;
    %scatter(ax, landmarks(1,:), landmarks(2,:), 'r^')
    scatter(ax, 0.3, 0, 'r^')
    scatter(ax, 0.3, 0.3, 'r^')
    scatter(ax, 0.6, 0.3, 'r^')
    scatter(ax, 0.9, 0.0, 'r^')
    
    
    scatter(ax, xhat{k}(4:2:end), xhat{k}(5:2:end), 'b*')
    scatter(ax, xhat{k}(1), xhat{k}(2), 'rx')
    %plot(ax, poseGT(1, 1:k), poseGT(2,1:k), 'r-o','markerindices',10:10:k);
    plot(ax, cellfun(@(x) x(1), xhat(1:k)), cellfun(@(x) x(2), xhat(1:k)), 'b-o','markerindices',10:10:k);
    
    %xlim([xhat{k}(1)-windx, xhat{k}(1)+windx])
    %ylim([xhat{k}(2)-windy, xhat{k}(2)+windy])
    
%     
%     if k > 1 % singular cov at k = 1
%         el = ellipse(xhat{k}(1:2),Phat{k}(1:2,1:2),5,200);
%         plot(ax,el(1,:),el(2,:),'b');
%     end
%     
%     for ii=1:((size(Phat{k}, 1)-3)/2)
%        rI = squeeze(Phat{k}(3+[1,2]+(ii-1)*2,3+[1,2]+(ii-1)*2)); 
%        el = ellipse(xhat{k}(3 + (1:2) + (ii-1)*2),rI,5,200);
%        plot(ax, el(1,:),el(2,:),'b');
%     end
    
    title(ax, sprintf('k = %d',k))
    grid(ax, 'on');
    pause(pauseTime);
end

