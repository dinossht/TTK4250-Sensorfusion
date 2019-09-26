% trajectory generation

% scenario parameters
x0 = [pi/2, -pi/100];
Ts = 0.05;
K = round(40/Ts);

% constants
g = 9.81;
l = 1;
a = g/l;
d = 0;
S = 5;

% disturbance PDF
fpdf = makedist('uniform', 'lower', -S, 'upper', S); % disturbance PDF

% dynamic function
modulo2pi = @(x) [mod(x(1) + pi, 2 * pi) - pi; x(2)]; % loop theta to [0, 2pi] 
contPendulum = @(x) [x(2); -d*x(2) - a*sin(x(1))]; % continuous dynamics
discPendulum = @(x, v, Ts) modulo2pi(x + Ts * contPendulum(x) + Ts*[0; v]); % euler discretize with noise

% sample a trajectory
x = zeros(2, K);
x(:, 1) = x0;
for k = 1:(K-1)
    v = random(fpdf);
    x(:, k + 1) = discPendulum(x(:,k), v, Ts);
end

% vizualize
figure(1);clf;
subplot(2,1,1);
plot(x(1,:))
xlabel('Time step')
ylabel('\theta')
subplot(2,1,2)
plot(x(2,:))
xlabel('Time step')
ylabel('d/dt \theta')
%%
% measurement generation

% constants
Ld = 4;
Ll = 0;
r = 0.25;

% noise pdf
hpdf = makedist('Triangular','a',-r,'b',0,'c',r); % measurement PDF

% measurement function
h = @(x) sqrt((Ld - l * cos(x(1)))^2 + (l * sin(x(1)) - Ll)^2 ); % measurement function

Z = zeros(1, K);
for k  = 1:K
    w = random(hpdf);
    Z(k) = h(x(:,k)) + w;
end

% vizualize
figure(2); clf;
plot(Z)
xlabel('Time step')
ylabel('z')
%%
% Task: Estimate the pendulum state using a particle filter 
% -- FILL IN THE DOTS

% number of particles to use (tuning)
N = 1000;

% initialize particles, pretend you do not know where the pendulum starts
px = [2* pi *( rand(1,N) - 0.5) ; randn(1,N) * pi /4];

% initial weights
w = ones(N,1)/N;

% allocate a variable for resampling particles
pxn = zeros(size(px));

% PF transition PDF: SIR proposal, or something you would like to test
% (tuning)
PFfpdf = makedist('uniform', 'lower', -S, 'upper', S);

% initialize a figure for particle animation.
figure(4); clf; grid on; hold on;
set(gcf,'Visible','on')
plotpause = 0;

% estimate
for k = 1:K
    % weight update
    for n = 1:N
        w(n) =  pdf ( hpdf , Z( k) - h( px (: , n) )); % write help pdf
    end
    w = w/sum(w); % normalize
    
    % resample
    cumw = cumsum (w) ;
    i = 1;
    for n = 1:N
        % find a particle i to pick
        % algorithm in the book, but there are other options as well
        u = (n - 1) /N; 
        while u > cumw (i)
            i = i + 1;
        end
        pxn(:, n) = px(:, i);
    end
    
    % trajecory sample prediction
    for n = 1:N
        px(:, n) = discPendulum ( pxn (: , n) , random ( PFfpdf ) , Ts );
    end
    
    %plot
    clf; grid on; hold on;
    scatter(l * sin(pxn(1,:)), -l * cos(pxn(1,:)),'b.');
    sh = scatter(l * sin(x(1, k)), -l * cos(x(1,k)), 'rx');
    axis([-l,l,-l,l]*1.5)
    xlabel('x')
    ylabel('y')
    title('theta mapped to x-y')
    legend('particl','\theta true')
    drawnow;
    pause(plotpause);
end
%}