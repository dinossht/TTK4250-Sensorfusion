% initialize the values
x_bar = zeros(2,1)
P = 25 * eye(2)

H_r = eye(2)
H_c = eye(2)

R_r = [79, 36; 36, 36]
R_c = [28, 4; 4, 22]

z_c = [2; 14]
z_r = [-4; 6]

% set up for plotting ellipses
Npts = 100; % number of points around the circle
circle = [cos( 0:(2*pi/Npts):(2*pi) ); sin( 0:(2*pi/Npts):(2*pi) )];

% initial 1 sigma ellipses
figure(1); clf; hold on; grid on;
data = x_bar + chol(P)' * circle;
plot(data(1,:), data(2, :), 'DisplayName','prior')

% measurements 
scatter(z_c(1), z_c(2), 'DisplayName', 'z_c')
scatter(z_r(1), z_r(2), 'DisplayName', 'z_r')

legend()

%%
% You can make some functions to ease the conditionioning
% FILL IN THE DOTS ...
condition_mean = @(x_bar, z, P, H, R) x_bar+P*H'*inv(H*P*H'+R)*(z-H*x_bar);
condition_cov = @(P, H, R) P-P*H'*inv(H*P*H'+R)*H*P';
%%
% task 5 (f)
% FILL IN FOR THE DOTS ...

% condition on camera
x_bar_c = condition_mean(x_bar,z_c,P,H_c,R_c);
P_c = condition_cov(P,H_c,R_c);

% condition on radar
x_bar_r = condition_mean(x_bar,z_r,P,H_r,R_r);
P_r = condition_cov(P,H_r,R_r);

% Plot 1 sigma ellipses
figure(2); clf; hold on; grid on;
data = x_bar + chol(P)' * circle;
plot( data(1,:), data(2, :) , 'DisplayName','prior')

data = x_bar_c + chol(P_c)' * circle;
plot( data(1,:), data(2, :) , 'DisplayName', 'c')

data = x_bar_r + chol(P_r)' * circle; 
plot( data(1,:), data(2, :) , 'DisplayName', 'r')

% measurements
scatter(z_c(1), z_c(2), 'DisplayName', 'z_c')
scatter(z_r(1), z_r(2), 'DisplayName', 'z_r')
legend()

%%
% task 5 (g)

% condition the already camera conditioned on the radar
x_bar_cr = condition_mean(x_bar_c,z_r,P_c,H_r,R_r);
P_cr = condition_cov(P_c,H_r,R_r);

% condition the already radar conditioned on the camera
x_bar_rc = condition_mean(x_bar_r,z_c,P_r,H_c,R_c);
P_rc = condition_cov(P_r,H_c,R_c);

% Plot 1 sigma ellipses
figure(3); clf; hold on; grid on;

... VISUALIZE ALL ELLIPSES HERE ...

data = x_bar_c + chol(P_c)' * circle;
plot( data(1,:), data(2, :) , 'b', 'DisplayName', 'c')

data = x_bar_r + chol(P_r)' * circle; 
plot( data(1,:), data(2, :) , 'm', 'DisplayName', 'r')

data = x_bar_rc + chol(P_rc)' * circle; 
plot( data(1,:), data(2, :) ,'r', 'DisplayName', 'rc')

data = x_bar_cr + chol(P_cr)' * circle; 
plot( data(1,:), data(2, :) ,'k', 'DisplayName', 'cr')

% meausrements
scatter(z_c(1), z_c(2), 'b', 'DisplayName', 'z_c')
scatter(z_r(1), z_r(2), 'm', 'DisplayName', 'z_r')
scatter(x_bar_rc(1), x_bar_rc(2), 'r', 'DisplayName', 'x_r_c')
scatter(x_bar_cr(1), x_bar_cr(2), 'k', 'DisplayName', 'x_c_r')

legend()