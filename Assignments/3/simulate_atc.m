function [x,z] = simulate_atc(q, r, K, init, theCase)
Ts = 2.5;
if theCase
    rng('default');
    rng(8); % 8 good
end

models = cell(2,1);
models{1}.f = @(x, Ts) [x(1:4);0] + [Ts * x(3:4); zeros(3,1)];
models{1}.Q = @(x, Ts) [q(1) * [Ts^3/3, 0,      Ts^2/2,     0;
                        0,      Ts^3/3, 0,          Ts^2/2;
                        Ts^2/2,  0,     Ts,        0;
                        0,      Ts^2/2, 0,         Ts], zeros(4,1);
                        zeros(1, 4), 1e-20];
models{1}.sqrtQ = @(x, Ts) chol(models{1}.Q(x, Ts))';


models{2}.f = @(x, Ts) f_m2_withT(x, Ts);
models{2}.Q = @(x, Ts) [q(1) * [   Ts^3/3, 0,      Ts^2/2,	0;
                    0,      Ts^3/3, 0,    	Ts^2/2 ;
                    Ts^2/2, 0,      Ts,   	0,     ;
                    0,      Ts^2/2, 0,    	Ts,    ],   zeros(4,1);
                    zeros(1, 4),                        q(2) * Ts];
models{2}.sqrtQ = @(x, Ts) chol(models{2}.Q(x, Ts))';

 
h = @(x) x(1:2);
R = r * eye(2);
sqrtR = chol(R)';
n = 5;

xzero = init.x + chol(init.P)' * randn(n, 1);

x = zeros(size(xzero,1),K);
m = size(R,1);
z = zeros(m,floor(K));

x(:,1) = xzero;
z(:,1) = h(x(:,1)) + sqrtR*randn(m,1);

s = 1;

for k=2:K
    
    if k == 25  %First turn
        x(5, k - 1) = (pi/2)/(18 * Ts);
        s = 2;
    elseif k == 43 %Straight movement
        x(5, k - 1) = 0;
        s = 1;
    elseif k == 68 % second turn
        x(5, k - 1) = (-2*pi/3)/(15*Ts);
        s = 2;
    elseif k == 74 %Straight movement
        x(5, k - 1) = 0;
        s = 1;
    end
    
    x(:,k) = models{s}.f(x(:, k-1), Ts) + models{s}.sqrtQ(x(:, k - 1), Ts) * randn(n,1);
    z(:, k) = h(x(:,k))+ sqrtR*randn(m,1);
end
end

function xout = f_m2_withT(x,T)
    if(abs(x(5)) > 0.0001)
        xout = [x(1) + sin(T * x(5)) * x(3) / x(5) - (1 - cos(T * x(5))) * x(4) / x(5);
                x(2) + (1 - cos(T * x(5))) * x(3) / x(5) + sin(T * x(5)) * x(4) / x(5);
                cos(T * x(5)) * x(3) - sin(T * x(5)) * x(4);
                sin(T * x(5)) * x(3) + cos(T * x(5)) * x(4);...
                x(5)];
    else
        xout = [x(1) + T*x(3); x(2) + T*x(4); x(3); x(4); 0];
    end
end