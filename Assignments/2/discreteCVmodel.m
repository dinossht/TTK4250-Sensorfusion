function model = discreteCVmodel(q, r)
    % returns a structure that implements a discrete time CV model with
    % continuous time accelleration covariance q and positional
    % measurement with noise with covariance r, both in two dimensions.
    
    % discrete prediction function
    model.f = @(x, Ts) [1 0 Ts 0;
                        0 1 0 Ts;
                        0 0 1 0;
                        0 0 0 1]*x;
                    
    % jacobian of prediction function
    model.F = @(x, Ts) [1 0 Ts 0;
                        0 1 0 Ts;
                        0 0 1 0;
                        0 0 0 1]; 
                    
    % additive discrete noise covariance
    model.Q = @(x, Ts) q * [Ts^3 / 3,   0,          Ts^2 / 2,   0;
                            0,          Ts^3 / 3,   0,          Ts^2 / 2;
                            Ts^2 / 2,   0,          Ts,         0;
                            0,          Ts^2 / 2,   0,          Ts];

	% measurement function                        
    model.h = @(x) [1, 0, 0, 0;
                    0, 1, 0, 0] * x; 

    % measurement function jacobian
    model.H = @(x) [1, 0, 0, 0;
                    0, 1, 0, 0];  
    % additive measurement noise covariance
    model.R = r; 
end