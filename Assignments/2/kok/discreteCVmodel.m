function model = discreteCVmodel(q, r)
    % returns a structure that implements a discrete time CV model with
    % continuous time accelleration covariance q and positional
    % measurement with noise with covariance r, both in two dimensions. 
    
    model.f = @(x, Ts) [1, 0, Ts, 0;
                        0, 1, 0, Ts;
                        0, 0, 1, 0;
                        0, 0, 0, 1] * x;
    model.F = @(x, Ts) [1, 0, Ts, 0;
                        0, 1, 0, Ts;
                        0, 0, 1, 0;
                        0, 0, 0, 1];
    % in the CV model, assuming q = sigma^2 eye(2)
    model.Q = @(x, Ts) q * [Ts^3 / 3,   0,          Ts^2 / 2,   0;
                            0,          Ts^3 / 3,   0,          Ts^2 / 2;
                            Ts^2 / 2,   0,          Ts,         0;
                            0,          Ts^2 / 2,   0,          Ts];
                    
    model.h = @(x) [1, 0, 0, 0;
                    0, 1, 0, 0] * x;
    model.H = @(x) [1, 0, 0, 0;
                    0, 1, 0, 0];
    model.R = r;
end