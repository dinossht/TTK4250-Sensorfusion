classdef EKF
    % FILL IN THE DOTS
    properties
        model
        
        f % discrete prediction function
        F % jacobian of prediction function
        Q % additive discrete noise covariance
        
        h % measurement function
        H % measurement function jacobian
        R % additive measurement noise covariance
    end
    methods 
        function obj = EKF(model)
            obj = obj.setModel(model);
        end
        
        function obj = setModel(obj, model)
           % sets the internal functions from model
           obj.model = model;
           
           obj.f = model.f;
           obj.F = model.F;
           obj.Q = model.Q;
           
           obj.h = model.h;
           obj.H = model.H;
           obj.R = model.R;
        end
        
        function [xp, Pp] = predict(obj, x, P, Ts)
            % returns the predicted mean and covariance for a time step Ts
            Fk = obj.F(x, Ts);
            
            xp = obj.f(x, Ts);
            Pp = Fk * P * (Fk') + obj.Q(x, Ts);
        end

        function [vk, Sk] = innovation(obj, z, x, P)
            % returns the innovation and innovation covariance
            Hk = obj.H(x);
            
            % Assuming z = z_{k}
            % Assuming x = x_{k|k-1}
            vk = z - obj.h(x); 
            % Assuming P = P_{k|k-1}
            % Assuming R is implemented as such
            Sk = Hk * P * (Hk') + obj.R;
        end

        function [xupd, Pupd] = update(obj, z, x, P)
            % returns the mean and covariance after conditioning on the
            % measurement
            
            % Same assumptions as above
            [vk, Sk] = obj.innovation(z, x, P);
            Hk = obj.H(x);
            I = eye(size(P));
            
            Wk = P * (Hk') / Sk;

            xupd = x + Wk * vk;
            Pupd = (I - Wk * Hk) * P;
        end

        function NIS = NIS(obj, z, x, P)
            % returns the normalized innovation squared
            [vk, Sk] = obj.innovation(z, x, P);
            
            NIS = vk / Sk * vk;
        end

        function ll = loglikelihood(obj, z, x, P)
            % returns the logarithm of the marginal mesurement distribution
            [vk, Sk] = obj.innovation(z, x, P);
            NIS = obj.NIS(z, x, P);

            ll = -0.5 * (NIS + log(det(2 * pi * Sk)));
        end

    end
end