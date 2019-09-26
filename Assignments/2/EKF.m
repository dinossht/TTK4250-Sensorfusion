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
            xp = obj.f(x,Ts);  % line2,alg2          
            Fk = obj.F(x,Ts);  % line3,alg2
            Pp = Fk*P*Fk'+obj.Q(x,Ts);  % line4,alg2
        end

        function [vk, Sk] = innovation(obj, z, x, P)
            % returns the innovation and innovation covariance
            vk = z-obj.h(x);
            Hk = obj.H(x);
            Sk = Hk*P*Hk'+obj.R;
        end

        function [xupd, Pupd] = update(obj, z, x, P)
            % returns the mean and covariance after conditioning on the
            % measurement
            [vk, Sk] = obj.innovation(z, x, P);
            Hk = obj.H(x);
            Wk = P*Hk'/Sk;
            xupd = x+Wk*vk;
            I = eye(size(P));
            Pupd = (I-Wk*Hk)*P; 
        end

        function NIS = NIS(obj, z, x, P)
            % returns the normalized innovation squared
            [vk, Sk] = obj.innovation(z, x, P);
            NIS = vk'*inv(Sk)*vk;
        end

        function ll = loglikelihood(obj, z, x, P)
            % returns the logarithm of the marginal mesurement distribution
            % This is equivalent to implement the logarithm of the distribution of the innovation,
            % log(N (vk; 0, HP HT + R));
            [vk,Sk] = obj.innovation(z,x,P);
            NIS = obj.NIS(z,x,P);
            ll = -0.5*(NIS+log(det(2*pi*Sk)));
        end

    end
end