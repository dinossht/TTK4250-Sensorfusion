classdef EKFSLAM
    properties
        Q % 
        R % støy i landmark måling
        doAsso
        alpha
        sensOffset
    end
    methods
        function obj = EKFSLAM(Q, R, doAsso, alphas, sensorOffset)
            obj.Q = Q;
            obj.R = R;
            
            if nargin < 3
                doAsso = false;
            end
            obj.doAsso = doAsso;
            
            if nargin < 4
                alphas = [0.001, 0.0001];
            end
            obj.alpha = alphas;
            
            if nargin < 5
                sensorOffset = zeros(2,1);
            end
            obj.sensOffset = sensorOffset(:);
        end
        
        function xpred = f(~, x, u)
            % takes a pose and odometry and predicts it to the next time step
            % (eq: 11.7)
            xpred = x + [rotmat2d(x(3))*u(1:2); u(3)];
            xpred(3) = wrapTo2Pi(xpred(3));
        end
        
        function Fx = Fx(~, x, u)
            % Jacobian of the above function with respect to the pose
            % (eq: 11.13)
            Fx = [
                1   0  -u(1)*sin(x(3))-u(2)*cos(x(3));
                0   1   u(1)*cos(x(3))-u(2)*sin(x(3));
                0   0   1                              
                ];
            
            % check that jacobian is correct, remove for speed
            %if norm(Fx - jacobianFD(@(X) obj.f(X, u), x, 1e-5), 'fro') > 1e-3
            %    error('some error in pred Jac')
            %end
        end
        
        function Fu = Fu(~, x, u)
             % that is the Jacobian of the above function with respect to the odometry
             % (eq: 11.14)
             
              Fu = [
                  cos(x(3)) -sin(x(3)) 0;
                  sin(x(3))  cos(x(3)) 0;
                  0          0         1
                  ];
              % check that the jacobian is correct, remove for speed
            %if norm(F - jacobianFD(@(U) obj.f(x, U), x, 1e-5), 'fro') > 1e-3
            %    error('some error in pred Jac')
            %end            
        end
        
        function [etapred, P] =  predict(obj, eta, P, zOdo)
            x = eta(1:3); % pose
            m = eta(4:end); % map
            
            % takes the state (pose and map) , its covariance and 
            % odometry to predict the state and covariance
            % (eq: 11.19)
            xpred   = obj.f(x, zOdo);
            Fx      = obj.Fx(x, zOdo);
            Fu      = obj.Fu(x, zOdo);
            
            % in place for performance
            P(1:3, 1:3) = Fx*P(1:3, 1:3)*Fx' + Fu*obj.Q*Fu';
            P(1:3, 4:end) = Fx*P(1:3, 4:end);
            P(4:end, 1:3) = P(1:3, 4:end)'; 
            
            % concatenate pose and landmarks again
            etapred = [xpred; m];
            
            % check that the covariance makes sense
            %if any(eig(Ppred) <= 0) % costly, remove when tested
            %    warn('EKFpredict got cov not PSD')
            %end
        end
        
        function zpred = h(obj, eta)
            % that predicts the measurements of all
            % the landmarks in the state, corresponding to equations (11.10) - (11.11)
            x = eta(1:3); % pose 
            m = reshape(eta(4:end), 2 ,[]); % map (2 x m now)
            
            Rot = rotmat2d(-x(3)); % rot from world to body
            
            % cartesian measurement in world
            z_c = m - x(1:2) - Rot' * obj.sensOffset;
            
            % in body
            z_b = Rot*z_c;
            
            % polar (use maybe cart2pol)
            for i=1:size(m,2)
                zpred(:,i) = [norm(z_c(:,i),2); atan2(z_b(2,i),z_b(1,i))];
            end
            %zpred(1:2,:) = cart2pol(z_b(1,:),z_b(2,:));
            
            % make column again
            zpred = zpred(:); 
        end
        
        function H = H(obj, eta)
            % the jacobian of the above measurement function, EKFSLAM.H(obj, eta), 
            % given by equa-tions (11.16) - (11.18)
            x = eta(1:3); % pose
            m = reshape(eta(4:end), 2 ,[]); % map
            
            numM = size(m, 2);  % number of landmarks
            
            Rot = rotmat2d(x(3));
            
            m_minus_rho = m - x(1:2);
            z_c = m_minus_rho - Rot * obj.sensOffset;
            
            zpred = reshape(obj.h(eta), 2, []);
            zr = zpred(1,:);
            
            Rpihalf = [0, -1; 1, 0]; 
            
            % allocate
            Hx  = zeros(2 * numM, 3); % pose columns
            Hm = zeros(2 * numM, 2 * numM); % map columns (the rest)
            
            I2 = eye(2);
            
            for i = 1:numM
                inds = 2*(i - 1) + [1; 2];
                
                jac_z_b = [-I2 -Rpihalf*m_minus_rho(:,i)];
                
                Hx(inds(1), :) = z_c(:,i)'/(zr(i))*jac_z_b;  % jac z_r 
                Hx(inds(2), :) = z_c(:,i)'*Rpihalf'/(zr(i)^2)*jac_z_b;  % jac z_phi
            
                Hm(inds, inds) = -Hx(inds,1:2);  %... should be negative of the two first colums of Hx
            end
        
            % concatenate the H matrix
            H = [Hx, Hm];
            
            % check that it is done correctly, remove for speed
            %if norm(H - jacobianFD(@(X) obj.h(X), eta, 1e-5), 'fro') > 1e-3
            %    error('some error in meas Jac')
            %end
        end
        
        function [etaadded, Padded] = addLandmarks(obj, eta, P, z)
            % Implement the function that inverts the masurement function and creates new landmarks and their
            % covariances i
            n = size(P, 1);
            numLmk = numel(z)/2;
            
            % allocate
            lmnew = zeros(size(z));
            Gx = zeros(numLmk * 2, 3);
            Rall = zeros(numLmk * 2, numLmk * 2);
            
            I2 = eye(2);
            R_plus_pihalf = rotmat2d(eta(3)+ pi/2);

            for j = 1:numLmk
                % find indeces and the relevant measurement
                inds = 2 * (j - 1) + [1, 2];
                zj = z(inds);
                
                rot = rotmat2d(zj(2) + eta(3));

                lmnew(inds) = zj(1)*rot(:,1) + eta(1:2) + rotmat2d(eta(3))*obj.sensOffset;
                %B2W = rotmat2d(eta(3));
                %lmnew(inds) = eta(1:2) + B2W*([zj(1)*cos(zj(2)); zj(1)*sin(zj(2))] + obj.sensOffset);
                
                Gx(inds, :) =  [I2 zj(1)*rot(:,2) + R_plus_pihalf*obj.sensOffset]; % jac h^-1 wrt. x
                Gz =  rot*diag([1,zj(1)]);% jac h^-1 wrt. z
                
                Rall(inds, inds) = Gz*obj.R*Gz'; % the linearized measurement noise

            end
            
            % augment state
            etaadded = [eta;lmnew];
            
            % add covariances
            Padded = blkdiag(P, Gx*P(1:3,1:3)*Gx' + Rall);
            Padded((n+1):end, 1:n) = Gx*P(1:3,:);
            Padded(1:n, (n+1):end) = Padded((n+1):end, 1:n)';

            % sanity check, remove for speed
            %if any(eig(Padded) <= 0) % costly, remove when tested
            %    warning('EKFupdate got cov not PSD after adding a landmark');
            %end
        end
        
        function [z, zpred, H, S, a] = associate(obj, z, zpred, H, S)
            if obj.doAsso
                % associate
                a = JCBB(z, zpred, S, obj.alpha(1), obj.alpha(2));

                % extract associated measurements
                zinds = false(size(z));
                zinds(1:2:end) = a > 0;
                zinds(2:2:end) = zinds(1:2:end);
                z = z(zinds);

                % extract and rearange predicted measurements and cov
                zbarinds = reshape([2*a(a>0) - 1, 2*a(a>0)]', [], 1);
                zpred = zpred(zbarinds);
                S = S(zbarinds, zbarinds);
                H = H(zbarinds, :);
            % else  % if no association is to be done, assume that all measurements are there and in order
            end
        end
        
        function [etaupd, Pupd, NIS, a] = update(obj, eta, P, z)    
            % takes the prior mean and
            % covariance and a set of measurements to update the map and pose estimates as well as calculating
            % NIS and creating new landmarks. The data association is done for you so you only need to make
            % make the EKF update.
            numLmk = (numel(eta) - 3)/2; % number of landmarks
            if numLmk > 0
                % prediction and innovation covariance
                zpred = obj.h(eta);
                H = obj.H(eta);
                
                
                S = H*P*H' + kron(eye(numLmk),obj.R);
                z = z(:); % vectorize
                
                % perform data association if it is asked for
                [za, zpred, H, S, a] = obj.associate(z, zpred, H, S);

                % create the associated innovation
                v = za(:) - zpred;
                v(2:2:end) = wrapToPi(v(2:2:end)); % angles are in [-pi, pi]

                % Kalman update
                W = P*H'/S;
                etaupd = eta + W*v;
                NIS = v'*inv(S)*v;
                Pupd = (eye(size(P,1)) - W*H)*P;
                
                % sanity check, remove for speed
                %if any(eig(Pupd) <= 0) % costly, remove when tested
                %    warn('EKFupdate got cov not PSD');
                %end
            else % all measurements are new landmarks
                a = zeros(size(z, 2), 1);
                z = z(:);
                NIS = 0;
                etaupd = eta;
                Pupd = P;
            end
            
            % create new landmarks if any is available
            if obj.doAsso
                isNewLmk = (a == 0);
                if any(isNewLmk)
                    % extract unassociated measurements
                    zNewInds = false(size(z));
                    zNewInds(1:2:end) = isNewLmk;
                    zNewInds(2:2:end) = isNewLmk;
                    znew = z(zNewInds);
                    
                    % create new landmarks
                    [etaupd, Pupd] = obj.addLandmarks(eta, P, znew);
                end
            end  
        end
    end
end
