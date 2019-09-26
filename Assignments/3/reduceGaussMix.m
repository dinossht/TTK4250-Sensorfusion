function [xmix, Pmix] = reduceGaussMix(w, x, P)
    % calculates the mean, xmix, and covariance, Pmix, of a mixture 
    %   p(y) = sum_i w(i)*N(y; x(:, i), P(:, :, i))
    %
    % w (numel(w) x #mix): weights of the mixture
    % x (dim(state) x #mix): means of the mixture
    % P (dims(state) x dim(state) x #mix): covariances of the mixture
    %
    % xmix (dim(state) x #mix): total mean
    % Pmix (dim(state) x dim(state) x #mix): total covariance
     
    w = w(:);
    M = numel(w);
    n = size(x, 1);

    %% implementation
    % allocate
    xmix = zeros(n, 1);
    Pmix = zeros(n, n);
    
    %% mean %    
    for i=1:M
        xmix = xmix+w(i)*x(:,i);
    end
    
    %% covariance
    
    % first term in covariance of gaussian mixture
    Pmix1 = zeros(n,n);
    for i=1:M
        Pmix1 = Pmix1+w(i)*P(:,:,i);
    end
    
    % second term, spread of the innovations term
    Pmix2 = zeros(n,n);
    for i=1:M
        Pmix2 = Pmix2+w(i)*x(:,i)*x(:,i)';
    end
    Pmix2 = Pmix2-xmix*xmix';

    Pmix = Pmix1+Pmix2;
    
end