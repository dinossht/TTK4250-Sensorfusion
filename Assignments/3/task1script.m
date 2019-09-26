% fill in for mu, sigma and w. This script is only made to deal with 3
% compent mixtures.
mus = [0,2,4.5]; mus = reshape(mus, [1,3]);
sigmas = [1,1,1]; sigmas = reshape(sigmas, [1,1,3]);% note std and not variance as used in mixture reduction!
w = [1/3,1/3,1/3]; w = w(:)/sum(w(:));

% get overall mean and sigma for finding interval of plotting
[totMean, totSigma2] = reduceGaussMix(w, mus, sigmas.^2);
plotNsigmas = 3;
x = totMean + plotNsigmas * sqrt(totSigma2) * ((0:1000) - 500) / 500;

% plot individual distributions
figure(1);clf; grid on; hold on;
pdftot = zeros(size(x));
for i = 1:3
    N{i} = makedist('Normal', 'mu', mus(i), 'sigma', sigmas(i)); % create component pdf
    Nx{i} = pdf(N{i}, x); % evaluate component pdf
    pdftot = pdftot + w(i)*Nx{i}; % add up total pdf
    plot(x, w(i)*Nx{i}, 'DisplayName', sprintf('component %d', i))
end
legend()

% plot the total distribution before and after merging.
figure(2); clf; grid on; hold on;
plot(x, pdftot, 'DisplayName', 'original')
k = 1;
mucomb = zeros(3,1);
sigma2comb = zeros(3, 1);
for i = 1:2 % index of first to merge
    for j = (i+1):3 % index of second to merge
        kother = 4 - k; % the index of the non merging (only works for 3 components)
        
        % merge components
        wcomb(k) = w(i) + w(j);
        [mucomb(k), sigma2comb(k)] = reduceGaussMix(w([i; j])/wcomb(k), mus([i; j]), sigmas([i; j]).^2);
        
        % create and evaluate pdfs
        Ncomb{k} = makedist('Normal','mu',mucomb(k),'sigma',sqrt(sigma2comb(k)));
        Ncombx{k} = pdf(Ncomb{k}, x);
        pdftotcomb{k} = Nx{kother}*w(kother) + Ncombx{k}*wcomb(k);
        
        % plot
        plot(x, pdftotcomb{k}, 'DisplayName', sprintf('combining %d %d',i, j))
        k = k + 1;
    end
    
end
legend()

mucomb
sigma2comb
sigmacomb = sqrt(sigma2comb)