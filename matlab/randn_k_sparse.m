function [xh] = random_k_sparse(n, k, sigma2)

supp = randperm(n-1,k);
xh = zeros(1,n);
xh(supp) = sqrt(sigma2)*randn(1,k);

