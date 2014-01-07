function [xh] = random_k_sparse(n, k, L)

supp = randperm(n-1,k);
xh = zeros(1,n);
xh(supp) = 2*L*(rand(1,k)-0.5);

