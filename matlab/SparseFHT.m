function [Y, S, U, I] = SparseFHT(X, K, B, C, L, T)
% SparseFHT: syntax: [Y, S, U, I] = SparseFHT(X, K, B, C, L, T)
%
% Wrapper for the Sparse Fast Hadamard Transform
% 
% Input arguments:
% X: input vector  (size n)
% K: the sparsity (and size of y)
% B: number of buckets
% C: oversampling factor
% L: maximum number of iterations of decoder
% T: Type of algorithm to use ('Random' / 'Deterministic' / 'Optimized')
% 
% Output arguments: 
% Y: output vector (size K)
% S: support vector (size K)
% U: the number of unsatisfied checks (optional)
% I: the number of loops run (optional)

if (nargin ~= 6)
  error('Number of input arguments must be 6.');
end

if (nargout < 2 | nargout > 4)
  error('Number of output arguments must be between 2 and 4.');
end

if (strcmp(T, 'Random'))
  type = 1;
  seed = randi(2^53-1);
elseif (strcmp(T, 'Deterministic'))
  type = 2;
  seed = 0;
elseif (strcmp(T, 'Optimized'))
  type = 3;
  seed = 0;
else
  error('Error: the SparseFHT type must be ''Random'', ''Deterministic'', ''Optimized''.');
end

if (nargout == 2)
  [Y, S] = SparseFHT_mex(X, K, B, C, L, type, seed);
  U = 0;
  I = 0;
elseif (nargout == 3)
  [Y, S, U] = SparseFHT_mex(X, K, B, C, L, type, seed);
  I = 0;
elseif (nargout == 4)
  [Y, S, U, I] = SparseFHT_mex(X, K, B, C, L, type, seed);
end

