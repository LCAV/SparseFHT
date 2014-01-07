
clear;

addpath('../C');

% parameters
n = 5:23;
L = 10;
Loops  = 1000;
Subloop = 10;
Warm   = 2;
Body   = 10;
MaxMag = 500;

% Output array
b = cell(1, length(n));
C = cell(1, length(n));
Tsfht = cell(1, length(n));
Tfht = cell(1, length(n));
intersection = zeros(size(n));

% randomness!
seed = 0; % RNG seed (fixing the seed fixes the sparsity pattern)
seed = sum(100*clock); % randomness!
rng(seed, 'twister');

for k=1:Loops/Subloop
  fprintf('Subloop %d/%d : ', k, Loops/Subloop);
  for i=1:length(n)
    fprintf('%d ', n(i));
    b{i} = 1:n(i)-2;
    N = 2^n(i);
    K = 2.^b{i};
    B = K;
    R = [Subloop Warm Body MaxMag];

    C{i} = zeros(size(b{i}));
    C{i}(1:floor(n(i)/3)) = n(i)./(b{i}(1:floor(n(i)/3)));
    C{i}(ceil(n(i)/3):floor(2*n(i)/3)) = 3;
    C{i}(ceil(2*n(i)/3):end) = ceil(n(i)./(n(i) - b{i}(ceil(2*n(i)/3):end)));

    % experiment!
    [tfht tsfht] = HadamardBenchmark(N, K, B, C{i}, L, R, randi(2^53-1));

    % Store in cell array (unit is in ms)
    Tfht{i} = [Tfht{i} ; tfht];
    Tsfht{i} = [Tsfht{i} ; tsfht];

  end
  fprintf('\n');
end

% save the file.
save('TimingSim.mat');

