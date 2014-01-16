
clearvars -except 'Repetitions' 'Type';

addpath('../C');

if (~exist('Repetitions'))
  Repetitions = 1000;
end
if (~exist('Type'))
  Type = 'Deterministic';
end

% Global parameters
b = 17;
n = 22;
C = 4;
beta = 1:0.25:4;
alpha = (log2(beta)+b)/n;
L = 20;
N = 2^n;
B = 2^b;

alen = length(alpha);

% initialize RNG
seed = sum(100*clock); % randomness!
rng(seed, 'twister');

% Arrays
Error     = zeros(alen, Repetitions);
BitError  = zeros(alen, Repetitions);
Unsat     = zeros(alen, Repetitions);
Loops     = zeros(alen, Repetitions);
SuppSize  = zeros(alen, Repetitions);
Success   = zeros(alen, Repetitions);

%%%%%%%%%%%%%%%%%%
% RUN EXPERIMENT %
%%%%%%%%%%%%%%%%%%

if (matlabpool('size') == 0)
  matlabpool open;
end

parfor ia=1:alen

  % number of non-zero coefficients
  K = round(2^(alpha(ia)*n));

  for i=1:Repetitions

    % random sparse spectrum (magnitude zero mean unit variance normal)
    xh = randn_k_sparse(N, K, 100);

    % compute time-domain (TD) signal
    x = FastHadamard(xh')/sqrt(N);


    % SFHT!
    [zhat, zsup, us, lp] = SparseFHT(x, K, B, C, L, Type);

    % Accumulate error
    I = find(zsup >= 0);

    e = zeros(size(xh));
    e(zsup(I)+1) = zhat(I);
    Error(ia, i)     = sum((xh - e).^2);

    SuppSize(ia, i)  = length(I);
    BitError(ia, i)  = numel(setdiff(zsup(I)+1, find(xh)));
    Unsat(ia, i)     = us;
    Loops(ia, i)     = lp;
    Success(ia, i)   = (us == 0);

  end
end

%%%%%%%%%%%%%%%%%%%%%%%%
% DO SMTH WITH RESULTS %
%%%%%%%%%%%%%%%%%%%%%%%%

% save the data
save('LessSparseSim.mat');

