
clearvars -except 'Repetitions' 'Type';

addpath('../C');

if (~exist('Repetitions'))
  Repetitions = 1000;  % repetition
end
if (~exist('Type'))
  Type = 'Deterministic';
end

% Global parameters
n = 22;        % log2 of signal size
N = 2^n;       % signal size
L = 20;        % Maximum number of loop

%%%%%%%%%%%%%%%%%
% SELECT HASHES %
%%%%%%%%%%%%%%%%%

b = 1:n-2;
C = 1:12;

blen = length(b);
clen = length(C);

% initialize RNG
seed = sum(100*clock); % randomness!
rng(seed, 'twister');

% Arrays
Error     = zeros(blen, clen, Repetitions);
BitError  = zeros(blen, clen, Repetitions);
Unsat     = zeros(blen, clen, Repetitions);
Loops     = zeros(blen, clen, Repetitions);
SuppSize  = zeros(blen, clen, Repetitions);
Success   = zeros(blen, clen, Repetitions);

%%%%%%%%%%%%%%%%%%
% RUN EXPERIMENT %
%%%%%%%%%%%%%%%%%%

if (matlabpool('size') == 0)
  matlabpool open;
end

for ic = 1:clen
  fprintf('%d ', C(ic));
  parfor i=1:Repetitions
    for ib = 1:blen

      % current parameters
      B = 2^b(ib);
      K = 2^b(ib);

      % random sparse spectrum (magnitude zero mean unit variance normal)
      xh = randn_k_sparse(N, K, 100);

      % compute time-domain (TD) signal
      x = FastHadamard(xh')/sqrt(N);

      % SFHT!
      [zhat, zsup, us, lp] = SparseFHT(x, K, B, C(ic), L, Type);

      % Accumulate error
      I = find(zsup >= 0);

      e = zeros(size(xh));
      e(zsup(I)+1) = zhat(I);
      Error(ib, ic, i)     = sum((xh - e).^2);

      SuppSize(ib, ic, i)  = length(I);
      BitError(ib, ic, i)  = numel(setdiff(zsup(I)+1, find(xh)));
      Unsat(ib, ic, i)     = us;
      Loops(ib, ic, i)     = lp;
      Success(ib, ic, i)   = (us == 0);

    end
  end
end
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%
% DO SMTH WITH RESULTS %
%%%%%%%%%%%%%%%%%%%%%%%%

% save the data
save('ErrorSim.mat');

