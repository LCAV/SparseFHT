function [m ci] = median_and_ci(X);
% computes median and 0.95% confidence interval.

X = sort(X);
n = length(X);

if mod(n,2) == 1
  % if n is odd, take central element
  m = X((n+1)/2);
else
  % if n is even, average the two central elements
  m = 0.5*(X(n/2) + X(n/2+1));
end

% This table is taken from the Performance Evaluation lecture notes by J-Y Le Boudec
% available at: http://perfeval.epfl.ch/lectureNotes.htm
CI = [  1  6 ;   1 7 ;  1  7 ;  2  8 ;  2  9 ;  2 10 ;  3 10 ;  3 11 ;  3 11 ; 4 12 ; ...
        4 12 ;  5 13 ;  5 14 ;  5 15 ;  6 15 ;  6 16 ;  6 16 ;  7 17 ;  7 17 ; 8 18 ; ...
        8 19 ;  8 20 ;  9 20 ;  9 21 ; 10 21 ; 10 22 ; 10 22 ; 11 23 ; 11 23 ; ...
       12 24 ; 12 24 ; 13 25 ; 13 26 ; 13 27 ; 14 27 ; 14 28 ; 15 28 ; 15 29 ; ...
       16 29 ; 16 30 ; 16 30 ; 17 31 ; 17 31 ; 18 32 ; 18 32 ; 19 33 ; 19 34 ; 19 35 ; ...
       20 35 ; 20 36 ; 21 36 ; 21 37 ; 22 37 ; 22 38 ; 23 39 ; 23 39 ; 24 40 ; ...
       24 40 ; 24 40 ; 25 41 ; 25 41 ; 26 42 ; 26 43 ; 26 44 ; 27 44 ];

if (n < 6)
  % If we have less than 6 samples, we cannot have a confidence interval
  ci = [0 0];
elseif (n <= 70)
  % For 6 <= n <= 70, we use exact values from the table
  j = CI(n-5,1);
  k = CI(n-5,2);
  ci = [X(j)-m X(k)-m];
else
  % For 70 < n, we use the approximation for large sets
  j = floor(0.5*n - 0.98*sqrt(n));
  k = ceil(0.5*n + 1 + 0.98*sqrt(n));
  ci = [X(j)-m X(k)-m];
end
