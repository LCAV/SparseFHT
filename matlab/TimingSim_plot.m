
% Plot parameters
nplot = 15;
in = find(n == 15);

% plot with error bars
m_fht = zeros(1, size(Tfht{in}, 2));
ci_fht = zeros(2, size(Tfht{in}, 2));
m_sfht = zeros(1, size(Tsfht{in}, 2));
ci_sfht = zeros(2, size(Tsfht{in}, 2));

% Compute median and confidence intervals on it
[m_fht ci_fht] = median_and_ci(Tfht{in});
for i=1:size(Tsfht{in}, 2);
  [m_sfht(i) ci_sfht(:,i)] = median_and_ci(Tsfht{in}(:,i));
end

figure();
%errorbar(b{in}/n(in), m_fht*ones(size(b{in})), ci_fht(1)*ones(size(b{in})), ci_fht(2)*ones(size(b{in})), 'b-');
plot(b{in}/n(in), m_fht*ones(size(b{in})), 'b-');
hold on;
%errorbar(b{in}/n(in), m_sfht, ci_sfht(1,:), ci_sfht(2,:), 'r-');
plot(b{in}/n(in), m_sfht, 'r.-');
hold off;
title(['Median runtime in [ms] -- N=2^{' num2str(n(in)) '}']);
xlabel('\alpha');
set(gca, 'YGrid','on');
set(gca,'XTick', [1 2]/3);
set(gca, 'XTickLabel', ['1/3' ; '2/3']);
axis([b{in}(1)/n(in)-0.01 b{in}(end)/n(in)+0.01 0-0.01 1.5+0.01]);

% Plot the largest b such that sfht is faster than fht
for i=1:length(n)
  mtfht = median(Tfht{i}, 1);
  mtsfht = median(Tsfht{i}, 1);
  I = find(diff(sign(mtfht-mtsfht)) ~= 0);
  if (numel(I) == 0 & mtfht(1) > mtsfht(1))
    intersection(i) = b{i}(end);
  elseif (mtfht(1) < mtsfht(1))
    intersection(i) = 0;
  else
    intersection(i) = min(b{i}(I));
  end
end
figure();
plot(n, intersection./n, 'b.-');
set(gca, 'YGrid','on');
xlabel('n = log_2(N)');
ylabel('\alpha^*')
title('\alpha^* = max(\alpha : T_{fht}(n) > T_{sfht}(\alpha'',n), \forall \alpha''\leq \alpha)');
axis([n(1)-0.1 n(end)+0.1 0-0.01 0.8+0.01]);
