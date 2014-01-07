% Plot the density evolution equation of Figure 5.

C = 3;
beta = [ 3 2 1 ];
pj = 0:0.01:1;

colors = [ 'rbg' ];

figure();
plot(pj, pj, 'k--');
hold on;
for i=1:length(beta(:))
  % This comes from Eq. 8 in the paper.
  % The functions lambda and rho are given in Proposition 8.
  plot(pj, (1 - exp(-beta(i)*pj)).^(C-1), colors(i));
end
hold off;

title('Figure 5. -- Density Evolution equation');
xlabel('p_j');
ylabel('p_{j+1}');
set(gca, 'XTick', 0:0.2:1);
set(gca, 'YTick', 0:0.2:1);
axis square;

leg = strvcat('p_{j+1} = p_j', [ repmat(['\beta='], length(beta(:)), 1) num2str(beta(:)) ]);
legend(leg, 'Location', 'NorthWest');
