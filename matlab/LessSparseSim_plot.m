
figure();

beta = round(2.^(alpha*n))/B;

figure();
plot(beta, mean(Success, 2));

xlabel('\beta');
title(['Probability of success -- N=2^{' num2str(n) '}']);
set(gca, 'YGrid', 'on');
set(gca, 'XTick', 1:4);
set(gca, 'YTick', 0:0.2:1);
axis([1-0.01 4+0.01 0-0.01 1+0.01]);
