
alpha = b/n;

Ps = mean(Success, 3)';

figure();
imagesc(flipud(Ps));
%colormap('pink');
xlabel('\alpha');
ylabel('C');
title(['Probability of success - N=2^{' num2str(n) '} - Algo: ' Type]);
set(gca,'XTick', n/3*[1 2]);
set(gca, 'XTickLabel', ['1/3' ; '2/3']);
cticks = sort(clen - (3:3:clen) + 1);
set(gca,'YTick', cticks);
set(gca, 'YTickLabel', clen - cticks + 1);

Cth = zeros(1,blen);
Cth(1:floor(n/3)) = n./b(b <= n/3);
Cth(ceil(n/3):floor(2*n/3)) = 3;
Cth(ceil(2*n/3):end) = n./(n - b(b >= ceil(2*n/3)));

hold on;
H = plot(b, clen - Cth + 1, 'k');
set(H, 'LineWidth', 2);
hold off;

colormap('Gray');
colorbar;
