% deltas = [0,10,50,80,100,300,500,700,1000];
% ginis = [0.194, 0.197, 0.209, 0.209, 0.218, 0.261, 0.294, 0.282, 0.282];
% sums = [6720, 7792, 9292, 9863,10240, 10864, 11112, 11140, 11140];
% 
% figure (1)
% yyaxis left
% plot(deltas, -ginis)
% xlabel('Deltas')
% ylabel('-Gini coeffs')
% yyaxis right
% plot(deltas, sums)
% ylabel('Sum of community flow')

%%
deltas = linspace(0, 500, 50);
% l = length(deltas);
% ginis = zeros(1, l);
% sums = zeros(1,l);
% 
% 
% for i = 1:l
%     [gini, sum] = vertiportselecfunc(deltas(i), 17, 0, 0);
%     ginis(i) = gini;
%     sums(i) = sum;
% end

load('ginis.mat');
load('sums.mat');

f = figure ('Position', [50 50 950 1000]);
ax = gca;
set(ax, 'FontName', 'Times', 'FontSize', 22);
ax.LineWidth = 1.5;
yyaxis left
plot(deltas, -ginis, 'linewidth', 2)
title('Fairness/Efficiency Tradeoff')
xlabel('Delta')
ylabel('-Gini Coefficients')
yyaxis right
plot(deltas, sums, 'linewidth', 2)
ylabel('Sum of Community Flow')