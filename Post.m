figure;
hold on;
%Ñ¡Ò»¸ö
% plot(QSet*1e5,results(:,1:3)*1e3);
plot(QSet*1e5,results*1e3);
exp=csvread('ExpInKate.csv',1,0);
scatter(exp(:,1),exp(:,2));
scatter(exp(:,4),exp(:,5));
legend();
%legend('Kate-exp','Arak-exp','location','northwest');
xlabel('Q(*1e-5 m3/h)');
ylabel('R_j(mm)');
title('R_j  Q');
grid;