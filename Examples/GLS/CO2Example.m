load CO2.data;
plot(CO2);
n = rows(CO2);
trend = (1:n)';
x = [ones(n,1) trend];
[b, junk, e] = mc_ols(CO2,x);
plot(trend(n-36:n,:), e(n-36:n,:));
legend("residuals", "fit");
print("CO2Residuals.png","-dpng");
