# true parameters

load fa.out;
data = fa;
theta1 = data(:,3:14);
printf("Rows: %d\n", rows(data));
dstats(theta1);
theta1 = theta1(:,5:8) - theta1(:,1:4); 
bias = mean(theta1);
rmse = sqrt(mean(theta1 .^2));
printf("BIAS\n");
disp(bias);
printf("RMSE\n");
disp(rmse);

