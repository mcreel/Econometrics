# true parameters

load ar1.out;
data = ar1;
theta = data(:,3:8);
printf("Rows: %d\n", rows(data));
dstats(theta);
theta = theta(:,3:6) - repmat(theta(:,1:2),1,2); 
bias = mean(theta);
rmse = sqrt(mean(theta .^2));
bias = reshape(bias, 2, 2);
rmse = reshape(rmse, 2, 2);
printf("BIAS\n");
disp(bias);
printf("RMSE\n");
disp(rmse);

