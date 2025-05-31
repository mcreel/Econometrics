# true parameters

load art.out;
data = art;
theta1 = data(:,3:11);
printf("Rows: %d\n", rows(data));
dstats(theta1);
theta1 = theta1(:,4:9) - repmat(theta1(:,1:3),1,2); 
bias = mean(theta1);
rmse = sqrt(mean(theta1 .^2));
bias = reshape(bias, 3, 2);
rmse = reshape(rmse, 3, 2);
printf("BIAS\n");
disp(bias);
printf("RMSE\n");
disp(rmse);

