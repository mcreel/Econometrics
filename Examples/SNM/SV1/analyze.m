e# true parameters

load sv.out;
data = sv;
theta1 = data(:,3:11);
theta1(:,2) = 2*slog(theta1(:,2));
theta1(:,5) = 2*slog(theta1(:,5));
theta1(:,8) = 2*slog(theta1(:,8));
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

