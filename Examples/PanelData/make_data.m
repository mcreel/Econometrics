n = 20;
T = 10;
k = 2;
[y, x] = dgp(n,T,k);
data = [y x];
save -ascii 'data.txt' data;
