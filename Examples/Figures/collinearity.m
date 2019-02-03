# this program illustrates the effect of collinearity on the 
# least squares objective function


n=1000;
x=randn(n,2);
c=[1 0.9;0.9 1];
p=chol(c);
x2=x*p;   % the correlated regressors

y=randn(n,1);  % true coefficients are zero
b1=-5:0.2:5;
b2=b1;
p = columns(b1);
e=zeros(p,p);
e2=e;
for i=1:p
   for j=1:p
      ee=y-x(:,1)*b1(i)-x(:,2)*b2(j);
      e(i,j)=ee'*ee;
      ee=y-x2(:,1)*b1(i)-x2(:,2)*b2(j);
      e2(i,j)=ee'*ee;

   end
end

figure(1);
contour(b1,b2,e/n+10);
print("nocollin.eps", "-depsc2");

figure(2);
contour(b1,b2,e2/n+10);
print("collin.eps", "-depsc2");



