
x = randn(256,1);
A= normc(randn(128,256));
y = A*x;

x1= pinv(A)*y;
x2 = A \ y;

x_c=zeros(256,3);
x_c(:,1)=x;
x_c(:,2)=x1;
x_c(:,3)=x2;