% exercise 2

m = 128;
n = 256;
S = 12;

x = randn(256,1);
x_12 = zeros(256,1);
A = normc(randn(128,256));

x1 = randsample(x,12*1,false);
index =zeros(12,1);

for i=1:12
    [index(i),a]=find(x==x1(i));
end

for j=1:256
   switch j
        case index(1)
           x_12(j)=x(j);
        case index(2)
           x_12(j)=x(j);
        case index(3)
           x_12(j)=x(j);
        case index(4)
           x_12(j)=x(j);
        case index(5)
           x_12(j)=x(j);
        case index(6)
           x_12(j)=x(j);
        case index(7)
           x_12(j)=x(j);
        case index(8)
           x_12(j)=x(j);
        case index(9)
           x_12(j)=x(j);
        case index(10)
           x_12(j)=x(j);
        case index(11)
           x_12(j)=x(j);
        case index(12)
           x_12(j)=x(j);
       otherwise
           x_12(j)=0;
   end
end

y=A*x_12;
x1= pinv(A)*y;
x2 = A \ y;

%------greedy algorithms£¨OMP£©--------------------------------------------
S_OMP=[];
x_OMP=zeros(256,1);
yr_OMP=y;
As_OMP=[];
x_OMP1=[];
err_OMP=[];
for i=1:12
    H1_OMP = (A.')*yr_OMP;
    [z,supp_OMP] = max(abs(H1_OMP),[],1);           %find index supp
    S_OMP = [S_OMP,supp_OMP];                       %add index
    As_OMP(:,i) = A(:,supp_OMP);                    %get the value of as 
    x_OMP1= ((As_OMP.')*As_OMP)^(-1)*(As_OMP.')*y;  %calculate the x
    for j=1:i
        x_OMP(S_OMP(j),:)=x_OMP1(j,:);              %set other elements 0
    end
    err_OMP =[err_OMP,sum((y-A*x_OMP),1)] ;
    yr_OMP=y-A*x_OMP;
end
%---------greedy algorithms (SP)------------------------------------------
H1_SP=(A.')*y;
[B,I]=maxk(abs(H1_SP),12);                        %I respesent the supp index
As_SP=[];
s_hat=[];
bs_SP=zeros(256,1);
x_SP=zeros(256,1);
err_SP=[];
x_SPall=[];
k=1;

for i=1:12
    As_SP(:,i) = A(:,I(i));
end
%yr_SP = y-(As_SP*pinv(As_SP))*y;         
yr_SP = y-(As_SP*((As_OMP.')*As_OMP)^(-1)*(As_OMP.'))*y;
err_SP = abs(sum(yr_SP,1));
s_hat=[s_hat,I];

for i=1:50
    Hs_SP = (A.')*yr_SP;                        %expand support
    [B,I]=maxk(abs(Hs_SP),12);  
    s_hat = cat(1,s_hat,I);                         %add I to index s_hat
    for j=1:24
        As_SP(:,j) = A(:,s_hat(j));              %calculate the as with 2s sprace
    end
    bs=pinv(As_SP)*y;                           %calculate bs
    for j=1:24
        bs_SP(s_hat(j))=bs(j);
    end
    [B,s_hat]=maxk(abs(bs_SP),12); 
    bs_SP=zeros(256,1);
    As_SP =[];
    for j=1:12
        As_SP(:,j) = A(:,s_hat(j));              %calculate the as with s sprace
    end
    p_SP =pinv(As_SP)*y ;
    x_SP=zeros(256,1);
    for j=1:12
        x_SP(s_hat(j))=p_SP(j);
    end
    x_SPall=[x_SPall,x_SP];
    if i==1
        err_SP =cat(1,1000,err_SP);
        k=k+1;
    end
    err_SP =cat(1,err_SP,sum((y-A*x_SP),1));
    yr_SP=y-A*x_SP;
    if err_SP(k)>err_SP(k-1) 
        break
    end
    k=k+1;
end
%-----------greedy algorithms(IHT)----------------------------------------
x_IHT=zeros(256,1);
err_IHT = [];
k1=1;

for i=1:100
    y_IHT = (A.')*(y-A*x_IHT);
    x_IHT=x_IHT + y_IHT;
    [B,s_IHT]=maxk(abs(x_IHT),12); 
    x_IHT1 = x_IHT;
    x_IHT = zeros(256,1);
    
    for j=1:256
        for k=1:12
            if j==s_IHT(k)
                x_IHT(j)=x_IHT1(j);
            end
        end
    end
    
    if i==1
        err_IHT =cat(1,1000,err_IHT);
        k1=k1+1;
    end
    err_IHT =cat(1,err_IHT,sum((y-A*x_IHT),1));
    if err_IHT(k1)>err_IHT(k1-1) 
        break
    end
    k1=k1+1;
    
end
compare =[x_12,x1,x2,x_OMP,x_SP,x_IHT];
compare1=[find(x_12~=0),find(x_OMP~=0),find(x_SP~=0),find(x_IHT~=0)];

