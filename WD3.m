clear;
s=[3:1:63];
successOMP=[];
successSP=[];
successIHT=[];

for t=1:61
    tic;
    S=s(t);
    success_OMP=0;
    success_SP=0;
    success_IHT=0;
    
    for tt=1:500
        
        x = randn(256,1);               %create x
        x_12 = zeros(256,1);
        A = normc(randn(128,256));      %create A
        
        x1 = randsample(x,S*1,false);
        index =zeros(S,1);
        
        
        for i=1:S
        [index(i),a]=find(x==x1(i));    %get the index of non-zero value
        end
        
        for i=1:256
            for j=1:S
                if i==index(j)
                    x_12(i)=x(i);
                end
            end
        end                             %create x_12(the sparse matrix with S dim)
        
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
      for i=1:S
    
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
    [B,I]=maxk(abs(H1_SP),S);                        %I respesent the supp index
    As_SP=[];
    s_hat=[];
    bs_SP=zeros(256,1);
    x_SP=zeros(256,1);
    err_SP=[];
    x_SPall=[];
    k=1;

    for i=1:S
        As_SP(:,i) = A(:,I(i));
    end        
    yr_SP = y-(As_SP*((As_OMP.')*As_OMP)^(-1)*(As_OMP.'))*y;
    err_SP = abs(sum(yr_SP,1));
    s_hat=[s_hat,I];

    for i=1:100
        Hs_SP = (A.')*yr_SP;                        %expand support
        [B,I]=maxk(abs(Hs_SP),S);  
        s_hat = cat(1,s_hat,I);                         %add I to index s_hat
        for j=1:(2*S)
            As_SP(:,j) = A(:,s_hat(j));              %calculate the as with 2s sprace
        end
        bs=pinv(As_SP)*y;                           %calculate bs

        for j=1:(2*S)
            bs_SP(s_hat(j))=bs(j);
        end
        [B,s_hat]=maxk(abs(bs_SP),S); 
        bs_SP=zeros(256,1);
        As_SP =[];
        for j=1:S
            As_SP(:,j) = A(:,s_hat(j));              %calculate the as with s sprace
        end
        p_SP =pinv(As_SP)*y ;
        x_SP=zeros(256,1);
        for j=1:S
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
    [B,s_IHT]=maxk(abs(x_IHT),S); 
    x_IHT1 = x_IHT;
    x_IHT = zeros(256,1);
    
    for j=1:256
        for k=1:S
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
        
    compare =[x_12,x_OMP,x_SP,x_IHT];
    
    e_OMP=sum(abs(x_OMP-x_12).^2)/sum(abs(x_12));
    e_SP=sum(abs(x_SP-x_12).^2)/sum(abs(x_12));
    e_IHT=sum(abs(x_IHT-x_12).^2)/sum(abs(x_12));       
    if e_OMP<10^(-6)
        success_OMP=success_OMP+1;
    end
    if e_SP<10^(-6)
        success_SP=success_SP+1;
    end
    if e_IHT<10^(-6)
        success_IHT=success_IHT+1;
    end
    
    end
    successOMP=[successOMP,success_OMP];
    successSP=[successSP,success_SP];
    successIHT=[successIHT,success_IHT];
    toc;
end
plot(s,successOMP/500,s,successSP/500,s,successIHT/500);
xlabel('S');
ylabel('error_rate');
title('comparison of OMP,SP,IHT');
legend('OMP','SP','IHT');