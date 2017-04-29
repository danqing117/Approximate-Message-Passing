function [ S, M, T, H, Z, MAI, stopnorm ] = AMP_real_tau_fun(y, x, A, N, delta, varx, I, msetol, stoptol, sigma2, x0)
        
% ------------------------ Iterative Estimator, the AMP algorithm -------------------------

x_bar = zeros(N,I+1);
n = length(y);
z = zeros(n,I+1);
stop = zeros(1,I+1);

% t=1 in fact means time 0, where h_bar(1,:)=0, z(1,:)=y, mse(0) = 1
% thresholds
tau(1) = varx;
z(:,1) = y; zsquare(1) = norm(z(:,1))^2/n;
r(:,1) = y;

% initail guess 
x_bar(:,1) = x0;
stop(1) = norm(y-A*x_bar(:,1))/ norm(y);
truemse(1) = sum(abs(x_bar(:,1)-x).^2)/N;
MAI(:,1) = (A'*A - eye(N)) * (x_bar(:,1) - x);

x_bar(:,2) = etafunreal( A'*z(:,1)+x_bar(:,1), tau(1) );
stop(2) = norm(y-A*x_bar(:,2))/ norm(y);
truemse(2) = sum(abs(x_bar(:,2)-x).^2)/N;
t = 2; 

while ((t<I+1))
    
    tau(t) = sqrt( sigma2 + (1/delta * tau(t-1) * mean( etafunprimereal(A'* z(:,t-1) + x_bar(:,t-1), tau(t-1)) ))^2);
    %tau(t) = sqrt( sigma2+(tau(t-1) * mean((A.^2) * etafunprimereal(A'* z(:,t-1) + x_bar(:,t-1), tau(t-1))))^2);
    
    z(:,t) = y - A * x_bar(:,t) + 1/delta * z(:,t-1) * mean( etafunprimereal(A'* z(:,t-1) + x_bar(:,t-1), tau(t-1)) ) ;
    %z(:,t) = y - A * x_bar(:,t) +  ((A.^2) * etafunprimereal(A'* z(:,t-1) + x_bar(:,t-1), tau(t-1))).* z(:,t-1); 
    zsquare(t) = norm(z(:,t))^2/n;
    r(t) = norm(y - A*x_bar(:,t))/n;
    
%     for a = 1:n
%         z(a,t) = y(a) - A(a,:) * x_bar(:,t) + z(a,t-1) * (A(a,:).^2) * etafunprimereal(A'* z(:,t-1) + x_bar(:,t-1), tau(t-1));
%     end


    MAI(:,t) = (A'*A - eye(N)) * (x_bar(:,t) - x);
    x_bar(:,t+1) = etafunreal( A'*z(:,t)+x_bar(:,t), tau(t) );
    truemse(t+1) = sum(abs(x_bar(:,t+1)-x).^2)/N;
    %truemse(t+1) = norm(x_bar(:,t+1) - x)/norm(x);
    truemse(t+1);
     
    if (truemse(t+1) <= msetol)
       break;
    end
   
 % stopping rule   
    stop(t+1) = norm(y-A*x_bar(:,t+1))/ norm(y);
%     if stop<1e-3
%         break;
%     end
%     
   if zsquare(t) >=2||(t>5 && var(zsquare(t-5:t))<=stoptol)
        a = find(zsquare == min(zsquare));
        if a(1)==1
            truemse(t+1) = truemse(a(1)+1);
        else
            truemse(t+1) = truemse(a(1));
        end
        break;
    end

    
    t = t+1; 
end
    
% success 
   if (truemse(length(truemse)) <= msetol)
       S = 1;
   else
       S = 0;
   end


 M=truemse;
 T = tau';
 H = x_bar;
 Z = zsquare;
stopnorm = stop(1:t-1).';



   