function [ FEVDOUT] = FEVD(HORIZON,N,L,betam,A0NEWm  )
FEVDOUT=0;
   
        
    %forecast error descomposition
    yhatall=zeros(HORIZON,N*N);
    totvar=0;
    jjj=1;
   for iii=1:N
        yhat=zeros(HORIZON+L,N);
        vhat=zeros(HORIZON+L,N);
        vhat(L+1,iii)=1; 

for j=L+1:HORIZON+L
    xhat=[];
    for jj=1:L
        xhat=[xhat yhat(j-jj,:)];
    end
    xhat=[xhat 0];
 yhat(j,:)=xhat*reshape(betam,N*L+1,N)+vhat(j,:)*A0NEWm;
end
yhatall(:,jjj:jjj+N-1)=yhat(L+1:end,:);
totvar=totvar+cumsum(yhat(L+1:end,:).^2);
jjj=jjj+N;
   end
    
  FEVDOUT=cumsum(yhatall.^2)./repmat(totvar,1,N);  

end

