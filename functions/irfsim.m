function y=irfsim(b,n,l,v,s,t);
%b= VAR coefs
%n=number of variables
%l=lag length
%v=A0 matrix
%s=shock vector
%t=horizon

e=zeros(t+l,n);

e(l+1,:)=s;

y=zeros(t+l,n);

for k=l+1:t;

x=[];


for i=1:l

for j=1:n

x=[x y(k-i,j)];

end
end


y(k,:)=([x 0]*reshape(b,n*l+1,n))+(e(k,:)*v);

end

 y=y(l+1:rows(y)-l,:);
 
 
 
 
 


