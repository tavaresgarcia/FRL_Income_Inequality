function xx=checkr(x,e,timei)

out=[];
n=cols(x);

for i=1:n
x1=x(:,i);
x2=e(i);
x3=timei(i)+1;
out1=0;
if x2 ~= 0
   
    x1neg=sum(x1(1:x3)<0) == x3;
    x1pos=sum(x1(1:x3)>0)==x3;
if x1neg && x2<0
    
out1=1;
end;
if x1pos && x2>0
out1=1;
end
end
out=[out out1];
end



xx=sum(out,2)==sum(abs(e),2);
