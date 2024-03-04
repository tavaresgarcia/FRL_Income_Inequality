function a0new=geta0(a0old,pattern,timemat,beta,L)

n=rows(a0old);
out=[];

tempid= sum(abs(pattern),2)>0;
patterni=pattern(tempid,:);
timei=timemat(tempid,:);

nshocks=rows(patterni);



for MM=1:nshocks

patternj=patterni(MM,:); %signs to check
timej=timei(MM,:); %horizon to check for
maxtime=max(timej)+1;


zmp=-999;

for j=1:n
z1=a0old(j,:);
%impulse response
shock=zeros(1,n);
shock(j)=1;
irf=irfsim(beta,n,L,a0old,shock,maxtime+L);


em1=checkr(irf,patternj,timej);

if em1==1;

   zmp=1*z1;
else

em1=checkr(irf*-1,patternj,timej);
if em1==1;
   zmp=-1*z1;
end;
end

end
 
 if zmp==-999;
out=[out;nan(1,n)];
 else
out=[out;zmp];
 end

end


if sum(sum(isnan(out)))>0;
a0new=zeros(n,n);
else
% out
    outx=[];
	  
      for i=1:n
	   zA=a0old(i,:);
eam=[];
	  
       for j=1:rows(out)
           tempcheck=(sum(zA == out(j,:)) || sum(-zA == out(j,:)));
	   eam=[eam tempcheck];
       end
	   eALL=sum(eam,2);
	     if eALL==0
	     outx=[outx;zA];
         end
      end
%find appropriate row to insert shocks

a0new=zeros(n,n);
i=1;
j=1;
jj=1;
for i=1:n
    tempcheck=sum(abs(pattern(i,:)),2)==0;
    if tempcheck==1
a0new(i,1:n)=outx(j,:);
j=j+1;
    else
a0new(i,1:n)=out(jj,:);
jj=jj+1;
    end
end

end
	 