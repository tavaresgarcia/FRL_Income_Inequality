function [ydata, laggedydata] = makelags(ydata0,nL);

[Ty Ny]=size(ydata0);

if Ny > Ty; 
    ydata0 = ydata0';
    [Ty Ny]=size(ydata0);
end;

ydata = ydata0(nL+1:Ty,:);

for i = 1 : nL;
    laggedydata(:,(i-1)*Ny+1:i*Ny) = ydata0(nL+1-i:Ty-i,:);
end;