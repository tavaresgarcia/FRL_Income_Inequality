function DecompVAR = var_historical_decomp(Aols,IdentMat,datamat,...
    errormat,const)
%%
if size(Aols,1) > size(Aols,2)
    Aols = Aols';
end
if size(datamat,1) > size(datamat,2)
    datamat = datamat';
end
[dy,pdy] = size(Aols(:,const+1:end));
p = pdy/dy;
[x,lx] = makelags(datamat',p);
% xx = [x lx(:,1:(p-1)*dy)];
T = size(x,1);
Acomp = varcompanion(Aols(:,const+1:end));
TermD = zeros(pdy,dy);
TermD(1:dy,:) = IdentMat;
Jmat = [eye(dy) zeros(dy,(p-1)*dy)];

xfd = zeros(p*dy,T+1,dy);
yfd = zeros(dy,T+1,dy);

for j = 1 : dy;
    termA  = zeros(dy,T+1);
%      size(errormat(j,:))
%      size(termA(j,2:end))

    termA(j,2:end) = errormat(j,:);
    for i = 2 : T+1;
        xfd(:,i,j) = Acomp * xfd(:,i-1,j)+ TermD * termA(:,i);
        yfd(:,i,j) = Jmat * xfd(:,i,j);
    end;
end;


initeffecty=zeros(dy, T+1);
initeffectx=zeros(p*dy, T+1);
initeffectx(:,1)=lx(1,:)';
initeffecty(:,1)=Jmat*initeffectx(:,1);
if const==1;
    CC=zeros(p*dy,1);
    CC(1:dy,:)=Aols(:,1);
else
    CC=zeros(p*dy,1);
end
for i = 2 : T+1;
    initeffectx(:,i)=CC+Acomp*initeffectx(:,i-1);
    initeffecty(:,i) = Jmat * initeffectx(:,i);
end;

yf = initeffecty + sum(yfd,3);

yfs = zeros(T+1,dy,dy);
for i = 1 : dy;
    for j = 1 : dy;
        yfs(:,j,i)=yfd(i,:,j)';
    end;
end;
DecompVAR.yfd = yfd;
DecompVAR.yfs = yfs; %yfs(:,:,1) will be T by nshocks 
DecompVAR.yf = yf;
DecompVAR.initeffecty = initeffecty;

