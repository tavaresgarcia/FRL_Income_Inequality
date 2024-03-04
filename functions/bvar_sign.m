function [fsave,histmat,basemat,fevdmat,emat]=bvar(data,pattern,timemat,Lbig,REPS,BURN,HORZ,Update,maxtrys,lamdaP,tauP,epsilonP,identification,mreps,title)


Y=data;
L=Lbig;
N=cols(Y);
%take lags
X=[];
for j=1:L
    X=[X lag0(data,j) ];
end
X=[X ones(rows(X),1)];
Y=Y(L+1:end,:);
X=X(L+1:end,:);

% Additional priors for VAR coefficients
muP=mean(Y)';
sigmaP=[];
deltaP=[];
e0=[];
for i=1:N
    ytemp=Y(:,i);
    xtemp=[lag0(ytemp,1) ones(rows(ytemp),1)];
    ytemp=ytemp(2:end,:);
    xtemp=xtemp(2:end,:);
    btemp=xtemp\ytemp;
    etemp=ytemp-xtemp*btemp;
    stemp=etemp'*etemp/rows(ytemp);
    if abs(btemp(1))>1
        btemp(1)=1;
    end
    deltaP=[deltaP;btemp(1)];
    sigmaP=[sigmaP;stemp];
    e0=[e0 etemp];
end
if lamdaP~=0
    [yd,xd] = create_dummies(lamdaP,tauP,deltaP,epsilonP,L,muP,sigmaP,N);

    Y0=[Y;yd];
    X0=[X;xd];
else
    Y0=Y;
    X0=X;
end

T=rows(Y0);

fsave=zeros(REPS-BURN,N,HORZ,N);

% Display empty line (to separate samples in the screen shot)
% disp(sprintf(' '))

%conditional mean of the VAR coefficients
mstar=vec(X0\Y0);  %ols on the appended data
xx=X0'*X0;
ixx=xx\eye(cols(xx));  %inv(X0'X0) to be used later in the Gibbs sampling algorithm
sigma=eye(N); %starting value for sigma
beta0=vec(X0\Y0);
igibbs=1;
jgibbs=1;

% Burn in window
burn_w = waitbar(0,strcat(title, ' draw burn-in'));
close_burn_w = 0;
close_rep_w = 1;

while jgibbs<REPS-BURN+1

    if mod(igibbs,Update)==0 && ~close_burn_w
        waitbar(igibbs / (BURN))
    end

    % Display progress:
    if mod(igibbs,Update)==0 && (igibbs > BURN)
        % disp(sprintf(' Replication %s of %s.', ...
        %      num2str(jgibbs), num2str(REPS-BURN)) );
        waitbar(jgibbs / (REPS-BURN))
    end


    %step 1: Sample VAR coefficients
    [ beta2,PROBLEM] = getcoef( mstar,sigma,ixx,maxtrys,N,L );
    if PROBLEM
        beta2=beta0;
    else
        beta0=beta2;
    end

    %draw covariance
    e=Y0-X0*reshape(beta2,N*L+1,N);
    scale=e'*e;
    sigma=iwpq(T,inv(scale));
    A0hat=chol(sigma);

    % close burn in window
    if igibbs>BURN && ~close_burn_w
        close(burn_w);
        rep_w = waitbar(0,strcat(title, ' draws'));
        close_burn_w = 1;
    end
    if igibbs>BURN && ~PROBLEM
        beta3=reshape(beta2,N*L+1,N);
        Aols=[beta3(end,:)' beta3(1:end-1,:)'];

        datamat=data';

        %find A0 matrix
        if identification==1
            A0NEWx=A0hat;
            for jj=1:N
                shock=zeros(1,N);
                shock(jj)=1;
                irfmat(jj,:,:)=irfsim(beta2,N,L,A0NEWx,shock,HORZ+L);

            end
            fsave(jgibbs,:,:,:)=irfmat;
            %fevd
            FEVDOUT= fevd(HORZ+L,N,L,beta2,A0NEWx  );
            %hd
            IdentMat=A0NEWx';

            errormat=e(1:size(Y,1),:)*invpd(A0NEWx);
            const=1;
            DecompVAR=var_historical_decomp(Aols,IdentMat,datamat,...
                errormat',const);
            tempH=DecompVAR.yfs;
            tempB=DecompVAR.initeffecty';
            histmat(jgibbs,:,:,:)=tempH;
            basemat(jgibbs,:,:)=tempB';
            fevdmat(jgibbs,:,:)=FEVDOUT;
            jgibbs=jgibbs+1;
        else
            if mreps==1
                chck=-1;
                trys=1;
                while chck<0 && trys<maxtrys
                    K=randn(N,N);
                    Q=getqr(K);
                    A0hat1=(Q*A0hat);  %candidate draw
                    A0NEWx=geta0(A0hat1,pattern,timemat,beta2,L);

                    if sum(sum(A0NEWx))~=0
                        chck=10;
                    else
                        trys=trys+1;
                    end
                end
            else
                A0big=zeros(mreps,N*N);
                parfor m=1:mreps
                    %%
                    chckx=-1;
                    A0NEWxx = [];
                    while chckx<0
                        K=randn(N,N);
                        Q=getqr(K);
                        A0hat1=(Q*A0hat);  %candidate draw
                        A0NEWxx=geta0(A0hat1,pattern,timemat,beta2,L);

                        if sum(sum(A0NEWxx))~=0
                            chckx=10;
                        end
                    end

                    A0big(m,:)=vec(A0NEWxx)';
                end

                MA0=median(A0big); %find the median
                DIFFA0=sum((A0big-repmat(MA0,mreps,1)).^2,2); %distance from the median
                [junk,id]=min(DIFFA0);
                A0NEWx=reshape(A0big(id,:),N,N); %this is the A0 matrix closest to the median
                chck=1;
            end


            if chck>0
                irfmat=zeros(N,HORZ,N);

                %       A0NEWx'*A0NEWx-sigma
                %       A0NEWx
                for jj=1:N
                    shock=zeros(1,N);
                    shock(jj)=1;
                    irfmat(jj,:,:)=irfsim(beta2,N,L,A0NEWx,shock,HORZ+L);

                end
                fsave(jgibbs,:,:,:)=irfmat;
                %fevd
                FEVDOUT= fevd(HORZ+L,N,L,beta2,A0NEWx  );
                %hd
                IdentMat=A0NEWx';
                errormat=e(1:size(Y,1),:)*invpd(A0NEWx);
                const=1;
                DecompVAR=var_historical_decomp(Aols,IdentMat,datamat,...
                    errormat',const);
                tempH=DecompVAR.yfs;
                tempB=DecompVAR.initeffecty';
                histmat(jgibbs,:,:,:)=tempH;
                basemat(jgibbs,:,:)=tempB';
                fevdmat(jgibbs,:,:)=FEVDOUT;
                emat(jgibbs,:,:)=errormat;
                jgibbs=jgibbs+1;
            end
        end

    end
    igibbs=igibbs+1;


end
close(rep_w);

