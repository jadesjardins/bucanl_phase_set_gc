% dN is the number of shifts in each bin
% R2 is the number of shifts in previous M windows of length W
% L1 is the likelihood values
% tempM is the value of M which optimizes AIC
dN=cell(1,NumChan);
R2=dN;
LikelihoodM=zeros(NumChan,length(Mvec));
theAIC = LikelihoodM;
theM=zeros(1,NumChan);
ParametersM=cell(NumChan,length(Mvec));

% The values of M that are considered
...%Mvec=3:6;

 for MM=1:length(Mvec);
     M=Mvec(MM);
     % Calculate dN and R2 for each channel
     for p=1:NumChan
         [dN{p},R2{p}]=PSGCplugin_getPoints(Data,p,M,boundInd);
     end
     K=length(dN{1});
     % R combines an intercept term and the values of R2 from
     % each channel
     R=zeros(K,NumChan*M+1);
     R(:,1)=1;
     for ii=1:NumChan
         R(:,((ii-1)*M+2):(ii*M+1)) = R2{ii};
     end
     % Calculate likelihood and parameter estimates for each
     % channel
     for p=1:NumChan
         [LikelihoodM(p,MM),Gamma]=PSGC_noIC(dN{p},R);
         ParametersM{p,MM}=Gamma;
     end
 end
 % AIC is 2*(number of parameters) - 2*loglikelihood
 for MM=1:length(Mvec)
     M=Mvec(MM);
     theAIC(:,MM)=2*(1+NumChan*M) - 2*LikelihoodM(:,MM);
 end
 % Find M which minimizes AIC
 indexM = 1:length(Mvec);
 theMM = theM;
 for i=1:NumChan
     temp=theAIC(i,:);
     theM(i) = Mvec(temp==min(temp));
     theMM(i) = indexM(temp==min(temp));
 end
 eval(sprintf('save %s LikelihoodM theAIC theM theMM ParametersM', SaveString))
