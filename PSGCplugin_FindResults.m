
% Output variables
ReducedLikelihood=zeros(NumChan);
GammaEstimates=cell(1,NumChan);

% Pick out the likelihood and parameter estimates from the calculations
% done in PSGC_FindModel

FullLike=zeros(1,NumChan);
for jj=1:NumChan
    FullLike(jj)=LikelihoodM(jj,theMM(jj));
    GammaEstimates{jj}=ParametersM{jj,theMM(jj)};
end
FullLikelihood=repmat(FullLike, NumChan, 1);


% Temporary variables for each participant
dN=cell(1,NumChan);
R2=dN;

           

  
for MM=1:length(Mvec)
    M=Mvec(MM);
    % If any channels have this value of M, then find dN and R2
    if sum(theM==M)>0
        for p=1:NumChan
            [dN{p},R2{p}]=PSGCplugin_getPoints(Data,p,M,boundInd);
        end
        K=length(dN{1});
        % R combines an intercept term and the value of R2 from each
        % channel
        R=zeros(K,1+NumChan*M);
        R(:,1)=1;
        for ii=1:NumChan
            R(:,((ii-1)*M+2):(1+ii*M)) = R2{ii};
        end
        for Chan=1:NumChan
            % If the channel has order M, then compute the reduced
            % likelihood for each other channel
            if theM(Chan)==M
                Chan
                for Chan2=1:NumChan
                    if Chan~=Chan2
                        % R4 is R with the contribution of Chan2
                        % removed
                        R4=R;
                        R4(:,((Chan2-1)*M+2):(Chan2*M+1))='';
                        % x0 is the initial conditions for the
                        % algorithm based on results from
                        % PSGC_FindModel
                        x0 = GammaEstimates{Chan};
                        x0(((Chan2-1)*M+2):(Chan2*M+1))='';
                        [ReducedLikelihood(Chan2,Chan),Gamma]=PSGC(dN{Chan},R4,x0);
                    else
                        % No analysis when Chan=Chan2
                        ReducedLikelihood(Chan2,Chan)=0;
                    end
                end
                eval(sprintf('save %s LikelihoodM theAIC theM ParametersM FullLikelihood ReducedLikelihood GammaEstimates', SaveString))
                PSGCi=1;
                if isfield(EEG,'PSGC');
                    PSGCi=length(EEG.PSGC)+1;
                end
                EEG.PSGC(PSGCi).label=SaveString;
                EEG.PSGC(PSGCi).LikelihoodM=LikelihoodM;
                EEG.PSGC(PSGCi).AIC=theAIC;
                EEG.PSGC(PSGCi).M=theM;
                EEG.PSGC(PSGCi).ParamM=ParametersM;
                EEG.PSGC(PSGCi).fullLikelihood=FullLikelihood;
                EEG.PSGC(PSGCi).ReducedLlihood=ReducedLlihood;
                EEG.PSGC(PSGCi).GammaEstimates=GammaEstimates;
                    
            end
        end
    end
end

Statistic=-2*(ReducedLikelihood-FullLikelihood);
eval(sprintf('save %s LikelihoodM theAIC theM ParametersM FullLikelihood ReducedLikelihood GammaEstimates Statistic', SaveString))


