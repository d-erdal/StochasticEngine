% Function: trains a GPE-proxy model and uses a monte carlo (brute force)
% method to find a highly uncertain parameter set to be added to the
% training data.
% This is version 1 used in Erdal, Xiao, Nowak & Cirpka 2020 .......
% Part of: Stochastic sampling scheme version 4
%
% Created by: Daniel Erdal, Uni Tübingen, Nov. 2019
% Last edited by: N/A

% Input:
% R         --> parameter matrix (scaled 0-1), nrX x nrParam
% obsMTRX 	--> observation matrix, nrX x nrObs
% targets   --> cell structure with the observation targets that the sampling
%       is based on. One cell per observation
% nrParam   --> number of parameters
% proxyData --> proxyModel.data, a free format stucture for prox-model inputs
% doParalleTraining --> should the training be done in parallel

% Output:
% newR      --> new parameter set to be evaulated with the full model
% obsExp    --> expected value of the observation accrding to the proxy
% doALFinished  --> if true, the active learning is ending now
% proxyData     --> same as input


% [newR,obsExp,proxySetup.doALFinished,proxySetup.data]=proxySetup.ActiveLearning.run(Rsave,obsMTRX,modelSetup.targets,modelSetup.nrParam,proxySetup.data);
function  [newR,obsExp,doALFinished,proxyData]=getApprovedR_activeLearning_GPE(R,obsMTRX,targets,nrParam,proxyData,doParalleTraining)

% sort out the matrixes and a counter
R(isnan(obsMTRX(:,1)),:)=[];
obsMTRX(isnan(obsMTRX(:,1)),:)=[];
nrObs=size(obsMTRX,2);


% Setup some training information
interV_TMP=proxyData.estimateGPEparameterInterval;
proxyData.estimateGPEparameterInterval=proxyData.AL_estimateGPEparameterInterval;

% Initilize and train the GPE based on the current set of R and obs
GPE=cell(nrObs,1); data_TMP=cell(nrObs,1);
if doParalleTraining
parfor i=1:nrObs
    [GPE{i},data_TMP{i}]=trainGPE(R,obsMTRX(:,i),proxyData,[],i);
end
for i=1:nrObs
    proxyData.paramX{i}=data_TMP{i}.paramX{i};
    proxyData.trainParamCounter{i}=data_TMP{i}.trainParamCounter{i};
end
else
for i=1:nrObs
    [GPE{i},proxyData]=trainGPE(R,obsMTRX(:,i),proxyData,[],i);
end
end
% Reset the esitmation iterval to the original
proxyData.estimateGPEparameterInterval=interV_TMP;
%% Find a highly uncertain parameter set
% Use a Monte Carlo style approach with a large sample

% do nrL repetitions of 1,000 simulations each 
nrL=round(proxyData.AL_MCsampleSize/1e3);
st=zeros(nrL,nrParam+1);
prJsave=zeros(nrL,1e3);
RX=rand(1e3,nrParam,nrL);
for j=1:nrL   
    rCand=RX(:,:,j);
    prX=zeros(size(rCand,1),length(targets));
    for i=1:nrObs 
        % do a precition with the STK toolbox
        obsPred = stk_predict (GPE{i}.M_post ,rCand);
        % Compute the probability of exceedance
        prX(:,i)=normcdf(targets{i}.targetTRUE,obsPred.mean,sqrt(obsPred.var));
        
        if ~targets{i}.targetBELOW
            prX(:,i)=1-prX(:,i);
        end
    end
    % compute the joint probabilities
    prJoint=prod(prX,2);
    % and save
    prJsave(j,:)=prJoint;
    % compute the missclassification (0.5 is the worst we can get, then we do not 
    % know at all if we are behavioral or non-bahvioral)
    pr05=abs(0.5-prJoint);
    
    % and store in the st-matrx
    st(j,:)=[rCand(min(pr05)==pr05,:),pr05(min(pr05)==pr05)];
    
end

% now find the minumim; i.e. the point closest to 0.5, i.e. the most
% uncertain point of them all
newR=st(min(st(:,end))==st(:,end),1:nrParam);

% get the expected value for each of the targets for the new R
obsExp=zeros(length(targets),1);
for i=1:nrObs
     obsPred = stk_predict (GPE{i}.M_post ,newR);
     obsExp(i)=obsPred.mean;
end

%% judge if the booster is done with its work.....
% This if currently rather ad-hoc as all more appropriate statistics never
% seemed to converge for the tested cases

doALFinished=false;
Q=reshape(prJsave,numel(prJsave),1).*(1-reshape(prJsave,numel(prJsave),1));
disp(['    Current confusion: ' num2str(max(Q))])
if max(Q) < proxyData.AL_maxConfusion ||size(R,1) > proxyData.AL_maxDataSet || proxyData.STaccRate > proxyData.AL_maxShortTermAccRate
    doALFinished=true;
    
    disp('Active learning is finishing.')
    if max(Q) < proxyData.AL_maxConfusion
        disp(['  --> confusion target ' num2str(proxyData.AL_maxConfusion) ' reached'])
    end
    if size(R,1) > proxyData.AL_maxDataSet
        disp(['  --> max dataset size target ' num2str(proxyData.AL_maxDataSet) ' reached'])
    end
    if proxyData.STaccRate > proxyData.AL_maxShortTermAccRate
      disp(['  --> short term acceptance target ' num2str(proxyData.AL_maxShortTermAccRate) ' reached'])  
    end
    
end



