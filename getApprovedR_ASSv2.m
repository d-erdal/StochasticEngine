
% Function: Checks whether the rCand parameter sets are
% acceptable according to the active subspace sampler
% This is version 2 used in Erdal & Cirpka: Technical Note: Improved
% Sampling.....
% Part of: Stochastic sampling scheme version 4
%
% Created by: Daniel Erdal, Uni Tübingen, Nov. 2019
% Last edited by: N/A

% Input:
% rCand --> input matrix nrX x nrParam with candidate parameter sets
% ASS   --> active subspace proxy models, in a cell structure with one ASS
%       model per observation
% ARCH  --> Archive of not-good-enough-but-still-pretty-reasonable parameter
%       sets that could be approved with an updated proxy model. Used on if
%       the corresponding flag is turned on
% targets   --> cell structure with the observation targets that the sampling
%       is based on. One cell per observation
% nrParam   --> number of parameters
% doArchive --> use the ARCH feature
% data      --> proxyModel.data, a free format stucture for prox-model inputs

% Output:
% prBIN --> boolean vector nrX x 1 for each parameters set if it passed the
%           stage-1 test or not
% prX --> numerical "probability" vector nrX x 1 of each parameter set. 
%       No limits, but must be positive and higher must be better


function [prBIN,prX,data]=getApprovedR_ASSv2(rCand,ASS,~,data)

data.accIndex(isnan(data.Rsave(1:size(data.accIndex,1),1)),:)=[];
data.Rsave(isnan(data.Rsave(:,1)),:)=[];



% loop through all the obsevations and for each find the neighbours
% that are closest in the active subspace and check their stage-2
% acceptance
prX=zeros(size(rCand,1),length(ASS));
for j=length(ASS):-1:1
    % active subspace location for old parameter sets
    av=ASS{j}.W1'*data.Rsave';
    % active subpace location for candidate parameter sets
    avNEW=ASS{j}.W1'*rCand';
    % minimum search radius
    lim=data.ASradius*mean(max(av,[],2)-min(av,[],2));
    % compute the distance between all candidates and all old
    % parameters
    d=squeeze(sqrt(sum((repmat(avNEW(:,:),1,1,size(data.Rsave,1))-permute(repmat(av,1,1,size(rCand,1)),[1,3,2])).^2,1)));
    [dSort,ix1]=sort(d,2);
    % loop trhough each candidate set and select the closest
    % neighbours
    for i=1:size(rCand,1)
        dSortX=dSort(i,:);
        ix=ix1(i,:);
        %  [dSort,ix]=sort(d);
        dSortX(dSortX>lim)=[];
        
        ix=ix(2:max(data.nrNeighb+1,length(dSortX)));
        % compute the mean of the neighbours acceptances
        prX(i,j)=mean(data.accIndex(ix,j));
    end
end

% final evaluation statistics is based on the user set limit;
% convert prX to this true or false scale
prXlog=false(size(prX));
prXlog(prX>data.acceptLimit)=true;

% a parameter set is accepted if all observations are accepted
prBIN=all(prXlog,2);
% hence a "numeric value" is the sum over all observations
prX=sum(prXlog,2);

