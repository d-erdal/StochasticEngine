% Function: Checks whether the rCand parameter sets are
% acceptable according to the GPE-proxy model
% This is version 1 used in Erdal, Xiao, Nowak & Cirpka 2020 .......
% Part of: Stochastic sampling scheme version 4
%
% Created by: Daniel Erdal, Uni Tübingen, Nov. 2019
% Last edited by: N/A

% Input:
% rCand --> input matrix nrX x nrParam with candidate parameter sets
% GPE   --> Gaussian process emulator models, in a cell structure with one GPE
%       model per observation
% targets   --> cell structure with the observation targets that the sampling
%       is based on. One cell per observation
% data      --> proxyModel.data, a free format stucture for prox-model inputs

% Output:
% prBIN --> boolean vector nrX x 1 for each parameters set if it passed the
%           stage-1 test or not
% prX --> numerical probability vector nrX x 1 of each parameter set.
%       Contains the actual probabilities


function [prBIN,prX,data]=getApprovedR_GPE(rCand,GPE,targets,data)



% check if it is accepted accodring the the GPE and the targets
%%
prEachObs=zeros(size(rCand,1),length(targets),1);
for i=1:length(GPE)
    % get the GPE-predicted observation from the STK
    obsPred = stk_predict (GPE{i}.M_post , rCand);
    
    % Compute the probability of this parameter set being on the right side
    % of its target
    prEachObs(:,i)=normcdf(targets{i}.targetTRUE,obsPred.mean,sqrt(obsPred.var));
    
    if ~targets{i}.targetBELOW
        % if we are looking for values below the target: invert
        prEachObs(:,i)=1-prEachObs(:,i);
    end
    
end
%%

% Compute the joint probability of all observations
prX=prod(prEachObs,2);
% and choose to stage-1 accept or reject
if isnan(data.probAccept)
    % if nan, acceptance should be at random
    prBIN=prX>rand(size(prJoint));
else
    % otherwise it should be a given limit
    prBIN=prX>data.probAccept;
end











