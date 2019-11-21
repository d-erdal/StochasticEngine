
% Function: Checks whether the rCand parameter sets are
% acceptable according to the active subspace sampler
% This is version 1 used in Erdal & Cirpka: Global sensitivity analysis and adaptive stochastic sampling of a subsurface-flow model using active subspaces
% Part of: Stochastic sampling scheme version 4
%
% Created by: Daniel Erdal, Uni Tübingen, Nov. 2019
% Last edited by: N/A

% Input:
% rCand --> input matrix nrX x nrParam with candidate parameter sets
% ASS   --> active subspace proxy models, in a cell structure with one ASS
%       model per observation
% targets   --> cell structure with the observation targets that the sampling
%       is based on. One cell per observation
% data      --> proxyModel.data, a free format stucture for prox-model inputs

% Output:
% prBIN --> boolean vector nrX x 1 for each parameters set if it passed the
%           stage-1 test or not
% prX --> numerical "probability" vector nrX x 1 of each parameter set. 
%       No limits, but must be positive and higher must be better


% [prBIN,prX,data]=proxySetup.getApprovedR(rCand,ProxyModel,targets,proxySetup.data);
function [prBIN,prX,data]=getApprovedR_ASSv1(rCand,ASS,targets,data)

% check if it is accepted accodring the the ASS and the targets
acc=true(size(rCand,1),length(ASS));
for i=length(ASS):-1:1
    
    if ASS{i}.RsquareMax>data.RsquareLim
        
        if data.asdSetup.gradientModel==1 && data.asdSetup.doSubSampleR==false
            % we only have on dimension to use
            obs=ASS{i}.fitresult{3}((ASS{i}.W1(:,1)'*rCand')');
        else
            % we can use 2 subspace dimensions
            obs=ASS{i}.fitresult3D((ASS{i}.W1'*rCand')');
        end
        % decide if the full model should be run or not
        r=(obs-targets{i}.outer)/(targets{i}.target-targets{i}.outer);
        ixx= r<rand;
        acc(ixx,i)=false;
    end
end

% a parameter set is accepted if all observations are accepted
prBIN=all(acc,2);
% hence a "numeric value" is the sum over all observations
prX=sum(acc,2);



