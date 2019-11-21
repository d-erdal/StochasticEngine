% --In between sorting function--
% Function: Train an active subspace surrogate model for the input
% observations and parameters. 
% Part of: Stochastic sampling scheme version 4
%
% Created by: Daniel Erdal, Uni Tübingen, Nov. 2019
% Last edited by: N/A

% Input: 
% R     --> parameter matrix (scaled 0-1),
% obs 	--> observation vector, 
% data  --> proxySetup.data: proxy model setup structure
% ActLearnRes --> resulting R and obs from active learning (empty if not AL)
% iobs  --> index of the current observation

% Output: 
% GPE   --> the GPE-model for the current obs and R;
% data  --> the proxySetup.data structure including the stats and counters
%           to be passed back and forth between this and the main code



function [GPE,data]=trainGPE(R,obs,data,ActLearnRes,iobs)

if ~isempty(ActLearnRes)
    % if there has been an active learning phase, inlcude those results
    R=[ActLearnRes.Rfinal;R];
    obs=[ActLearnRes.obsFinal(:,iobs);obs];
end

% make sure to not take any crashed runs along in the training
R(isnan(obs),:)=[];
obs(isnan(obs),:)=[];

if ~isfield(data,'trainParamCounter') || length(data.trainParamCounter)<iobs
    data.trainParamCounter{iobs}=0;
end


% dimension
dim = size(R,2);

% Parameters are preset in the proxySetup-stucture
% SIGMA2 = 1.0;  % variance parameter
% NU     = 2.0;  % regularity parameter
% RHO1   = 0.4*ones(dim,1);  % scale (range) parameter
param0 = log ([data.SIGMA2; 1./data.RHO1]);

% initialte the GPE-model
GPE.model = stk_model (data.covarianceFunction, dim);

% model.lm = stk_lm_constant;
GPE.model.lm = stk_lm_affine;

tmp=(data.estimateGPEparameterInterval(1,data.estimateGPEparameterInterval(2,:)<=size(R,1)));
CURRENTInterval=tmp(end);
% estimate the parameters only on regular itervals
if ~isfield(data,'paramX') || length(data.paramX)<iobs ||  all(data.paramX{iobs}==0) || size(R,1)<100 || mod(data.trainParamCounter{iobs},CURRENTInterval)==0 
    % update the parameteres
   data.paramX{iobs}=stk_param_estim (GPE.model, R, obs, param0);
end
GPE.model.param = data.paramX{iobs};

% Setup the model
GPE.M_post = stk_model_gpposterior(GPE.model, R, obs);

% update the training counter
data.trainParamCounter{iobs}=data.trainParamCounter{iobs}+1;

% % % compute the kriging prediction
% % output_predict = stk_predict (model, input_train, output_train, Input);