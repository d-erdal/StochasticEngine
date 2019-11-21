% Master script for starting a stochastic run, locally or on a linux based
% cluster

% The scrit consists of three parts: setting up the stochastic engine,
% setting up the full model and setting up the proxy model.

% This is the template version

%% Stochastic Engine Setup

% Name to save the results under
stochSetup.saveName='templateTest';
% How many stage-1 acceptances between the saving?
stochSetup.saveEveryXruns=10;
% Sample everything from scratch?
stochSetup.doNewStart=true;
% or do a special start --> needs input into the main code!!!
stochSetup.doSpecialStart=false;
% Obs: if none of the two are active the currentStatus is attempted to be
% loaded

% number of runs to do before we start asking the proxy model for advice
stochSetup.nrFirstSamples=100;
% should we wait for all runs in the first sample to finish before
% submitting any new jobs?
stochSetup.waitForFirstSample=false;
% total number of stage-2 accepted samples we want before we are happy
stochSetup.nrFinalS2AccRuns=300;
% or, total number of finished samples (first of these two to be reached
% will cause the program to finish)
stochSetup.nrMaximumRuns=5000;

% Is the model very fast --> true (running locally and can be run sequentially
% without a cluster) or needs a cluster --> false
stochSetup.modelisImidiate=true;
% Number of parallel jobs to consider on a cluster (modelisImidiate=false)
stochSetup.currentRunsMaxSize=100;
% What time should we wait when the queue is full (modelisImidiate=false)
stochSetup.pauseInterval=3600/5;

% Should rejected runs from the proxy model be kept in memory and saved
stochSetup.storeRejected = false;
% When saving, should also the proxy model be saved?
% (for GPE models this can mean rather large files)
stochSetup.alsoSaveProxy = false;
% When saving, should also the modelSetup.data be saved?
% (if a GPE is the truth, this should be false)
stochSetup.alsoSaveModel = false;
% how often do we print info
stochSetup.printEveryXit=1;
% over what interval should short term acceptance be computed
% (for printing and Active learning)
stochSetup.shortTermACrange=100;

% There is also an option to run a non-standard function every time a new
%   parameter set is to be evaluated. 
% OBS: this is not a handle, but a string to be used with "eval"
% example: stochSetup.nonstandardFunction='[R,X]=messWithR(R,obs);'
% example of use: change the gradient model for ASS depending on the sample
%   size
% HINT: don't forget the ; inside the string!!
stochSetup.nonstandardFunction=[]; %

%% True Model Setup

% Function handle to running the model
% call:  [pid,obsX]=modelSetup.run(newR,c)
% where obsX is only used if stochSetup.modelisImidiate=true
% modelSetup.run = @runHGS_KBlite_V1;      
modelSetup.run = @runPROXYasTRUTH; % 

% Function handle to get the observations.
% Only if stochSetup.modelisImidiate=false;
% call: modelSetup.getObs [obs,crashed]=modelSetup.getObs(modelSetup,currentRuns(i,2));
% modelSetup.getObs = @getObsHGS_KBlite_V1;


% stucutre modelSetup.data is a free format that should include whatever
% the model might need as input.
modelSetup.data=[];
modelSetup.data.name ='GPE'; % AS or GPE for active subspace or GPE for proxy as truth
stupidCodeToGetSTKtoRun
% load FinalEndOfRun_GPEonTheFly GPE
modelSetup.data.GPE=GPE;% % for GPE as truth
% load EndOFRun_3000_ASS_Org_GPEdata fitresult3D W1
% modelSetup.data.fitresult3D=fitresult3D; % for AS as truth
% modelSetup.data.W1=W1; % for AS as truth

% On the cluster, what should be directories be named (will be this + a
% runing number)
% modelSetup.data.runDir='templateTest_run';
% modelSetup.data.HGSbaseFolder='version2BASE_StT1';

% number of model parmaeters to be randomized
modelSetup.nrParam=32;

% Number of observations and their targets
% Useage: name - not used; target - the traget to be considered by the
% sampler; targetTRUE - used for stage-2 acceptance (most of the time the
% same as target, but can be different if we want to give the samler some
% slack); outer - outer point for the ASSv1. OBS: also for the other
% sampling schemes this property is needed to determin if the aim is to be 
% above or below the target. This could be changed of course.
modelSetup.nrObs=6;
targets{1}.name='TopDrain';
targets{1}.target=-0.002; targets{1}.outer=-0.004; targets{1}.targetTRUE=-0.002;
targets{2,1}.name='GageC';
targets{2}.target=-0.005; targets{2}.outer=-0.003; targets{2}.targetTRUE=-0.005;
targets{3}.name='Stretch6';
targets{3}.target=-0.001; targets{3}.outer=-0.002; targets{3}.targetTRUE=-0.001;
targets{4}.name='Schoenbrunnen';
targets{4}.target=-0.000001; targets{4}.outer=-0.0000005; targets{4}.targetTRUE=-0.000005;
targets{5}.name='RiverRatioLower';
targets{5}.target=0.25; targets{5}.outer=0.2; targets{5}.targetTRUE=0.25;
targets{6}.name='RiverRatioUpper';
targets{6}.target=0.6; targets{6}.outer=0.75; targets{6}.targetTRUE=0.6;

for i=1:modelSetup.nrObs
    if targets{i}.targetTRUE>targets{i}.outer
        targets{i}.targetBELOW=false; targets{i}.mpl=1;
    else
        targets{i}.targetBELOW=true; targets{i}.mpl=-1;
    end
end
modelSetup.targets=targets; 
clear targets

% maximum wall-clock run time any model instance
modelSetup.maxRunTime=3600*20; % [s]


%% Proxy model setup
% Generally: the proxySetup.data is a stucture for all proxy-model specific
% inputs and will be passed back and forth with the training and samling

% at what interval should the training be done
% (stucture: first row values, second row corresponding ensemble sizes)
proxySetup.trainInterval=[100,200,300;1,1800,3000];

% Should the sampling be proceeded by an active learning session
proxySetup.doActiveLearning = false;
    % Function handle for the active learning call
    proxySetup.data.AL_run = @getApprovedR_activeLearning_GPE;   
    % at what interval should the training be done
    % (stucture: first row values, second row corresponding ensemble sizes)
    proxySetup.data.AL_estimateGPEparameterInterval=[10,200;1,100];%[1,50, 200 ; 1, 100, 1000];
    % What size should the MC-sample be (please keep to even 1e3s)
    proxySetup.data.AL_MCsampleSize=1e4;
    % maximum confusion for convergance (this rarely worked for me)
    proxySetup.data.AL_maxConfusion = 0.05;
    % maximum sample size before the AL is terminated
    proxySetup.data.AL_maxDataSet = 200;
    % should we exit the active learning when the short term acceptance
    % rate has reached a certain level? Should either be 1 (not active) or
    % around 0.5 when active (we should land here when things are stable)
    proxySetup.data.AL_maxShortTermAccRate=0.52;
    
% Function handle for getting a new stage-1 approved R
 % call: [newR,ARCH,rejR,proxySetup.data]=proxySetup.getApprovedR(proxySetup,ARCH,stochSetup.targets,modelSetup.nrParam)
%  proxySetup.getApprovedR = @getApprovedR_ASSv1;
%  proxySetup.getApprovedR = @getApprovedR_ASSv2;
 proxySetup.getApprovedR = @getApprovedR_GPE;
 
% Function handle for training the proxy model with on a data set
% call [proxyModel{i},proxySetup.data]=proxySetup.train(Rsave(1:c-1,:),obsMTRX(1:c-1,i),proxySetup,ActLearnResults);
%  proxySetup.train = @trainASS;
 proxySetup.train = @trainGPE; 
 
 % Accept a candidate point at a random ratio rand<proxySetup.randAcceptLevel
 % without asking the proxy model before
 proxySetup.doRandAccept = true;
 proxySetup.randAcceptLevel = 0.1;
 
 % Should the archiving feature be activeted (saves the best of the
 % rejeceted runs for re-testing after the proxy is updated)
 proxySetup.doArchive = true;
 proxySetup.maxARCHsize=10000;
 
% should the training be done in parallel
proxySetup.doParalleTraining=false;
proxySetup.parpoolSize=4;
 
%% for the GPE-proxy, use this section:===================================
% % % % % % initialize the GPE-code

proxySetup.data.probAccept=0.5; % accept a joint probability above this
      % % % % % % level. If random acceptance is wanted (>rand) set probAccept=nan;
proxySetup.data.covarianceFunction='stk_materncov32_aniso';
proxySetup.data.SIGMA2 = 1.0;  % variance parameter
proxySetup.data.NU     = 2.0;  % regularity parameter
proxySetup.data.RHO1   = 0.4*ones(modelSetup.nrParam,1);  % scale (range) parameter

% % % % at what interval should the parameter estimation be done
% % % % (stucture: first row values, second row corresponding ensemble sizes)
% % % % OBS: first row is "estimate every X training" (not sample size)
proxySetup.data.estimateGPEparameterInterval=[1,2,3,4,5;1,500,1000,1500,2000];

proxySetup.data.STKdir='stk2/';
addpath(proxySetup.data.STKdir);
stk_init;



%==========================================================================
%% for the Active Subspace Sampler V1, use this section:==================
% proxySetup.data.RsquareLim=0.7; % only sample if the r2 is higher than this limit
% % % % % AS-parameters
% proxySetup.data.asdSetup.gradientModel=3;
% proxySetup.data.asdSetup.doSubSampleR=false;
% proxySetup.data.asdSetup.nrSamplesToUse=NaN;
%==========================================================================

%% for the Active Subspace Sampler V2, use this section:==================
% proxySetup.data.acceptLimit=0.6; % only accept a sample if the mean of its neighbours are higher than
% proxySetup.data.nrNeighb=5;    % how many neighbors should minium be considered
% proxySetup.data.ASradius=0.01; % what is the search radius in the active suspace within which all neighbours are taken in
% % % % % % AS-parameters
% proxySetup.data.asdSetup.gradientModel=3;
% proxySetup.data.asdSetup.doSubSampleR=false;
% proxySetup.data.asdSetup.nrSamplesToUse=NaN;
%==========================================================================
%% small compatibility check

if isfield(proxySetup.data,'estimateGPEparameterInterval') && proxySetup.data.estimateGPEparameterInterval(2,1)~=1
    warning('proxySetup.data.estimateGPEparameterInterval(2,1) should be 1')
end
if proxySetup.trainInterval(2,1)~=1
    warning('proxySetup.trainInterval(2,1) should be 1')
end
if modelSetup.nrObs ~= length(modelSetup.targets)
    warning('modelSetup.nrObs ~= length(targets)')
end
if proxySetup.doParalleTraining && isempty(gcp('nocreate'))
    parpool(proxySetup.parpoolSize)
end

%% run the sampler

stochasticEngineMaster_version4_MASTER(stochSetup,modelSetup,proxySetup)






