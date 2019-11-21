% Function: Main function for the Stochastic Engine with purpose to sample
% behavioral parameters from a wast parameter space, utilizing a proxy
% model to apply 2-stage rejection sampling.
% 
% Handles all interactions with the proxy model as well as with the full
% model. If run on a linux cluster, this function also manages the active
% jobs.
% This is version 4 used in Erdal & Cirpka 2019, Erdal, Xiao, Nowak & 
% Cirpka 2020 and Erdal & Cirpka 2020 (as this is written 2019, the 2020
% publications are not yet accepted)
% Part of: Stochastic sampling scheme version 4
%
% Created by: Daniel Erdal, Uni Tübingen, Nov. 2019
% Last edited by: N/A
%
% For information about content of the input, please see
% Start_stochEng_v4_TEMPLATE.m
%
% Input:
% stochSetup    --> structure containing all information needed to run the
%           this main part of the  stochastic engine
% modelSetup    --> stucture containing all information needed to run the
%               full (heavy/slow) model
% proxySetup    --> stucture containing all information needed to setup and
%               run the proxy (surrogate/emulator) model
% 
% Output:
% this function has no output since all data is being saved at the end of
% the run

function stochasticEngineMaster_version4_MASTER(stochSetup,modelSetup,proxySetup)

%%
tidtagare=tic;
disp(datetime)
disp('===========================================================')
disp('Welcome to the stochastic engine. This is version 4.')
disp('This program was written by Daniel Erdal, Uni Tübingen, nov 2019')
disp('You are using the following setup: ')
disp(['Full model: ' func2str(modelSetup.run)])
disp(['Proxy model: ' func2str(proxySetup.getApprovedR)])
disp(['  with training: ' func2str(proxySetup.train)])
if proxySetup.doActiveLearning
    disp(['  with active learning : ' func2str(proxySetup.data.AL_run)])
end
if proxySetup.doArchive
    disp('  with achive function')
end
disp(['The sampler will run until ' num2str(stochSetup.nrFinalS2AccRuns) ' stage-2 accepted samples, or a maximum sample size of ' num2str(stochSetup.nrMaximumRuns)])
disp(['Results will be save under: ' stochSetup.saveName])
disp('===========================================================')
disp('Now we start')
disp(date)



if stochSetup.doNewStart
    % Initial sample drawn with a latin hypercube design
    Rinitial=lhsdesign(stochSetup.nrFirstSamples,modelSetup.nrParam);
    save([stochSetup.saveName 'Rinitial'],'Rinitial')
    % Matrix of current active (i.e. running) jobs
    currentRuns=nan(stochSetup.currentRunsMaxSize,2);
    % Vector of time taken to complete each run
    runTimes=zeros(stochSetup.nrFinalS2AccRuns*10,1);
    % Vector of current wallclock times of the current active jobs
    currentActiveTime=zeros(stochSetup.currentRunsMaxSize,1);
    currentTIC=cell(stochSetup.nrFinalS2AccRuns*100,1);
    
    % Matrix of observations to be saved
    obsMTRX=nan(stochSetup.nrFinalS2AccRuns*10,modelSetup.nrObs);
    % Matrix of parameter sets
    Rsave=nan(stochSetup.nrFinalS2AccRuns*10,modelSetup.nrParam);
    
    % Counters etc
    c=1; stats.crashedRuns=[]; stats.finishedRuns=[];
    Rcrash=[]; stats.nrFinishedRuns=0; stats.nrFinishedRunsLASTTraining=-inf;
    stats.nrAcceptedRuns=0; stats.rejectCounter=1; stats.rc=1;
    stats.trainParamCounter=0; stats.trainTime=[]; ARCH=[];
    proxySetup.doALFinished=false; stats.ix=0;
    rejectedR=zeros(0,modelSetup.nrParam+1); cSave=1;
    
    % initilize the archive
    if proxySetup.doArchive
        ARCH.data=zeros(proxySetup.maxARCHsize,modelSetup.nrParam);
        ARCH.pr=zeros(proxySetup.maxARCHsize,1);
        ARCH.redo = true;
    end
    
    ActLearnResults=[]; proxyModel=cell(modelSetup.nrObs,1);    
    
elseif stochSetup.doSpecialStart
    % here we put in how this is supposed to work --> MANUALLY
    error('Nothing here yet.....')
else
    disp('========OBS OBS OBS OBS ===================================')
    disp(['Starting from: currentStatus_' stochSetup.saveName])
    disp('========OBS OBS OBS OBS ===================================')
    load(['currentStatus_' stochSetup.saveName])
end

cOLD=c;  TX=1;
% Start the outer loop to create the behavioral parameter-sets
while ~isempty(currentRuns)
    
    if stochSetup.modelisImidiate || (isnan(currentRuns(end,1)) && ~(stochSetup.waitForFirstSample && c==stochSetup.nrFirstSamples+1 && stats.nrFinishedRuns<stochSetup.nrFirstSamples))
        % print somethings
        proxySetup.data.accRate=stats.nrAcceptedRuns/stats.nrFinishedRuns;
        sti=max(1,stats.nrFinishedRuns-stochSetup.shortTermACrange);
        proxySetup.data.STaccRate=sum(all(stats.ix(sti:end,:),2))/(stats.nrFinishedRuns-sti+1);      
        if mod(c,stochSetup.printEveryXit)==0
            disp('----------------------------------------')
            disp(['Iteration: ' num2str(c) '. Accepted runs: ' num2str(stats.nrAcceptedRuns) '.'])
            disp(['Acc. rate: ' num2str(round(proxySetup.data.accRate,2)) '. Short term acc. rate: ' num2str(round(proxySetup.data.STaccRate,2))])
        end
        %% run any non-standard function
        if isfield(stochSetup,'nonstandardFunction') && ~isempty(stochSetup.nonstandardFunction)
            eval(stochSetup.nonstandardFunction)
        end
        
        %% re-train the surrogate model
        
        % get the training iterval (stucture: first row values, second row
        % corresponding ensemble sizes)
        tmp=(proxySetup.trainInterval(1,proxySetup.trainInterval(2,:)<=c));
        proxySetup.CURRENTtrainInterval=tmp(end);
        
        % check if there is time for a training session
        if  stats.nrFinishedRuns>=stochSetup.nrFirstSamples && c-stats.nrFinishedRunsLASTTraining > proxySetup.CURRENTtrainInterval && ~proxySetup.doActiveLearning
            disp('-->  Re-training the proxy-model')
            ticTrain=tic;
            % create one proxy model for each observation
            if proxySetup.doParalleTraining && stats.trainParamCounter~=0
                parfor i=1:modelSetup.nrObs
                    [proxyModel{i},dataTMP{i}]=proxySetup.train(Rsave(1:c-1,:),obsMTRX(1:c-1,i),proxySetup.data,ActLearnResults,i);
                end              
                % go through each field in dataTMP and  update the 
                %    proxySetup.data. Especially, update those that are
                %    cells and nrObs long (i.e. the outputs from training)
                ff=fieldnames(dataTMP{1});
                for i=1:length(ff)
                    tmp=getfield(dataTMP{1},ff{i});
                    if iscell(tmp) && length(tmp)==modelSetup.nrObs
                        for j=1:length(tmp)
                            tmp2=getfield(dataTMP{j},ff{i});
                            tmp{j}=tmp2{j};
                        end
                    end
                    setfield(proxySetup.data,ff{i},tmp);
                    
                end
                
                
            else
                for i=1:modelSetup.nrObs
                    [proxyModel{i},proxySetup.data]=proxySetup.train(Rsave(1:c-1,:),obsMTRX(1:c-1,i),proxySetup.data,ActLearnResults,i);
                end
            end
            % update the stats
            stats.trainTime=[stats.trainTime;toc(ticTrain)];
            stats.nrFinishedRunsLASTTraining=c;
            stats.trainParamCounter=stats.trainParamCounter+1;
            
            % set the archive recomputation flag to true
            ARCH.redo = true;
            
            
        end
        
        %%
        % we have a free slot, submit a new one, if we have not yet aquired
        % enough stage-2 accepted samples or crossed the maximum number of
        % allowed simulations
        if (stats.nrAcceptedRuns <= stochSetup.nrFinalS2AccRuns && stats.nrFinishedRuns<=stochSetup.nrMaximumRuns) || proxySetup.doActiveLearning
            
            % if we use the second version of the active subspace sampler
            % we need some additional input
            if strcmp(func2str(proxySetup.getApprovedR),'getApprovedR_ASSv2') && c>stochSetup.nrFirstSamples
                proxySetup.data.Rsave=Rsave;
                proxySetup.data.accIndex=stats.ix;
            end                      
            
            %% draw a new R that is accepted by the Proxy model
            
            if stats.nrFinishedRuns<=stochSetup.nrFirstSamples
                % sample comes from the initial latin hypercube
                if c <= size(Rinitial,1)
                    newR=Rinitial(c,:);
                else
                    % we are just waiting for the first sample here
                    newR=rand(1,modelSetup.nrParam);
                end
            elseif proxySetup.doActiveLearning
                % get the R by asking where we would like to sample to
                % reduce the uncertainty of our Proxy model the most
                [newR,obsExp,proxySetup.doALFinished,proxySetup.data]=proxySetup.data.AL_run(Rsave,obsMTRX,modelSetup.targets,modelSetup.nrParam,proxySetup.data,proxySetup.doParalleTraining);
                % and add the expected observation
                obsMTRX(c,:)=obsExp;
            else
                % Get a new parameter-set by finding a stage-1 accepted
                % parameter set given the current state of the model
                [newR,ARCH,rejR,proxySetup.data]=getApprovedR(proxyModel,proxySetup,ARCH,modelSetup.targets,modelSetup.nrParam);
                
                if stochSetup.storeRejected && ~isempty(rejR)
                    rejectedR=[rejectedR;rejR];
                end
                stats.rejectCounter=stats.rejectCounter+size(rejR,1);
            end
            
            if strcmp(func2str(proxySetup.getApprovedR),'getApprovedR_ASSv2')
                % to keep saving volumes down, delete this matrix once done
                proxySetup.data.Rsave=[];
                proxySetup.data.accIndex=[];
            end
            
            
            if size(rejectedR,1)>1e4
                save(['rejectedR_' num2str(stats.rc)],'rejectedR')
                stats.rc=stats.rc+1;
                rejectedR=zeros(0,size(rejectedR,2));
            end
            
            if any(sum(abs(Rsave-newR),2)==0)
                % we have a douplicate
                save DUPLICATEERROR
                error('newR is not unique....')
            end
            
            % save in Rsave
            Rsave(c,:)=newR;
            %%
            % start the real model with the new R
            [pid,obsX]=modelSetup.run(newR,modelSetup.data,c);
            currentRuns(end,:)=[pid,c];
            
            if stochSetup.modelisImidiate
                % if the model is so fast we do not submit anything
                obsMTRX(c,:)=obsX;
                stats.finishedRuns=[stats.finishedRuns;c];
            else
                disp(datetime)
                disp(['Submitted run ' num2str(c)])
            end
            c=c+1;
        else
            disp('All runs submitted. Removing 1 slot in the current runs')
            currentRuns(end,:)=[];
            stochSetup.currentRunsMaxSize=stochSetup.currentRunsMaxSize-1;
            disp(['  --> Current job size: ' num2str(stochSetup.currentRunsMaxSize)])
        end
        
        
    else
        % No free slots to be submitted: should we wait?
        if c~=cOLD; fprintf(1,'      Now we wait'); cOLD=c; end
        fprintf(1,'.');
        
        if TX>5; pause(stochSetup.pauseInterval); else; TX=TX+1; end    
        
        % check if any of the current runs are finished or have run too long
        for i=1:size(currentRuns,1)
            
            if ~isnan(currentRuns(i,2)) % i.e. is the job active
                
                % is the job done?
                if exist([modelSetup.data.runDir num2str(currentRuns(i,2)) '/done.txt'],'file')
                    % the run is finished
                    fprintf(1,'\n')
                    disp(['Run ' num2str(currentRuns(i,2)) ' is finished!'])
                    
                    try
                        % check the time (a bit unstable, thats why the try)
                        a=dir([modelSetup.data.runDir num2str(currentRuns(i,2)) '/start.txt']);
                        b=dir([modelSetup.data.runDir num2str(currentRuns(i,2)) '/done.txt']);
                        runTimes(currentRuns(i,2),1)=hours(datetime(b.date)-datetime(a.date));
                    end
                    
                    % add it to the stats.finishedRuns list
                    stats.finishedRuns=[stats.finishedRuns;currentRuns(i,2)];
                    stats.nrFinishedRuns=stats.nrFinishedRuns+1;
                    
                    % read the observations
                    [obs,crashed]=modelSetup.getObs(modelSetup,currentRuns(i,2));
                    
                    if crashed
                        % the job produced to output and should be removed
                        stats.crashedRuns=[stats.crashedRuns;currentRuns(i,2)];
                        Rcrash=[Rcrash;Rsave(currentRuns(i,2),:)];
                        Rsave(currentRuns(i,2),:)=nan;
                        obsMTRX(currentRuns(i,2),:)=nan;
                        disp(['OBS: Run ' num2str(currentRuns(i,2)) ' is a crash!'])
                    else
                        obsMTRX(currentRuns(i,2),:)=obs;
                    end
                    
                    % reset the current rusn to allow for new runs
                    currentRuns(i,:)=nan;
                    
                elseif isempty(currentTIC{currentRuns(i,2)})
                    % run is submitted but maybe not yet started
                    if exist([modelSetup.data.runDir num2str(currentRuns(i,2)) '/start.txt'],'file')
                        currentTIC{currentRuns(i,2)}.now=now;
                    end
                    
                elseif (now-currentTIC{currentRuns(i,2)}.now)*3600*24>modelSetup.maxRunTime
                    % the run has gone too long, kill it
                    system(['qdel ' num2str(currentRuns(i,1))]);
                    system(['rm ' modelSetup.data.runDir num2str(currentRuns(i,2)) '/scratch_debug']);
                    system(['rm ' modelSetup.data.runDir num2str(currentRuns(i,2)) '/KaesBachLite_v1o.observation_well_flow.ObsWell_nr*']);
                    system(['echo crash >' modelSetup.data.runDir num2str(currentRuns(i,2)) '/CRASH.txt'])
                    
                    stats.crashedRuns=[stats.crashedRuns;currentRuns(i,2)];
                    Rcrash=[Rcrash;Rsave(currentRuns(i,2),:)];
                    Rsave(currentRuns(i,2),:)=nan;
                    
                    % add it to the stats.finishedRuns list
                    stats.finishedRuns=[stats.finishedRuns;currentRuns(i,2)];
                    stats.nrFinishedRuns=stats.nrFinishedRuns+1;
                    
                    fprintf(1,'\n')
                    disp(['Run ' num2str(currentRuns(i,2)) ' has been removed.'])
                    currentRuns(i,:)=nan;
                    
                else
                    % just keep track of the rough sim-times
                    currentActiveTime(i)=(now-currentTIC{currentRuns(i,2)}.now)*24;
                end
            end
        end
    end
    
    %% is the active learning now finished?
    if proxySetup.doActiveLearning && proxySetup.doALFinished
        disp('========Active learning now completed===========')
        proxySetup.doActiveLearning=false;
%         Rsave(isnan(obsMTRX(:,1)),:)=[];
%         obsMTRX(isnan(obsMTRX(:,1)),:)=[];
%         
        ActLearnResults.RunsTotal=c;
        ActLearnResults.stats=stats;
        
        ActLearnResults.Rfinal=Rsave(stats.finishedRuns,:);
        ActLearnResults.obsFinal=obsMTRX(stats.finishedRuns,:);
        
        % save the current state of the system
        save([stochSetup.saveName 'ActLearnFinished'])
        
        if ~stochSetup.modelisImidiate
            % now kill all current runs and restart the sampler from scratch
            for ikill=1:size(currentRuns,1)
                system(['qdel ' num2str(currentRuns(ikill,1))]);
            end
            currentRuns=nan(size(currentRuns));
            % rename all the folders
            for i=1:c
                cm=['mv ' modelSetup.data.runDir num2str(i) ' ' modelSetup.data.runDir num2str(i) '_ActLearn'];
                system(cm);
            end
        end
        
        stats.crashedRuns=[]; stats.finishedRuns=[]; runTimes=zeros(stochSetup.nrFinalS2AccRuns*10,1);
        currentActiveTime=zeros(stochSetup.currentRunsMaxSize,1);
        currentTIC=cell(stochSetup.nrFinalS2AccRuns*100,1);
        obsMTRX=nan(stochSetup.nrFinalS2AccRuns*10,modelSetup.nrObs);
        Rsave=nan(stochSetup.nrFinalS2AccRuns*10,modelSetup.nrParam);
        c=1;
        Rcrash=[]; stats.nrFinishedRuns=0; stats.nrFinishedRunsLASTTraining=-inf; stats.nrAcceptedRuns=0;
        rejectedR=[]; stats.rejectCounter=1; stats.rc=1;
        stochSetup.nrFirstSamples=-1; % 0 does not work!!!
    end
    %%
    
    if ~stochSetup.modelisImidiate
        % resort the currentRuns matrix so that the nans comes at the end
        [~,ix]=sort(currentRuns(:,2));
        currentRuns=currentRuns(ix,:);
        
        % read the optionalCommandsFile.txt
        % this could be used to temper with the parameters, or to save. Most
        % commonly it is used to alter stochSetup.currentRunsMaxSize to alter
        % the number of jobs occupying a cluster
        try eval(textread('optionalCommandsFile.txt','%c')); end
        system('touch empty.txt');
        copyfile('empty.txt','optionalCommandsFile.txt');
        
        while isnan(currentRuns(end,1)) && size(currentRuns,1)>stochSetup.currentRunsMaxSize
            % we have too large of a system, remove
            currentRuns(end,:)=[];
        end
        
        if size(currentRuns,1)<stochSetup.currentRunsMaxSize
            currentRuns(end+1:stochSetup.currentRunsMaxSize,:)=nan;
        end
    else
        stats.nrFinishedRuns=c-1;
    end
    
    % check status of accepted runs
    % (not just from the proxy, but hard checks against the full model)
    stats.ix=false(stats.nrFinishedRuns,modelSetup.nrObs);
    for i=1:modelSetup.nrObs
        stats.ix((obsMTRX(1:stats.nrFinishedRuns,i)*modelSetup.targets{i}.mpl)>(modelSetup.targets{i}.targetTRUE*modelSetup.targets{i}.mpl),i)=true;
    end
    stats.ix(stats.crashedRuns)=false;
    stats.nrAcceptedRuns=sum(all(stats.ix,2));
    
    % save while running
    dataX=proxySetup.data;
    save(['currentRnObs_' stochSetup.saveName '_' num2str(round(now))],'Rsave','obsMTRX','stats','dataX')
    if mod(cSave,stochSetup.saveEveryXruns)==0 || c < 100
        if ~stochSetup.alsoSaveModel
            % temporaily clear the data to reduce saving size
            dataM_TEMP=modelSetup.data;
            modelSetup.data=[];
        end
        saveList=who; % list of all variables in the workspace
       ixNoSave=[];
        for iii=1:length(saveList)            
            if (~stochSetup.alsoSaveProxy && strcmp(saveList{iii},'proxyModel')) || strcmp(saveList{iii},'dataM_TEMP') || strcmp(saveList{iii},'tmp')
                % aremove them from the list
                ixNoSave=[ixNoSave;iii];
            end
        end
        saveList(ixNoSave)=[];
        
        % and do the saving
        save(['currentStatus_' stochSetup.saveName],saveList{:})
        clear iii ixNoSave saveList 
        
        if ~stochSetup.alsoSaveModel
            % reset the data structure
            modelSetup.data=dataM_TEMP;
            clear dataM_TEMP
        end
    end
    cSave=cSave+1;
    
end

totalTime=toc(tidtagare);
save(['EndOfRun_' stochSetup.saveName])

disp('=========================DONE================================')

disp(['Final results are saved in: EndOfRun_' stochSetup.saveName])
disp(['and located at ' pwd])
disp(['Total wall-clock run time was: ' num2str(round(totalTime/3600,3)) ' hrs'])























