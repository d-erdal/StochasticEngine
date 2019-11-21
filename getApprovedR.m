
% Function: Draws random parameter sets and checks whether they are
% acceptable according to the proxy model.
% This part handles drawing the sample and handling the archive.
% Part of: Stochastic sampling scheme version 4
%
% Created by: Daniel Erdal, Uni Tübingen, Nov. 2019
% Last edited by: N/A

% Input:
% ProyModel --> proxy model, in a cell structure with one model per observation
% ARCH  --> Archive of not-good-enough-but-still-pretty-reasonable parameter
%       sets that could be approved with an updated proxy model. Used on if
%       the corresponding flag is turned on
% targets   --> cell structure with the observation targets that the sampling
%       is based on. One cell per observation
% nrParam   --> number of parameters
% doArchive --> use the ARCH feature
% data      --> proxyModel.data, a free format stucture for prox-model inputs

% Output:
% newR  --> a row-vector of stage-1 approved parametes
% ARCH  --> updated archive (if active)
% rejectedR   --> matrix of rejected parameter-sets
% data      --> proxyModel.data, a free format stucture for prox-model inputs


function [newR,ARCH,rejectedR,data]=getApprovedR(ProxyModel,proxySetup,ARCH,targets,nrParam)

doArchive=proxySetup.doArchive;
data=proxySetup.data;

% minimum probability for a set of parameters
if proxySetup.doRandAccept && rand< proxySetup.randAcceptLevel
    % accept the first random parameter set
    rejectedR=[];
    newR=rand(1,nrParam);
else
    rejectedR=nan(1e4,nrParam+1);
    cc=1;
    while 1
        
        if  doArchive && ARCH.redo
            % new candidate is the old archive
            rCand=ARCH.data;
        else
            % draw 1,000 new Rs
            rCand=rand(1e3,nrParam);
        end
        
        % prBIN --> boolean vector for each parameters set if it passed the
        %   stage-1 test or not
        % prX --> numerical "probability" of each parameter set. No limits,
        %   but must be positive and higher must be better
      [prBIN,prX,data]=proxySetup.getApprovedR(rCand,ProxyModel,targets,proxySetup.data);
     
        % a parameter set is accepted if all observations are accepted
        if any(prBIN)
            % we have a new candidate
            ix=find(prBIN);
            newR=rCand(ix(1),:);
            
            if doArchive && ARCH.redo
                % just delete the candidate from the space
                ARCH.data(ix(1),:)=zeros(1,size(ARCH.data,2));
                ARCH.pr(ix(1),:)=0;
                if length(ix)==1
                    % if there are no more acceptable samples in the archive,
                    % skip testing again until after next training
                    ARCH.redo=false;
                end
            else
                rejectedR(cc:end,:)=[];
                if ix(1)>1
                   if isempty(rejectedR)
                        rejectedR=[rCand(1:ix(1)-1,:),prX(1:ix(1)-1)];                       
                   else
                        rejectedR=[rejectedR;[rCand(1:ix(1)-1,:),prX(1:ix(1)-1)]];
                   end
                end
            end
            break
            
        else
            if doArchive && ARCH.redo
                % there were no good samples in the archive
                ARCH.redo=false;
            else
                % store the rejected samples and the number of accepted
                % observations
                rejectedR(cc:cc+size(rCand,1)-1,:)=[rCand,prX];
                cc=cc+size(rCand,1);
            end
        end        
    end
    
    if doArchive && (~isempty(rejectedR) || all(isnan(rejectedR(:,1))))
        % Compare the stats of the archive with those in the rejected samples
        if min(ARCH.pr) < max(rejectedR(:,end))
            [~,ixAR]=sort(ARCH.pr);
            [~,ixRj]= sort(-rejectedR(:,end));
            for i=1:min(size(rejectedR,1),size(ARCH.pr,1))
                if rejectedR(ixRj(i),end)>ARCH.pr(ixAR(i))
                    % change them
                    ARCH.data(ixAR(i),:)=rejectedR(ixRj(i),1:end-1);
                    ARCH.pr(ixAR(i))=rejectedR(ixRj(i),end);
                else
                    % nothing more to change, break out of the loop
                    rejectedR(ixRj(1:i-1),:)=[]; % remove those moved to the ARCH
                    break
                end
            end
        end
    end
    
end
